# -*- coding: UTF-8 -*-
""" 
``kraken2``
--------------------------------

:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

.. Note:: This module was developed as part of a study led by Dr. Jacob Moran Gilad

A module for running ``kraken`` version 2:

Pass the full path to the ``kraken2`` executable in ``script_path``.

Merging of sample kraken reports in done with krona. See the section on Parameters that can be set.

You can follow this module with the ``kraken-biom`` module to create a biom table from the reports.

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* fastq files, either paired end or single:

    * ``sample_data[<sample>]["fastq.F"]``
    * ``sample_data[<sample>]["fastq.R"]``
    * ``sample_data[<sample>]["fastq.S"]``
    
    

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Puts the ``kraken`` output files in:  

    * ``self.sample_data[<sample>]["raw_classification"]``
    * ``self.sample_data[<sample>]["classification"]``
    * ``self.sample_data[<sample>]["kraken.report"]``
    
    * If ``ktImportTaxonomy_path`` parameter was passed, puts the krona reports in 

    * ``self.sample_data["project_data"]["krona"]``


    
Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "ktImportTaxonomy_path", "", "Path to ktImportTaxonomy. You can additional ``ktImportTaxonomy`` parameters at the end of the path. If not passed, the ``krona`` report will not be built."

    
Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    kraken1:
        module: kraken
        base: trim1
        script_path: {Vars.paths.kraken}
        qsub_params:
            -pe: shared 20
        ktImportTaxonomy_path: /path/to/ktImportTaxonomy  -u  http://krona.sourceforge.net
        redirects:
            --db: /path/to/kraken_std_db
            --preload: 
            --quick: 
            --threads: 20

References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Wood, D.E. and Salzberg, S.L., 2014. **Kraken: ultrafast metagenomic sequence classification using exact alignments**. *Genome biology*, 15(3), p.R46.

"""


import os
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept

from pkg_resources import resource_filename


__author__ = "Menachem Sklarz"
__version__ = "1.6.0"


class Step_kraken2(Step):
    """ A class that defines a pipeline step name (=instance).
    """
    
    def step_specific_init(self):
        """ Called on intiation
            Good place for parameter testing.
            Wrong place for sample data testing
        """
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        self.file_tag = ".kraken.out"

            
        # Checking this once and then applying to each sample:
        if "--db" not in list(self.params["redir_params"].keys()):
            raise AssertionExcept("--db not set.\n")

        if "scope" in self.params:
            if self.params["scope"] == "project":
                pass
            elif self.params["scope"] == "sample":
                for sample in self.sample_data["samples"]:  # Getting list of samples out of samples_hash
                    pass
            else:
                raise AssertionExcept("'scope' must be either 'sample' or 'project'")
        else:
            self.write_warning("No 'scope' specified. Using 'sample' scope")
            self.params["scope"] = "sample"

        # For backwards comaptibility:
        if "ktImportTaxonomy_path" in list(self.params):
            self.params["ktImportTaxonomy"] = dict()
            self.params["ktImportTaxonomy"]["path"] = self.params["ktImportTaxonomy_path"]
            self.params["ktImportTaxonomy"]["redirects"] = ""

    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        pass
        
    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        
        self.make_sample_file_index()   # see definition below
        
        # Get full path to kraken module:
        # kraken_path = os.path.dirname(os.path.realpath(__file__))
        merge_kraken_reports = resource_filename(__name__, 'merge_kraken_reports.R')

        self.script = self.get_setenv_part()

        ### Add code to create a unified krona plot

        if "ktImportTaxonomy" in list(self.params):
            if "redirects" in self.params["ktImportTaxonomy"]:
                redirects = " \\\n\t" + " \\\n\t".join([key+" "+val for key,val in self.params["ktImportTaxonomy"]["redirects"].items()])
            else:
                redirects = ""

            self.script += """# Running ktImportTaxonomy to create a krona chart for samples
{ktImportTaxonomy_path} \\
\t-o {out} \\
\t-q 1 \\
\t-t 2 \\
\t""".format(ktImportTaxonomy_path=self.params["ktImportTaxonomy"]["path"]+" "+redirects,
                      out=self.base_dir + self.sample_data["Title"] + "_krona_report.html")
            
            for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
                self.script += "%s.forKrona,%s \\\n\t" % (self.sample_data[sample]["raw_classification"],sample)
            # Remove final \\\n\t
            self.script = self.script.rstrip("\\\n\t")
            self.script += "\n\n"

            self.sample_data["project_data"]["krona"] = self.base_dir + "krona_report.html"
    
            self.stamp_file(self.sample_data["project_data"]["krona"])

    def build_scripts(self):
        """ This is the actual script building function
            Most, if not all, editing should be done here 
            HOWEVER, DON'T FORGET TO CHANGE THE CLASS NAME AND THE FILENAME!
        """

        if self.params["scope"] == "project":
            sample_list = ["project_data"]
        else:   #self.params["scope"] == "sample"
            sample_list = self.sample_data["samples"]

        for sample in sample_list:  # Getting list of samples out of samples_hash

            # Name of specific script:
            self.spec_script_name = self.set_spec_script_name(sample)
            self.script = ""

            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)

            # Define output filename 
            output_filename = ".".join([sample , self.file_tag])

            ######### Step 1, run kraken itself
            # Create script and write to SCRPT
            # Add parameters passed to main script by user:
            # If the following params are not supplied by the user, add the defaults...

            if "fastq.F" in self.sample_data[sample] and "fastq.R" in self.sample_data[sample]:
                reads = """\
--paired \\
\t{forward} \\
\t{reverse} """.format(forward=self.sample_data[sample]["fastq.F"],
                         reverse=self.sample_data[sample]["fastq.R"])
            elif "fastq.S" in self.sample_data[sample]:
                reads = self.sample_data[sample]["fastq.S"]
            else:
                self.write_warning("KRAKEN on mixed PE/SE samples is not defined. Using only PE data!\n")

            self.script = """
{const}--output {out} \\
\t--report {out}.report \\
\t--unclassified-out {out}.unclassified#.fq \\
\t--classified-out {out}.classified#.fq \\
\t{reads}
            """.format(out=use_dir+output_filename,
                       const=self.get_script_const(),
                       reads=reads)

            ######### Step 4, create krona report:
            if "ktImportTaxonomy" in list(self.params.keys()):
                self.script += """
# Create file for ktImportTaxonomy
if [ -e {krak_out} ]
then
    cut -f 2,3 {krak_out} \\
        > {krak_out}.forKrona
fi

""".format(krak_out=use_dir+output_filename)
    
            # Storing the output file in $samples_hash
            self.sample_data[sample]["raw_classification"] = "%s" % (sample_dir + output_filename)
            self.sample_data[sample]["unclassified"] = "%s.unclassified" % (sample_dir + output_filename)
            self.sample_data[sample]["classified"] = "%s.classified" % (sample_dir + output_filename)
            self.sample_data[sample]["kraken.report"] = "%s.report" % (sample_dir + output_filename)
                
            self.stamp_file(self.sample_data[sample]["raw_classification"])
            self.stamp_file(self.sample_data[sample]["unclassified"])
            self.stamp_file(self.sample_data[sample]["classified"])
            self.stamp_file(self.sample_data[sample]["kraken.report"])
            
            # Move all files from temporary local dir to permanent base_dir
            self.local_finish(use_dir,self.base_dir)
            self.create_low_level_script()

    def make_sample_file_index(self):
        """ Make file containing samples and target file names for use by kraken analysis R script
        """
        
        with open(self.base_dir + "kraken_files_index.txt", "w") as index_fh:
            index_fh.write("Sample\tkraken_report\n")
            for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
                index_fh.write("%s\t%s\n" % (sample,self.sample_data[sample]["kraken.report"]))
                
        self.sample_data["project_data"]["kraken_file_index"] = self.base_dir + "kraken_files_index.txt"
