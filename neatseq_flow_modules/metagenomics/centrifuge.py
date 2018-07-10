# -*- coding: UTF-8 -*-
""" 
``centrifuge``
--------------------------------

:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

.. Note:: This module was developed as part of a study led by Dr. Jacob Moran Gilad

A module for running ``centrifuge``:

Pass the full path to the ``centrifuge`` executable in ``script_path``.

Merging of sample centrifuge reports in done with krona. See the section on Parameters that can be set.

.. CHECK THIS WORKS!! You can follow this module with the ``kraken-biom`` module to create a biom table from the reports.

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* fastq files, either paired end or single:

    * ``sample_data[<sample>]["fastq.F"]``
    * ``sample_data[<sample>]["fastq.R"]``
    * ``sample_data[<sample>]["fastq.S"]``
    
    

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Puts the ``centrifuge`` output files in:  

    * ``self.sample_data[<sample>]["raw_classification"]``
    * ``self.sample_data[<sample>]["classification"]``
    * ``self.sample_data[<sample>]["classification_report"]``
    
* If ``ktImportTaxonomy_path`` parameter was passed, puts the krona reports in 

    * ``self.sample_data["krona"]``


    
Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "ktImportTaxonomy_path", "", "Path to ktImportTaxonomy. You can additional ``ktImportTaxonomy`` parameters at the end of the path. If not passed, the ``krona`` report will not be built."

    
Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    Centrifuge:
        module:         centrifuge
        base:           trim1
        script_path:    {Vars.paths.centrifuge}
        qsub_params:
            -pe:        shared 20
        ktImportTaxonomy_path: /path/to/ktImportTaxonomy  -u  http://krona.sourceforge.net
        redirects:
            --db:       /path/to/centrifuge_db
            --preload: 
            --quick: 
            --threads:  20

References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Kim, D., Song, L., Breitwieser, F. P., & Salzberg, S. L. (2016). **Centrifuge: rapid and sensitive classification of metagenomic sequences**. *Genome research*, 26(12), 1721-1729.


"""


import os, sys, re
from neatseq_flow.PLC_step import Step,AssertionExcept

from pkg_resources import resource_filename


__author__ = "Menachem Sklarz"
__version__ = "1.2.0"


class Step_centrifuge(Step):
    """ A class that defines a pipeline step name (=instance).
    """
    
    def step_specific_init(self):
        """ Called on intiation
            Good place for parameter testing.
            Wrong place for sample data testing
        """
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        self.file_tag = ".centrifuge.out"

            
        # Checking this once and then applying to each sample:
        if "-x" not in self.params["redir_params"].keys():
            raise AssertionExcept("You didn't pass a database with -x in redirects.\n")

            
        
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
        merge_centrifuge_reports = resource_filename(__name__, 'merge_centrifuge_reports.R')

        ### Add code to create a unified krona plot
        if "ktImportTaxonomy_path" in self.params.keys():
            
            sample_files = ""
            for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
                sample_files += "%s.forKrona,%s \\\n\t" % (self.sample_data[sample]["raw_classification"],sample)
            self.script = """
# Running ktImportTaxonomy to create a krona chart for samples
{ktImportTaxonomy_path} \\
    -o {output_fn}%s \\
    -q 1 \\
    -t 2 \\
    {sample_files}
""".format(ktImportTaxonomy_path=self.params["ktImportTaxonomy_path"],
           output_fn = (self.base_dir + self.sample_data["Title"] + "_krona_report.html"),
           sample_files = sample_files)
            
            self.script = re.sub("\\\\\s*$","\n\n",self.script)

            
            self.sample_data["krona"] = self.base_dir + "krona_report.html"
    
            self.stamp_file(self.sample_data["krona"])
            

        
    
    def build_scripts(self):
        """ This is the actual script building function
            Most, if not all, editing should be done here 
            HOWEVER, DON'T FORGET TO CHANGE THE CLASS NAME AND THE FILENAME!
        """
        
    
        
        # Each iteration must define the following class variables:
            # spec_script_name
            # script
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
            
            # Name of specific script:
            self.spec_script_name = self.set_spec_script_name(sample)
            self.script = ""

            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)
                
 
         
            # Define output filename 
            output_filename = "".join([use_dir , sample , self.file_tag])
        
        
            ######### Step 1, run centrifuge itself
            # Create script and write to SCRPT
            # Add parameters passed to main script by user:
            # If the following params are not supplied by the user, add the defaults...

            self.script += self.get_script_const()
            # self.script += "%s \\\n\t" % self.params["script_path"]
            # self.script += self.get_redir_parameters_script()
            self.script += "-S %s \\\n\t" % output_filename

            if "fastq.F" in self.sample_data[sample]:
                self.script += "-1 {m1} \\\n\t-2 {m2} \\\n\t".format(
                                    m1=self.sample_data[sample]["fastq.F"],
                                    m2=self.sample_data[sample]["fastq.R"])
            if "fastq.S" in self.sample_data[sample]:
                self.script += "-U {r} \n\n".format(r=self.sample_data[sample]["fastq.S"])
            self.script = re.sub("\\\\\n$","\n\n",self.script)
                
            # Find path to centrifuge scripts
            centrifuge_path = os.path.dirname(self.params["script_path"])
            if centrifuge_path:
                centrifuge_path = centrifuge_path + os.sep

            ######### Step 2, translate raw centrifuge into useful names
            kreport_params = ((self.params["centrifuge-kreport"] + " \\\n\t") if "centrifuge-kreport" in self.params else "") 
            self.script += """
# Create useful centrifuge output 
if [ -e {centrifuge_out} ]
then
{centrifuge_path}centrifuge-kreport \\
    -x {db} \\
    {params} {centrifuge_out} > \\
    {centrifuge_out}.report
fi
""".format(centrifuge_path = centrifuge_path, 
            centrifuge_out = output_filename,
            params        = kreport_params,
            db            = self.params["redir_params"]["-x"])

            ######### Step 3, create krona report:
            if "ktImportTaxonomy_path" in self.params.keys():
                self.script += """
# Create file for ktImportTaxonomy
if [ -e {centrifuge_out} ]
then
    cut -f 1,3 {centrifuge_out} \\
        > {centrifuge_out}.forKrona
fi

""".format(centrifuge_out = output_filename)
    
            # Storing the output file in $samples_hash
            self.sample_data[sample]["raw_classification"]        = "%s" % (sample_dir + os.path.basename(output_filename))
            self.sample_data[sample]["classification"]        = "%s.labels" % (sample_dir + os.path.basename(output_filename))
            self.sample_data[sample]["classification_report"] = "%s.report" % (sample_dir + os.path.basename(output_filename))
                
            self.stamp_file(self.sample_data[sample]["raw_classification"])
            self.stamp_file(self.sample_data[sample]["classification"])
            self.stamp_file(self.sample_data[sample]["classification_report"])
            
            # Move all files from temporary local dir to permanent base_dir
            self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
                        
            
            self.create_low_level_script()
                    

                    
    def make_sample_file_index(self):
        """ Make file containing samples and target file names for use by centrifuge analysis R script
        """
        
        with open(self.base_dir + "centrifuge_files_index.txt", "w") as index_fh:
            index_fh.write("Sample\tcentrifuge_report\n")
            for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
                index_fh.write("%s\t%s\n" % (sample,self.sample_data[sample]["classification_report"]))
                
        self.sample_data["centrifuge_file_index"] = self.base_dir + "centrifuge_files_index.txt"
        