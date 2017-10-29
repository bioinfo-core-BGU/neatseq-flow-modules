
""" 
Module ``kraken``
--------------------------------

A module for running ``kraken``:

Note that ``kraken`` executable must be in a folder together with ``kraken-translate`` and ``kraken-report``. This is the default for ``kraken`` installation. 

Pass the full path to the ``kraken`` executable in ``script_path``.

Merging of sample kraken reports in done in two ways: with an R script shipped with the module, and with krona. See the section on Parameters that can be set.

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
    * ``self.sample_data[<sample>]["classification_report"]``
    
* If ``ktImportTaxonomy_path`` parameter was passed, puts the krona reports in 

    * ``self.sample_data["krona"]``

* If ``merge_kraken_reports`` was passed, puts the merged reports in 

    * ``self.sample_data["biom_table"]``
    * ``self.sample_data["biom_table_summary"]``
    * ``self.sample_data["biom_table_tsv"]``


    
Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "merge_kraken_reports", "", "If set, will merge the reports with ``merge_kraken_reports.R``, which is shipped with the module. If ``Rscript`` is not in the path, you will have to use ``setenv`` to add it to the path."
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
        merge_kraken_reports:
        ktImportTaxonomy_path: /path/to/ktImportTaxonomy  -u  http://krona.sourceforge.net
        redirects:
            --db: /path/to/kraken_std_db
            --preload: 
            --quick: 
            --threads: 20

            
"""


import os
import sys
from PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "0.2.0"
class Step_kraken(Step):
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
        if "--db" not in self.params["redir_params"].keys():
            raise AssertionExcept("--db not set.\n")

            
        
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        pass
        
    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        
        self.make_sample_file_index()   # see definition below
        
        # Get full path to kraken module:
        kraken_path = os.path.dirname(os.path.realpath(__file__))
        try:
            self.params["merge_kraken_reports"]
        except KeyError:
            self.write_warning("You asked not to merge kraken reports...\n")
            self.script = ""
        else:
            self.script = "# Creating biom files from all kraken reports\n"
            self.script += "Rscript %s \\\n\t" % (kraken_path + os.sep + "merge_kraken_reports.R")  # self.params["kraken_merge_path"]
            self.script += "-f %s \\\n\t" % self.sample_data["kraken_file_index"]
            self.script += "--db %s \\\n\t" % self.params["redir_params"]["--db"]
            self.script += "-o %s\n\n" % self.base_dir

            self.sample_data["biom_table_tsv"] = "".join([self.base_dir,"kraken_biom.count.sum.tsv"])
            self.sample_data["biom_table_summary"] = "".join([self.base_dir,"kraken_biom.count.sum.biom.summary"])
            self.sample_data["biom_table"] = "".join([self.base_dir,"kraken_biom.count.sum.biom"])
        
            self.stamp_file(self.sample_data["biom_table_tsv"])
            self.stamp_file(self.sample_data["biom_table_summary"])
            self.stamp_file(self.sample_data["biom_table"])

        # # Create a qiime slot to store the biom table in case you want downstream analysis with QIIME
        # if "qiime" not in self.sample_data.keys():
            # self.sample_data["qiime"] = dict()
            
        ### Add code to create a unified krona plot
        if "ktImportTaxonomy_path" in self.params.keys():
            self.script += "# Running ktImportTaxonomy to create a krona chart for samples\n"
            self.script += "%s \\\n\t" % self.params["ktImportTaxonomy_path"]
            self.script += "-o %s \\\n\t" % (self.base_dir + self.sample_data["Title"] + "_krona_report.html")
            self.script += "-q 1 \\\n\t"
            self.script += "-t 2 \\\n\t"
            
            for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
                self.script += "%s.forKrona,%s \\\n\t" % (self.sample_data[sample]["raw_classification"],sample)
            # Remove final \\\n\t
            self.script = self.script.rstrip("\\\n\t")
            self.script += "\n\n"

            self.sample_data["krona"] = self.base_dir + "krona_report.html"
    
            self.stamp_file(self.sample_data["krona"])
            

        
    
    def build_scripts(self):
        """ This is the actual script building function
            Most, if not all, editing should be done here 
            HOWEVER, DON'T FORGET TO CHANGE THE CLASS NAME AND THE FILENAME!
        """
        
    
        
        if "--preload" not in self.params["redir_params"].keys():
            self.write_warning("Not setting --preload, but it IS recommended...\n")

        # Each iteration must define the following class variables:
            # spec_script_name
            # script
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
            
            # Name of specific script:
            self.spec_script_name = "_".join([self.step,self.name,sample])
            self.script = ""

            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)
                
 
         
            # Define output filename 
            output_filename = "".join([use_dir , sample , self.file_tag])
        
        
            ######### Step 1, run kraken itself
            # Create script and write to SCRPT
            # Add parameters passed to main script by user:
            # If the following params are not supplied by the user, add the defaults...

            self.script += self.get_script_const()
            # self.script += "%s \\\n\t" % self.params["script_path"]
            # self.script += self.get_redir_parameters_script()
            self.script += "--output %s \\\n\t" % output_filename

            if "PE" in self.sample_data[sample]["type"]:
                self.script += "--paired \\\n\t%s \\\n\t%s \n\n" % (self.sample_data[sample]["fastq.F"],self.sample_data[sample]["fastq.R"])
            elif "SE" in self.sample_data[sample]["type"]:
                self.script += "%s \n\n" % self.sample_data[sample]["fastq.S"]
            else:
                self.write_warning("KRAKEN on mixed PE/SE samples is not defined. Using only PE data!\n")
                
                

            ######### Step 2, translate raw kraken into useful names
            self.script += "# Create useful kraken output \n\n";
            self.script += "if [ -e %s ]\n" % output_filename;
            self.script += "then\n\t";
            self.script += os.path.dirname(self.params["script_path"]) + "/kraken-translate \\\n\t\t";
            self.script += "--db %s \\\n\t\t" % self.params["redir_params"]["--db"]

            self.script += output_filename + " \\\n\t\t";
            self.script += "> %s.labels \n\n" % output_filename
            self.script += "fi \n\n";

            ######### Step 3, create report from kraken output (=tabular report similar to metaphlan and QIIME):
            self.script += "# Create kraken report \n\n";
            self.script += "if [ -e %s ]\n" % output_filename
            self.script += "then\n\t";
            self.script += "%s/kraken-report \\\n\t\t" % os.path.dirname(self.params["script_path"])
            self.script += "--db %s \\\n\t\t" % self.params["redir_params"]["--db"]
            
            self.script += output_filename + " \\\n\t\t";
            self.script += "> %s.report \n\n" % output_filename
            self.script += "fi \n\n";
    
            ######### Step 4, create krona report:
            if "ktImportTaxonomy_path" in self.params.keys():
                self.script += """
# Create file for ktImportTaxonomy
if [ -e %(krak_out)s ]
then
    cut -f 2,3 %(krak_out)s \\
        > %(krak_out)s.forKrona
fi

""" % {"krak_out":output_filename}
    
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
        """ Make file containing samples and target file names for use by kraken analysis R script
        """
        
        with open(self.base_dir + "kraken_files_index.txt", "w") as index_fh:
            index_fh.write("Sample\tkraken_report\n")
            for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
                index_fh.write("%s\t%s\n" % (sample,self.sample_data[sample]["classification_report"]))
                
        self.sample_data["kraken_file_index"] = self.base_dir + "kraken_files_index.txt"
        