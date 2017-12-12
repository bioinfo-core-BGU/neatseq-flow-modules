# -*- coding: UTF-8 -*-
""" 
``kraken_biom``
--------------------------------

:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

.. Note:: This module was developed as part of a study led by Dr. Jacob Moran Gilad

A module for running ``kraken-biom`` (https://github.com/smdabdoub/kraken-biom)


Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Kraken reports:

    * ``sample_data[<sample>]["classification_report"]``


Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Puts the resulting *biom* output files in:  

    * ``self.sample_data["kraken_biom"]``
    * ``self.sample_data["biom_table"]``
    * ``self.sample_data["biom_table_tsv"]`` (if ``skip_tsv`` is not set)
    

    
Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "skip_tsv", "", "Set if you do not want to convert the report into tsv format."
    "biom_path", "/path/to/biom", "The path to biom. This is required for conversion to tsv"

    
Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    kraken_biom1:
        module:             kraken_biom
        base:               kraken1
        script_path:        '{Vars.paths.kraken_biom}'
        # skip_tsv:
        biom_path:          '{Vars.paths.biom}'
        redirects:
            --max:          D 
            --min:          S 
            --gzip:

References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

https://github.com/smdabdoub/kraken-biom

"""


import os
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept

from pkg_resources import resource_filename


__author__ = "Menachem Sklarz"
__version__ = "1.2.0"


class Step_kraken_biom(Step):
    """ A class that defines a pipeline step name (=instance).
    """
    
    def step_specific_init(self):
        """ Called on intiation
            Good place for parameter testing.
            Wrong place for sample data testing
        """
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        self.file_tag = ".kraken.out"

        if "skip_tsv" not in self.params:
            if "biom_path" not in self.params:
                self.write_warning("You did not set 'biom_path'. Using 'biom' from $PATH. This is not advisable. ")
                self.params["biom_path"] = ""
        else:
            if "biom_path" in self.params:
                self.write_warning("Skipping tsv creation even though you passed 'biom_path'")
        
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        
        for sample in self.sample_data["samples"]:
            if not "classification_report" in self.sample_data[sample]:
                raise AssertionExcept("Sample does not own a 'classification_report'", sample)
        
    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        

    
    def build_scripts(self):
        """ This is the actual script building function
            Most, if not all, editing should be done here 
            HOWEVER, DON'T FORGET TO CHANGE THE CLASS NAME AND THE FILENAME!
        """
        
        self.spec_script_name = "_".join([self.step,self.name,self.sample_data["Title"]])
        self.script = ""

        # This line should be left before every new script. It sees to local issues.
        # Use the dir it returns as the base_dir for this step.
        use_dir = self.local_start(self.base_dir)

        
        prefix = self.sample_data["Title"]
        # Decide on ".gz" suffix by whether --gzip is in redirects:
        dotgz = ".gz" if "--gzip" in self.params["redir_params"] else ""

        use_biom_path = "{base_dir}{prefix}.kraken.biom{suffix}".format(base_dir = use_dir,
                                                                        prefix = prefix,
                                                                        suffix = dotgz)
        
        self.script += self.get_script_const()
        self.script += "--output_fp {use_biom_path} \\\n\t".format(use_biom_path = use_biom_path)
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
            self.script += "{kraken_report} \\\n\t".format(kraken_report = self.sample_data[sample]["classification_report"])

        self.script.rstrip("\\\n\t")

        
        # Storing the output file in $samples_hash
        self.sample_data["kraken_biom"]        = "{base_dir}{prefix}.kraken.biom{suffix}".format(base_dir = self.base_dir,
                                                                                                 prefix = prefix,
                                                                                                 suffix = dotgz)
        self.sample_data["biom_table"]        = self.sample_data["kraken_biom"]
            
        self.stamp_file(self.sample_data["kraken_biom"])




        #######################################################################################
        ## 
        ## Step 2: Creating biom table summary
        
        # cmd_text = self.get_script_env_path() 
        if "skip_summary" not in self.params:
            cmd_text = """
    {biom_path} summarize-table \\
        -i {biom} \\
        -o {biom_summary} 
""".format(biom       = "{use_biom_path}".format(use_biom_path = use_biom_path),
            biom_summary  = "{use_biom_path}.summary.txt".format(use_biom_path = use_biom_path),
            biom_path = self.params["biom_path"])
           
            
            self.script += """
# Create biom table in table format

if [ -e {biom} ]
then
    {cmd_text}
fi

""".format(biom = "{use_biom_path}".format(use_biom_path = use_biom_path),
           cmd_text = cmd_text)
           
            
            # Store location of the biom_table summary:
            self.sample_data["biom_table_summary"] = "{base_dir}{prefix}.kraken.biom.summary.txt".format(base_dir = self.base_dir,prefix = prefix)

            
        ################################################################################################
        ## 
        ## Step 3: Creating biom table in table format
        
        if "skip_tsv" not in self.params:
            
            cmd_text = """
    {biom_path} convert \\
        -i {biom} \\
        -o {biom_tsv} \\
        --to-tsv \\
        --header-key taxonomy \\
        --output-metadata-id \"Consensus Lineage\"
""".format(biom       = "{use_biom_path}".format(use_biom_path = use_biom_path),
            biom_tsv  = "{use_biom_path}.tsv".format(use_biom_path = use_biom_path),
            biom_path = self.params["biom_path"])

            self.script += """
# Create biom table in table format

if [ -e {biom} ]
then
    {cmd_text}
fi

""".format(biom = "{use_biom_path}".format(use_biom_path = use_biom_path),
           cmd_text = cmd_text)
            
            # Store location of the tsv biom_table:
            self.sample_data["biom_table_tsv"] = "{base_dir}{prefix}.kraken.biom.tsv".format(base_dir = self.base_dir,prefix = prefix)


        
        # Move all files from temporary local dir to permanent base_dir
        self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
                    
        
        self.create_low_level_script()
                    

                    
    def make_sample_file_index(self):
        """ Make file containing samples and target file names for use by kraken analysis R script
        """
        
        pass