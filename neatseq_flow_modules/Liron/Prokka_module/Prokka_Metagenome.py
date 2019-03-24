# -*- coding: UTF-8 -*-
""" 
``Prokka Metagenome``
-----------------------

:Authors: Liron Levin
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.


Short Description
~~~~~~~~~~~~~~~~~~~~~
    Runs Prokka on Metagenome's Bins

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    * For each Sample, a bins directory with fasta format files in:
        ``self.sample_data[sample]["bins_dir"]``
    * For project level , a bins directory with fasta format files in:
        ``self.sample_data["project_data"]["bins_dir"]``
        

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    * For each Sample, puts the location of the Sample's GFF file in:
        ``self.sample_data[sample]["GFF"]``
    * For each Sample, puts the location of the Sample's identified genes files in:    
        ``self.sample_data[sample]["fasta.nucl_dir"]``
    * For each Sample, puts the location of the Sample's identified genes [translated] file in:    
        ``self.sample_data[sample]["fasta.prot_dir"]``
    * if generate_GFF_dir option exist, puts the directory location of all Samples GFFs in:
        ``self.sample_data["project_data"]["GFF_dir"]``

Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "generate_GFF_dir",  "", "Create GFF directory"

Comments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    Step_Name:                                  # Name of this step
        module: Prokka_Metagenome               # Name of the module to use
        base:                                   # Name of the step [or list of names] to run after [must be after a binning generating step]
        script_path:                            # Command for running Prokka 
        setenv:                                 # env parameters that needs to be in the PATH for running this module
        create_genus_db:                        # If protein DB is provided using the '--proteins' redirect argument, 
                                                # a blast database will be created and added to the genus db  ['--usegenus' redirect argument]
                                                # The current analysis will automatically  use this database!
                                                # IMPORTENT: you will need write permission in the prokka's installation location.
                                                # The database will be deleted in the end of the analysis!! 
        qsub_params:
            -pe:                                # Number of CPUs to reserve for this analysis
        redirects:
            --cpus:                             # parameters for running Prokka
            --force:                            # parameters for running Prokka
            --genus:                            # parameters for running Prokka
            --kingdom:                          # parameters for running Prokka
            --proteins:                         # Use the location of a protein DB [FASTA] for extra annotation or use "VFDB" to use the module VFDB built-in virulence/resistance DB  

References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Seemann, Torsten. "Prokka: rapid prokaryotic genome annotation." Bioinformatics 30.14 (2014): 2068-2069.‏
"""
import os
import sys
import re
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Liron Levin"
__version__= "1.2.0"

class Step_Prokka_Metagenome(Step):
    
    
    def step_specific_init(self):
        """ Called on intiation
            Good place for parameter testing.
            Wrong place for sample data testing
        """
        self.shell = "bash"
        self.file_tag = ".gff"
        import inspect
        self.module_location=os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))

    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        
        if "scope" not in list(self.params.keys()):
            self.params["scope"] = 'sample'
            
        if "project" in self.params["scope"]:   
            # Testing for existance of binning data
            assert "bins_dir" in list(self.sample_data["project_data"].keys()), \
                "In %s:\tThere are no project level binning results (bins_dir).\n" % (self.get_step_name())
        else:
            for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
                # Testing for existance of binning data
                assert "bins_dir" in list(self.sample_data[sample].keys()), \
                    "In %s:\tThere are no binning results (bins_dir) for sample %s.\n" % (self.get_step_name(), sample)
        pass
        

    def create_spec_preliminary_script(self):
        """ Add script to run BEFORE all other steps
        """
        #initiating new script 
        self.script = ""
        if ("create_genus_db" in list(self.params.keys()))&("--proteins" in list(self.params["redir_params"].keys())):
            if self.params["create_genus_db"]==None:
                self.params["create_genus_db"]='Proteins_db'
            self.params["create_genus_db"]=self.params["create_genus_db"].capitalize()
            self.script += "db_dir=$(which %s) \n\n" % self.params["script_path"]
            self.script += "db_dir=${db_dir%%/*}/../db/genus/%s \n\n" % self.params["create_genus_db"]
            self.script += "makeblastdb -max_file_sz 10G  -hash_index -dbtype prot -in %s -out $db_dir" % self.params["redir_params"]["--proteins"]
            self.params["redir_params"]["--genus"]=self.params["create_genus_db"]
            self.params["redir_params"]["--usegenus"]=None
            self.params["redir_params"].pop('--proteins', None)
        pass
        
    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        self.script = ""
        if ("create_genus_db" in list(self.params.keys()))&(self.params["redir_params"]!=None):
            self.script += "db_dir=$(which %s) \n\n" % self.params["script_path"]
            self.script += "db_dir=${db_dir%%/*}/../db/genus/%s \n\n" % self.params["create_genus_db"]
            self.script += "rm -f  $db_dir.* \n\n" 
            
            
    def build_scripts(self):
        """ This is the actual script building function
            Most, if not all, editing should be done here 
            HOWEVER, DON'T FORGET TO CHANGE THE CLASS NAME AND THE FILENAME!
        """
        if "project" in self.params["scope"]:
            self.build_scripts_byproject()
        else:
            self.build_scripts_bysample()
        
    def build_scripts_byproject(self):
        
        #Make a dir for the GFF files:
        GFF_dir = self.make_folder_for_sample("GFF")
        self.sample_data["project_data"]["GFF_dir"]=GFF_dir
        #Make a dir for the nucl fasta files:
        nucl_dir = self.make_folder_for_sample("nucl")
        self.sample_data["project_data"]["nucl_dir"]=nucl_dir
        #Make a dir for the prot fasta files:
        prot_dir = self.make_folder_for_sample("prot")
        self.sample_data["project_data"]["prot_dir"]=prot_dir

        # Name of specific script:
        self.spec_script_name = self.set_spec_script_name()
        self.script = ""
        # This line should be left before every new script. It sees to local issues.
        # Use the dir it returns as the base_dir for this step.
        use_dir = self.local_start(self.base_dir)
        
        self.script += "for fasta_file in $(ls -p " + self.sample_data["project_data"]["bins_dir"]  + " | grep -v / ); do \n\n"

             
        # Define output filename 
        output_filename = "".join([self.sample_data["Title"],'_','${fasta_file%.*}' ])

        self.script += self.get_script_const()
        if "--proteins VFDB" in self.script:
            if "Virulence_Resistance.fasta" in os.listdir(self.module_location):
                self.script=self.script.replace("--proteins VFDB","--proteins %s" % os.path.join(self.module_location,"Virulence_Resistance.fasta") )
            else:
                raise AssertionExcept("The file %s is not found in the Prokka module directory" % "Virulence_Resistance.fasta" )


        self.script += "--outdir %s \\\n\t"   % use_dir
        self.script += "--locustag %s \\\n\t" % output_filename
        self.script += "--strain %s \\\n\t"   % output_filename
        self.script += "--prefix %s \\\n\t"   % output_filename
        self.script += "%s \n\n" % os.path.join( self.sample_data["project_data"]["bins_dir"] , '$fasta_file' )
       
        self.script += "done \n\n"
        
        # Wrapping up function. Leave these lines at the end of every iteration:
        self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
        
        self.script += "mkdir -p %s && mv %%s  %%%%s \n\n" % prot_dir % os.path.join(self.base_dir,"*.faa") % prot_dir
        self.script += "mkdir -p %s && mv %%s  %%%%s \n\n" % nucl_dir % os.path.join(self.base_dir,"*.ffn") % nucl_dir
        self.script += "mkdir -p %s && mv %%s  %%%%s \n\n" % GFF_dir  % os.path.join(self.base_dir,"*.gff") % GFF_dir 
        
        
        self.create_low_level_script()
    
        
        
        
    def build_scripts_bysample(self):
        # Each iteration must define the following class variables:
            # spec_script_name
            # script
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
            
            # Name of specific script:
            self.spec_script_name = self.set_spec_script_name(sample)
            self.script = ""

            self.script +='for fasta_file in $(ls -p ' + self.sample_data[sample]["bins_dir"]  + ' | grep -v / ); do\n'

            
            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)
            
            
            #Make a dir for the GFF files:
            GFF_dir = os.path.join(sample_dir,"GFF")
            self.sample_data[sample]["GFF_dir"]=GFF_dir
            #Make a dir for the nucl fasta files:
            nucl_dir = os.path.join(sample_dir,"nucl")
            self.sample_data[sample]["nucl_dir"]=nucl_dir
            #Make a dir for the prot fasta files:
            prot_dir = os.path.join(sample_dir,"prot")
            self.sample_data[sample]["prot_dir"]=prot_dir
            
            
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)
                
                
            # Define output filename 
            output_filename = "".join([sample,'_','${fasta_file%.*}' ])

            self.script += self.get_script_const()
            if "--proteins VFDB" in self.script:
                if "Virulence_Resistance.fasta" in os.listdir(self.module_location):
                    self.script=self.script.replace("--proteins VFDB","--proteins %s" % os.path.join(self.module_location,"Virulence_Resistance.fasta") )
                else:
                    raise AssertionExcept("The file %s is not found in the Prokka module directory" % "Virulence_Resistance.fasta" )


            self.script += "--outdir %s \\\n\t"   % use_dir
            self.script += "--locustag %s \\\n\t" % output_filename
            self.script += "--strain %s \\\n\t"   % output_filename
            self.script += "--prefix %s \\\n\t"   % output_filename
            self.script += "%s \n\n" % os.path.join(self.sample_data[sample]["bins_dir"] , '$fasta_file' )
           
            self.script += 'done\n'
            
            # Wrapping up function. Leave these lines at the end of every iteration:
            self.local_finish(use_dir,sample_dir)       # Sees to copying local files to final destination (and other stuff)
            
            self.script += "mkdir -p %s && mv %%s  %%%%s \n\n" % prot_dir % os.path.join(sample_dir,"*.faa") % prot_dir
            self.script += "mkdir -p %s && mv %%s  %%%%s \n\n" % nucl_dir % os.path.join(sample_dir,"*.ffn") % nucl_dir
            self.script += "mkdir -p %s && mv %%s  %%%%s \n\n" % GFF_dir  % os.path.join(sample_dir,"*.gff") % GFF_dir
            
            
            self.create_low_level_script()
