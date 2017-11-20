# -*- coding: UTF-8 -*-
""" 

``qiime_chimera``
----------------------------


:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

.. Note:: This module was developed as part of a study led by Dr. Jacob Moran Gilad

A module for running QIIME's identify_chimeric_seqs.py:

The module can operate on the raw ``seqs.fna`` or on an aligned version. The latter is used for ChimeraSlayer and the former for usearch61


Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* A fasta file in:

    * ``sample_data["qiime.fasta"]``
    
* Alternatively, an aligned fasta file in:

    * ``sample_data["fasta.aligned"]``
    

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Puts the resulting list of chimeras in 

    * ``self.sample_data["chimeras"]``
    
* Puts the filtered fasta file in:

    * ``self.sample_data["fasta.chimera_removed"]``
    * ``self.sample_data["fasta.nucl"]``


    
.. Note:: When using ``parallel_identify_chimeric_seqs.py``, the module tries to build the scripts appropriately. **It is wise to check the parallel scripts before running them...**
  
Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "method", "usearch61 or ChimeraSlayer", "Method to use for the analysis (passed to the --chimera_detection_method of ``identify_chimeric_seqs.py``"

    
Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    q_chimera_usrch:
        module: qiime_chimera
        base: q_demult_1
        # script_path: '{Vars.qiime_path}/parallel_identify_chimeric_seqs.py'
        script_path: '{Vars.qiime_path}/identify_chimeric_seqs.py'
        method:         usearch61 # Or ChimeraSlayer. Will guess depending on existing files.
        redirects:
            # --jobs_to_start:              20
            --aligned_reference_seqs_fp:  /path/to/reference_files.otus_aligned
            --reference_seqs_fp:  /path/to/reference_files.otus
    
References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Caporaso, J.G., Kuczynski, J., Stombaugh, J., Bittinger, K., Bushman, F.D., Costello, E.K., Fierer, N., Peña, A.G., Goodrich, J.K., Gordon, J.I. and Huttley, G.A., 2010. "QIIME allows analysis of high-throughput community sequencing data". *Nature methods*, 7(5), pp.335-336.




"""



import os
import sys
import re
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "0.2.0"
class Step_qiime_chimera(Step):
   
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        self.file_tag = "qiime_chimera"
         
        
        
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
            Here you should do testing for dependency output. These will NOT exist at initiation of this instance. They are set only following sample_data updating
        """
        
        # # If does not exist 
        # try:
            # self.sample_data["qiime"]
        # except KeyError:
            # raise AssertionExcept("It seems like qiime_demult is the first qiime step. At the moment, it must come after qiime_prep...\n" )

        # try:
            # self.sample_data["qiime.fasta"]
        # except KeyError:
            # self.write_warning("fasta file does not exist. Did you intend to run a step before chimera searching?\n")

        if "otu_table" not in self.sample_data:
            self.write_warning("otu table does not exist.")


        # Guessing method to use:
        if "method" not in self.params.keys():
            self.write_warning("You did not define a method for %s\n\tTrying to guess what to use:\n")
            
            if "fasta.aligned" in self.sample_data.keys():
                # This step is being called after alignment. Will use ChimeraSlayer
                self.write_warning("Calling ChimeraSlayer on unfiltered alignment.\n\tYOU MUST MAKE SURE YOU FILTER THE ALIGNEMNT AFTER THIS STEP\n")
                
                self.method = "ChimeraSlayer";
            else:
                # This step is being called before alignment. Will use usearch61
                self.write_warning("Calling usearch61 on demulted fasta.\n")
                
                self.method = "usearch61"  
            
        else:
            self.method = self.params["method"]

                
     

        
    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        
        pass
      

    def build_scripts(self):
        

        
        # Name of specific script:
        self.spec_script_name = "_".join([self.step,self.name,self.sample_data["Title"]])

        self.script = ""

        
    
        # This line should be left before every new script. It sees to local issues.
        # Use the dir it returns as the base_dir for this step.
        use_dir = self.local_start(self.base_dir)

        ####################################################################################################
        ## 
        ## Step 1: Find chimeric sequences
        # if "env" in self.params.keys():         # Add optional environmental variables.
            # self.script += "env %s \\\n\t" % self.params["env"]
        # self.script += self.params["script_path"] + " \\\n\t"
        self.script += self.get_script_env_path()
        if re.search("parallel",self.params["script_path"]):
            # Adding a job_prefix 
            # Adding --poll_directly: Chimera determination step is automatically followed by a "filter_fasta" step. The parallel step must not return before final file is produced otherwise the filter_fasta step will throw an error. This does not affect performance since both steps are run in qsub and will not take up the terminal (if running with copy-paste through the terminal, remove the --poll_directly, if you like)
            self.script += "--job_prefix %s \\\n\t--poll_directly \\\n\t"  % self.spec_script_name


        if self.method == "ChimeraSlayer":
            if "fasta.aligned_unfiltered" in self.sample_data["qiime"].keys():
                self.script += "-i %s \\\n\t" % self.sample_data["fasta.aligned_unfiltered"]
            else:
                self.script += "-i %s \\\n\t" % self.sample_data["fasta.aligned"]

            # For ChimeraSlayer, the following must be set. If not set by user, setting to default:
            if not set(["-a","--aligned_reference_seqs_fp"]) & set(self.params["redir_params"].keys()):      # no template_fp was passed in redir_params:
                self.write_warning("--aligned_reference_seqs_fp not defined. Using default (find it with print_qiime_config.py\n")
            
            # For ChimeraSlayer, output is a file
            self.script += "-o %schimeras_CS.txt  \\\n\t" % use_dir
        else:   # Method is usearch61:
            self.script += "-i %s  \\\n\t" % self.sample_data["fasta.nucl"]
            # For usearch61, output is a folder
            self.script += "-o %s \\\n\t" % use_dir
    
        for key in self.params["redir_params"].keys():
            self.script += "%s %s \\\n\t" % (key,self.params["redir_params"][key] if self.params["redir_params"][key] else "")
        self.script += "-m %s \n\n" % self.method

        if self.method == "ChimeraSlayer":
            if "fasta.aligned_unfiltered" in self.sample_data:
                self.sample_data["chimeras"] = (self.base_dir + "chimeras_CS.txt")
            else:
                self.sample_data["chimeras"] = (self.base_dir + "chimeras_CS.txt")

        else:
            self.sample_data["chimeras"] = (self.base_dir + "chimeras.txt")


           

        ####################################################################################################
        ## 
        ## Step 2: Filter out chimeras
        
        self.script += "sleep 10\n\n"       # If running parallel chimeraslayer, it might take time for the output file to arrive in the correct location. Wait for a few seconds before continuing
        # filter_fasta.py -f seqs.fna -o seqs_chimeras_filtered.fna -s usearch_checked_chimeras/chimeras.txt -n
        # if "env" in self.params.keys():         # Add optional environmental variables.
            # self.script += "env %s \\\n\t" % self.params["env"]
        self.script += self.get_setenv_part()
        # Assuming filter_fasta.py is in the same location as chimera checking script used above...
        self.script += "%s \\\n\t" % os.sep.join([os.path.split(os.path.normpath(self.params["script_path"]))[0] , "filter_fasta.py"])
        self.script += "-f %s \\\n\t" % self.sample_data["fasta.nucl"]
        self.script += "-o " + use_dir + "seqs_chimeras_filtered.fna \\\n\t";
        self.script += "-s %s%s \\\n\t" %  (use_dir ,"chimeras_CS.txt" if self.method == "ChimeraSlayer" else "chimeras.txt")#self.sample_data["qiime"]["chimeras"]
        self.script += "-n \n\n";

        # Storing product in samples_hash, depending on the method used. At any case, storing product in active fasta
        if self.method == "ChimeraSlayer":
            if "fasta.aligned_unfiltered" in self.sample_data:
                self.sample_data["fasta.chimera_removed"] = self.base_dir + "seqs_chimeras_filtered.fna";
                self.sample_data["fasta.nucl"] = self.base_dir + "seqs_chimeras_filtered.fna";
            else:
                self.sample_data["fasta.chimera_removed"] = self.base_dir + "seqs_chimeras_filtered.fna";
                self.sample_data["fasta.nucl"] = self.base_dir + "seqs_chimeras_filtered.fna";

        else:
            self.sample_data["fasta.chimera_removed"] = self.base_dir + "seqs_chimeras_filtered.fna";
            self.sample_data["fasta.nucl"] = self.base_dir + "seqs_chimeras_filtered.fna";

       
        # Move all files from temporary local dir to permanent base_dir
        self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
     
            
        # self.stamp_dir_files(self.base_dir)
        
        self.create_low_level_script()
                    
            
            
     