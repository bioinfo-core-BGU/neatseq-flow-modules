# -*- coding: UTF-8 -*-
""" 
``trimmo`` :sup:`*`
-----------------------------------------------------------------

:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A module for running trimmomatic on fastq files

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
        * fastq files in at least one of the following slots:
        
            * ``sample_data[<sample>]["fastq.F"]``
            * ``sample_data[<sample>]["fastq.R"]``
            * ``sample_data[<sample>]["fastq.S"]``
        
Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
        * puts fastq output files in the following slots:
        
            * ``sample_data[<sample>]["fastq.F"|"fastq.R"|"fastq.S"]``


Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"

    "spec_dir", "path", "If *trimmomatic* must be executed within a particular directory, specify that directory here"
    "todo",     "LEADING:20 TRAILING:20", "The trimmomatic arguments"

Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::
 
    trim1:
        module: trimmo
        base: merge1
        script_path: java -jar trimmomatic-0.32.jar
        qsub_params:
            -pe: shared 20
            node: node1
        spec_dir: /path/to/Trimmomatic_dir/
        todo: LEADING:20 TRAILING:20
        redirects:
            -threads: 20

            
References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Bolger, A.M., Lohse, M. and Usadel, B., 2014. **Trimmomatic: a flexible trimmer for Illumina sequence data**. *Bioinformatics*, 30(15), pp.2114-2120.

"""

import os
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.1.0"


class Step_trimmo(Step):
    

    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        self.file_tag = "trimmo.fq"

    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        
        # Assert that all samples have reads files:
        for sample in self.sample_data["samples"]:    
            if not {"fastq.F", "fastq.R", "fastq.S"} & set(self.sample_data[sample].keys()):
                raise AssertionExcept("No read files defined\n",sample)

        pass
        
    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        
        pass
    
    def build_scripts(self):
        """ This is the actual script building function
            Most, if not all, editing should be done here 
            HOWEVER, DON'T FORGET TO CHANGE THE CLASS NAME AND THE FILENAME!
        """
        
       
        # Each iteration must define the following class variables:
            # spec_script_name
            # script
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash


            if "fastq.F" in self.sample_data[sample] and "fastq.R" in self.sample_data[sample]:
                #################### Start PE
                self.spec_script_name = "_".join([self.step,self.name,sample,"PE"])
                self.script = ""
                
                
                # This line should be left before every new script. It sees to local issues.
                # Use the dir it returns as the base_dir for this step.
                use_dir = self.local_start(self.base_dir)
                # Trimmomatic required CDing into specific dir:
                if "spec_dir" in self.params.keys():
                    self.script += "cd " + self.params["spec_dir"] + "\n\n";
                else:
                    self.write_warning("You did not supply a spec_dir param to %s. This may produce an error!\n")
                    
     
 
                # Add 'env' and 'script_path':
                self.script += self.get_script_env_path()
                

                # Here we do the script constructing for paired end
                # Define target filenames:
                basename_F = os.path.basename(self.sample_data[sample]["fastq.F"])
                basename_R = os.path.basename(self.sample_data[sample]["fastq.R"])
                # TODO: Remove ".fq" in middle of file name
                # Setting filenames before adding output arguments to script
                fq_fn_F = use_dir + ".".join([basename_F, self.file_tag])  #The filename containing the end result. Used both in script and to set reads in $sample_params
                fq_fn_R = use_dir + ".".join([basename_R, self.file_tag])  #The filename containing the end result. Used both in script and to set reads in $sample_params
                fq_fn_F_UP = use_dir + ".".join([basename_F, "unpaired",self.file_tag])   # The filename containing the end unpaired trimmo output
                fq_fn_R_UP = use_dir + ".".join([basename_R, "unpaired",self.file_tag])		#The filename containing the end unpaired trimmo output
                
                fq_fn_F_bn = os.path.basename(fq_fn_F);
                fq_fn_R_bn = os.path.basename(fq_fn_R);
                fq_fn_F_UP_bn = os.path.basename(fq_fn_F_UP);
                fq_fn_R_UP_bn = os.path.basename(fq_fn_R_UP);

                self.script +=  "PE \\\n\t";
                self.script += self.get_redir_parameters_script()
                self.script += "%s \\\n\t" % (" \\\n\t".join([self.sample_data[sample]["fastq.F"],  \
                                                              self.sample_data[sample]["fastq.R"],  \
                                                              fq_fn_F,                                      \
                                                              fq_fn_F_UP,                                   \
                                                              fq_fn_R,                                      \
                                                              fq_fn_R_UP]))
                # Add TODO line
                self.script +=  self.params["todo"] + "\n\n";
                                                                                      
                # Set current active sequence files to tagged files
                self.sample_data[sample]["fastq.F"]     = self.base_dir + fq_fn_F_bn
                self.sample_data[sample]["fastq.R"]     = self.base_dir + fq_fn_R_bn
                self.sample_data[sample]["fastq.F_UP"]  = self.base_dir + fq_fn_F_UP_bn
                self.sample_data[sample]["fastq.R_UP"]  = self.base_dir + fq_fn_R_UP_bn

                # md5sum-stamp output files:
                self.stamp_file(self.sample_data[sample]["fastq.F"])
                self.stamp_file(self.sample_data[sample]["fastq.R"])
                self.stamp_file(self.sample_data[sample]["fastq.F_UP"])
                self.stamp_file(self.sample_data[sample]["fastq.R_UP"])
                                                                                      
                                                                                          
 
                # Move all files from temporary local dir to permanent base_dir
                self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)


                if "spec_dir" in self.params.keys():
                    self.script += "cd " + self.pipe_data["home_dir"] + "\n\n";
                
                            
                self.create_low_level_script()
                #################### End PE

            if "fastq.S" in self.sample_data[sample]:
                #################### Start SE
                self.spec_script_name = "_".join([self.step,self.name,sample,"SE"])
                self.script = ""
                
                
                # This line should be left before every new script. It sees to local issues.
                # Use the dir it returns as the base_dir for this step.
                use_dir = self.local_start(self.base_dir)
                # Trimmomatic required CDing into specific dir:
                if "spec_dir" in self.params.keys():
                    self.script += "cd " + self.params["spec_dir"] + "\n\n";
                else:
                    self.write_warning("You did not supply a spec_dir param to %s. This may produce an error!\n")

                    
                # Add 'env' and 'script_path':
                self.script += self.get_script_env_path()
                

                # Here we do the script constructing for single end
                # Define target filenames:
                basename_S = os.path.basename(self.sample_data[sample]["fastq.S"])
                # TODO: Remove ".fq" in middle of file name
                
                fq_fn_S = use_dir + ".".join([basename_S, self.file_tag])          #The filename containing the end result. Used both in script and to set reads in $sample_params
                fq_fn_S_bn = os.path.basename(fq_fn_S);
                # TODO: use existing 
                self.script +=  "SE \\\n\t";
                self.script += self.get_redir_parameters_script()
                
                self.script += "%s \\\n\t" % (" \\\n\t".join([self.sample_data[sample]["fastq.S"],fq_fn_S]))
                # Add TODO line
                self.script +=  self.params["todo"] + "\n\n";
                                                                                      
                self.sample_data[sample]["fastq.S"] = self.base_dir + fq_fn_S_bn
                # md5sum-stamp output files:
                self.stamp_file(self.sample_data[sample]["fastq.S"])
                                                                                      
                                                                                          
 
                # Move all files from temporary local dir to permanent base_dir
                self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)


                if "spec_dir" in self.params.keys():
                    self.script += "cd " + self.pipe_data["home_dir"] + "\n\n";
                
                            
                self.create_low_level_script()
                #################### End SE
                
                
                
                
                
                
