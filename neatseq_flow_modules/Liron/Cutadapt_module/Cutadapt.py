# -*- coding: UTF-8 -*-
""" 
``Cutadapt``
-----------------------

:Authors: Levin Liron
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

Short Description
~~~~~~~~~~~~~~~~~~~~~~~~
    A module for running cutadapt on fastqc files

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    * fastq files in at least one of the following slots:
        ``sample_data[<sample>]["fastq.F"]``
        ``sample_data[<sample>]["fastq.R"]``
        ``sample_data[<sample>]["fastq.S"]``

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    * puts fastq output files in the following slots:
        ``sample_data[<sample>]["fastq.F"]``
        ``sample_data[<sample>]["fastq.R"]``
        ``sample_data[<sample>]["fastq.S"]``

Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "",  "", ""
    

Comments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    *  This module was tested on:
        ``Cutadapt v1.12.1``

Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    Step_Name:                       # Name of this step
        module: Cutadapt             # Name of the module used
        base:                        # Name of the step [or list of names] to run after [must be after a merge step]
        script_path:                 # Command for running the Cutadapt script
        paired:                      # Analyse Forward and Reverse reads together.
        Demultiplexing:              # Use to Demultiplex the adaptors, needs to be in the format of name=adaptor_seq
        qsub_params:
            -pe:                     # Number of CPUs to reserve for this analysis
        redirects:
            --too-short-output:      # will replace @ with the location of the sample dir  [e.g. @too_short.fq] 
            -a:                      # Use to trim poly A in SE reads [e.g. "A{100} -A T{100}"]

References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Martin, Marcel. "Cutadapt removes adapter sequences from high-throughput sequencing reads." EMBnet. journal 17.1 (2011): pp-10

"""
import os
import sys
import re
from copy import *
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Liron Levin"
__version__= "1.2.0"

class Step_Cutadapt(Step):

    

    def step_specific_init(self):
        self.shell = "bash"
        self.file_tag = ".fq"

        if "--output_dir" in self.params["redir_params"] or "-o" in self.params["redir_params"]:
            raise AssertionExcept("you should not give an output directory\n")

    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        
        # Assert that all samples have reads files:
        for sample in self.sample_data["samples"]:   
            if  not any(reads in ["fastq.F", "fastq.R", "fastq.S"] for reads in self.sample_data[sample].keys()):
                raise AssertionExcept("No read files\n",sample)
 
            
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
        original_self_params_redir_params = deepcopy(self.params["redir_params"])
        
        
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
            
            # Name of specific script:
            self.spec_script_name = self.set_spec_script_name(sample)
            self.script = ""
            output_prefix=sample+"cutadapt"
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(self.base_dir)
            
            for direction in self.sample_data[sample]["type"]:    # Iterate over Forward and single, if they exist in sample_data


                if direction == "PE":
                    # Here we do the script constructing for paired end
                    # raise AssertionExcept(self.sample_data)
                    # Define target filenames:
                    basename_F = os.path.basename(self.sample_data[sample]["fastq.F"])
                    basename_R = os.path.basename(self.sample_data[sample]["fastq.R"])
                    # TODO: Remove ".fq" in middle of file name
                    # Setting filenames before adding output arguments to script
                    if "Demultiplexing" in self.params.keys():
                        fq_fn_F = use_dir + "".join([re.sub("\.\w+$","",basename_F ), "{name}_cutadapt_1.fq"])  #The filename containing the end result. Used both in script and to set reads in $sample_params    
                        fq_fn_R = use_dir + "".join([re.sub("\.\w+$","",basename_R ), "{name}_cutadapt_2.fq"])  #The filename containing the end result. Used both in script and to set reads in $sample_params
                    else:
                        fq_fn_F = use_dir + "".join([re.sub("\.\w+$","",basename_F ), "_cutadapt_1.fq"])  #The filename containing the end result. Used both in script and to set reads in $sample_params    
                        fq_fn_R = use_dir + "".join([re.sub("\.\w+$","",basename_R), "_cutadapt_2.fq"])  #The filename containing the end result. Used both in script and to set reads in $sample_params
                    fq_fn_F_bn = os.path.basename(fq_fn_F)
                    fq_fn_R_bn = os.path.basename(fq_fn_R)

                    if ("paired" in self.params.keys())&("Demultiplexing" not in self.params.keys()):
                        # Add 'env' and 'script_path':
                        self.script += "("+self.get_script_env_path()
                        for key in original_self_params_redir_params.keys():
                            if isinstance(original_self_params_redir_params[key], str) and \
                                    original_self_params_redir_params[key].startswith("@"):
                                self.params["redir_params"][key]=os.path.join(use_dir,original_self_params_redir_params[key].replace("@",""))
                        self.script += self.get_redir_parameters_script()
                        self.script += "-o %s \\\n\t" % fq_fn_F 
                        self.script += "-p %s \\\n\t" % fq_fn_R
                        self.script += "%s \\\n\t" % self.sample_data[sample]["fastq.F"]
                        self.script += "%s \\\n\t"   % self.sample_data[sample]["fastq.R"]
                        self.script += "> %s.log ) >& %%s.out\n\n"  % os.sep.join([use_dir.rstrip(os.sep),sample+"_cutadapt"]) \
                                                                    % os.sep.join([use_dir.rstrip(os.sep),sample+"_cutadapt"])
                    else:
                        for files_types in ["fastq.F","fastq.R"]:
                            # Add 'env' and 'script_path':
                            self.script +="("+ self.get_script_env_path()
                            for key in original_self_params_redir_params.keys():
                                if original_self_params_redir_params[key].startswith("@"):
                                    self.params["redir_params"][key]=os.path.join(use_dir,original_self_params_redir_params[key].replace("@",files_types+"_"))

                            self.script += self.get_redir_parameters_script()
                            if files_types=="fastq.F":
                                self.script += "-o %s \\\n\t" % fq_fn_F
                                self.script += "%s \\\n\t" % self.sample_data[sample][files_types]
                                self.script += "> %s.log ) >& %%s.out\n\n"  % fq_fn_F \
                                                                            % fq_fn_F                            
                            else:
                                self.script += "-o %s \\\n\t" % fq_fn_R
                                self.script += "%s \\\n\t" % self.sample_data[sample][files_types]
                                self.script += "> %s.log ) >& %%s.out\n\n"  % fq_fn_R \
                                                                            % fq_fn_R                            
                elif direction=="SE":
                    # Here we do the script constructing for single end
                    # Define target filenames:
                    basename_S = os.path.basename(self.sample_data[sample]["fastq.S"])
                    if "Demultiplexing" in self.params.keys():
                        fq_fn_S = use_dir + "".join([re.sub("\.\w+$","",basename_S ), "{name}_cutadapt_trimmed.fq"])  #The filename containing the end result. Used both in script and to set reads in $sample_params    
                    else:
                        fq_fn_S = use_dir + "".join([re.sub("\.\w+$","",basename_S ), "_cutadapt_trimmed.fq"])          #The filename containing the end result. Used both in script and to set reads in $sample_params
                    fq_fn_S_bn = os.path.basename(fq_fn_S);
                    # Add 'env' and 'script_path':
                    self.script += "("+self.get_script_env_path()
                    for key in original_self_params_redir_params.keys():
                        if original_self_params_redir_params[key].startswith("@"):
                            self.params["redir_params"][key]=os.path.join(use_dir,original_self_params_redir_params[key].replace("@",""))
                    self.script += self.get_redir_parameters_script()
                    self.script += "-o %s \\\n\t" %  fq_fn_S
                    self.script += "%s \\\n\t" % self.sample_data[sample]["fastq.S"]
                    self.script += "> %s.log ) >& %%s.out\n\n"  % fq_fn_S \
                                                                % fq_fn_S

                else: # direction=="Reverse". Ignore (included in "Forward")
                    pass
                    
                    

                

                if direction == "PE":
                    if "Demultiplexing" not in self.params.keys():
                        #Set current active sequence files to tagged files
                        self.sample_data[sample]["fastq.F"] = self.base_dir + fq_fn_F_bn
                        self.sample_data[sample]["fastq.R"] = self.base_dir + fq_fn_R_bn
                    
                elif direction == "SE":
                    if "Demultiplexing" not in self.params.keys():
                        self.sample_data[sample]["fastq.S"] = self.base_dir + fq_fn_S_bn
                        
                    
            # Move all files from temporary local dir to permanent base_dir
            self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)


            if "spec_dir" in self.params.keys():
                self.script += "cd " + self.pipe_data["home_dir"] + "\n\n";
            
                        
            
            self.create_low_level_script()
                    
