# -*- coding: UTF-8 -*-
""" 
``Trim_Galore``
-----------------------

:Authors: Liron Levin
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

Short Description
~~~~~~~~~~~~~~~~~~~~~~
    A module for running Trim Galore on fastq files

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
    * puts unpaired fastq output files in the following slots:
        ``sample_data[<sample>]["fastq.F.unpaired"]``
        ``sample_data[<sample>]["fastq.R.unpaired"]``

Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "",  "", ""
    

Comments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    *  This module was tested on:
        ``Trim Galore v0.4.2``
        ``Cutadapt v1.12.1``


Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    Step_Name:                       # Name of this step
        module: Trim_Galore          # Name of the module used
        base:                        # Name of the step [or list of names] to run after [must be after a merge step]
        script_path:                 # Command for running the Trim Galore script
        qsub_params:
            -pe:                     # Number of CPUs to reserve for this analysis
        cutadapt_path:               # Location of cutadapt executable 
        redirects:
            --length:                # Parameters for running Trim Galore
            -q:                      # Parameters for running Trim Galore

References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    * Cutadapt:
        Martin, Marcel. "Cutadapt removes adapter sequences from high-throughput sequencing reads." EMBnet journal 17.1 (2011):pp-10
    * Trim Galore:
        Krueger F: Trim Galore. [http://www.bioinformatics.babraham.ac.uk/projects/]

"""
import os
import sys
import re
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Liron Levin"
__version__= "1.2.0"


class Step_Trim_Galore(Step):

    auto_redirs = "--paired".split(" ")

    def step_specific_init(self):
        self.shell = "bash"
        self.file_tag = "trim_galore.fq"
        if "--output_dir" in self.params["redir_params"] or "-o" in self.params["redir_params"]:
            raise AssertionExcept("You should not give output directory\n")

    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        
        # Assert that all samples have reads files:
        for sample in self.sample_data["samples"]:    
            if not {"fastq.F", "fastq.R", "fastq.S"} & set(self.sample_data[sample].keys()):
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
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
            
            # Name of specific script:
            self.spec_script_name = self.set_spec_script_name(sample)
            self.script = ""
            output_prefix=sample+"_trim_gal"
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(self.base_dir)

            if "fastq.F" in self.sample_data[sample] and "fastq.R" in self.sample_data[sample]:
                self.script += self.get_script_env_path()

                # Here we do the script constructing for paired end
                # Define target filenames:
                basename_F = os.path.basename(self.sample_data[sample]["fastq.F"])
                basename_R = os.path.basename(self.sample_data[sample]["fastq.R"])
                # TODO: Remove ".fq" in middle of file name
                # Setting filenames before adding output arguments to script
                fq_fn_F = use_dir + "".join([re.sub("\.\w+$", "", basename_F),
                                             "_val_1.fq"])  # The filename containing the end result. Used both in script and to set reads in $sample_params
                fq_fn_R = use_dir + "".join([re.sub("\.\w+$", "", basename_R),
                                             "_val_2.fq"])  # The filename containing the end result. Used both in script and to set reads in $sample_params
                fq_fn_F_UP = use_dir + "".join([re.sub("\.\w+$", "", basename_F),
                                                "_unpaired_1.fq"])  # The filename containing the end unpaired trimmo output
                fq_fn_R_UP = use_dir + "".join([re.sub("\.\w+$", "", basename_R),
                                                "_unpaired_2.fq"])  # The filename containing the end unpaired trimmo output
                fq_fn_F_bn = os.path.basename(fq_fn_F);
                fq_fn_R_bn = os.path.basename(fq_fn_R);
                fq_fn_F_UP_bn = os.path.basename(fq_fn_F_UP);
                fq_fn_R_UP_bn = os.path.basename(fq_fn_R_UP);

                self.script += self.get_redir_parameters_script()

                if "cutadapt_path" in self.params.keys():
                    self.script += "--path_to_cutadapt %s \\\n\t" % self.params["cutadapt_path"]
                self.script += "%s \\\n\t" % (" \\\n\t".join([self.sample_data[sample]["fastq.F"], \
                                                              self.sample_data[sample]["fastq.R"]]))
                # if "--paired" not in self.params["redir_params"].keys():
                self.script += "--paired \\\n\t"

                self.script += "-o %s \n\n" % use_dir
                # Set current active sequence files to tagged files
                self.sample_data[sample]["fastq.F"] = self.base_dir + fq_fn_F_bn
                self.sample_data[sample]["fastq.R"] = self.base_dir + fq_fn_R_bn
                self.sample_data[sample]["fastq.F.unpaired"] = self.base_dir + fq_fn_F_UP_bn
                self.sample_data[sample]["fastq.R.unpaired"] = self.base_dir + fq_fn_R_UP_bn

            elif "fastq.S" in self.sample_data[sample]:
                # Add 'env' and 'script_path':
                self.script += self.get_script_env_path()

                # Here we do the script constructing for single end
                # Define target filenames:
                basename_S = os.path.basename(self.sample_data[sample]["fastq.S"])
                # TODO: Remove ".fq" in middle of file name

                fq_fn_S = use_dir + "".join([re.sub("\.\w+$", "", basename_S),
                                             "_trimmed.fq"])  # The filename containing the end result. Used both in script and to set reads in $sample_params
                fq_fn_S_bn = os.path.basename(fq_fn_S);
                # # TODO: use existing
                # # Remove --paired and --retain_unpaired from redirects. Should not be passed if SE (Menachem)
                # for key in ["--paired","--retain_unpaired"]:
                # if key in self.params:
                # del self.params[key]
                self.script += self.get_redir_parameters_script()

                # for key in self.params["redir_params"].keys():
                # if key not in ["--paired","--retain_unpaired"]:
                # self.script += "%s %s \\\n\t" % (key,self.params["redir_params"][key])
                if "cutadapt_path" in self.params.keys():
                    self.script += "--path_to_cutadapt %s \\\n\t" % self.params["cutadapt_path"]
                self.script += "%s \\\n\t" % self.sample_data[sample]["fastq.S"]
                self.script += "-o %s \\\n\t" % use_dir
                # Generate log file:
                # self.script += ">& %s.log\n\n" % os.sep.join([use_dir.rstrip(os.sep),output_prefix])
                self.script += "\n\n"
                self.sample_data[sample]["fastq.S"] = self.base_dir + fq_fn_S_bn
            else:
                raise AssertionExcept("Unrecognized combination of fastq files")


            # Move all files from temporary local dir to permanent base_dir
            self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)


            if "spec_dir" in self.params.keys():
                self.script += "cd " + self.pipe_data["home_dir"] + "\n\n";
            
                        
            
            self.create_low_level_script()
                    
