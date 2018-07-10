# -*- coding: UTF-8 -*-
""" 
``NGSplot``
-----------------------


:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A module for running NGSplot:

Runs NGSplot on existing **sorted** BAM files. 

Please make sure the BAM is sorted, such as following the ``samtools`` module

If this is a ChIP-seq experiment and you have controls defined, it will also run NGSplot for the sample:control comparison.

At the moment, the module works only at the sample scope. (BAM files in the project scope are rare!)

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* BAM files in the following slots:

    * ``sample_data[<sample>]["bam"]``

    

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Puts output NGS reports in the following slots:

    * ``self.sample_data[<sample>]["NGSplot"]``

* For ChIP-seq data, puts comparison reports in
    
    * ``self.sample_data[<sample>]["NGSplot_vs_control"]``
    

Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "setenv", "NGSPLOT=/path/to/ngsplot", "Running NGSplot requires setting this EV."

Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    NGSplot_genebody:
        module:             NGSplot
        base:               sam_base
        script_path:        Rscript /path/to/ngsplot-2.61/bin/ngs.plot.r
        setenv:             NGSPLOT=/path/to/ngsplot-2.61
        redirects:
            -G:             mm10
            -R:             genebody
            -P:             20
            -GO:            hc
        qsub_params:
            -pe:            shared 20

References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Shen, L., Shao, N., Liu, X. and Nestler, E., 2014. **ngs.plot: Quick mining and visualization of next-generation sequencing data by integrating genomic databases**. *BMC genomics*, 15(1), p.284.

"""
import os
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.2.0"


class Step_NGSplot(Step):

   
    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        # self.file_tag = "Bowtie_mapper"
        
        # REquire -G  and -R parameters. 
        

    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """

        # # Assert there is mapping data and a sorted bam in particular:
        # for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
            # if not "mapping" in self.sample_data[sample]["fastq"]:
                # raise AssertionExcept("No mapping data defined\n", sample)
            # if not "unsorted_bam" in self.sample_data[sample]["fastq"]["mapping"]:
                # raise AssertionExcept("No sorted bam for sample %s", sample)
                # # If there is an unsorted BAM, then there is a sorted one as well :)


        # pass
        
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
            # self.spec_script_name
            # self.script
        
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash

            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)

            # Name of specific script:
            self.spec_script_name = self.set_spec_script_name(sample)
            self.script = ""
            
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)

            # Define location and prefix for output files:
            output_prefix = sample + "_NGSplot"
            
            # Get constant part of script:
            self.script += self.get_script_const()
            # Add bam:
            self.script += "-C %s \\\n\t" % self.sample_data[sample]["bam"]
            # Add sample name:
            self.script += "-O %s \n\n" % (use_dir + output_prefix)
            
            self.sample_data[sample]["NGSplot"] = sample_dir + output_prefix
    
        
            # Move all files from temporary local dir to permanent base_dir
            self.local_finish(use_dir,sample_dir)       # Sees to copying local files to final destination (and other stuff)
            
            self.create_low_level_script()
             
            # Sorted bam should exist by assert above. 
            # Add option to pass a configuration file. Check the NGSplot help pages
            if "Controls" in self.sample_data.keys() and sample in self.sample_data["Controls"].keys():   # This is a 'sample' (i.e. not a control)
                # Name of specific script:
                self.spec_script_name = self.jid_name_sep.join([self.step,self.name,sample,"vs_control"])
                self.script = ""
                use_dir = "_".join([use_dir,"vs_control"])
                output_prefix = "_".join([output_prefix,"vs_control"])
                # Name of control for this sample:
                control = self.sample_data["Controls"][sample]



                # Get constant part of script:
                self.script += self.get_script_const()
                # Add bam:
                self.script += "-C %s:%s \\\n\t" % (self.sample_data[sample]["bam"],
                                                    self.sample_data[control]["bam"])
                # Add sample name:
                self.script += "-O %s \n\n" % (use_dir + output_prefix)
                
                self.sample_data[sample]["NGSplot_vs_control"] = sample_dir + output_prefix
                
                # self.stamp_dir_files(sample_dir)
            
                # Move all files from temporary local dir to permanent base_dir
                self.local_finish(use_dir,sample_dir)       # Sees to copying local files to final destination (and other stuff)
                
                self.create_low_level_script()
                 

