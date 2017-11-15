# -*- coding: UTF-8 -*-
""" 
``CEAS``
-------------------


:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A module for running CEAS:

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Files in the following slots:

    * ``self.sample_data[<sample>]["peak_bed"]`` - Sample ``peak_bed`` file
    * ``self.sample_data[<sample>]["wig"]`` - An appropriate ``wig`` file 


Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Puts CEAS output files in the following slots:

    * ``sample_data[sample]["CEAS.xls"]``
    * ``sample_data[sample]["CEAS.R"]``
    * ``sample_data[sample]["CEAS.plots"]``

    
Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10





Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    CEAS1:
        module: CEAS
        base: UCSC_BW_to_wig
        script_path: /path/to/bin/ceas
        redirects:
            -g: /path/to/hg19.refGene

References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Shin, H., Liu, T., Manrai, A.K. and Liu, X.S., 2009. **CEAS: cis-regulatory element annotation system**. *Bioinformatics*, 25(19), pp.2605-2606.


"""

import os
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "0.2.0"
class Step_CEAS(Step):
    """ A class that defines a pipeline step name (=instance).
    """
    


    def step_specific_init(self):
        self.shell = "csh"      # Can be set to "bash" by inheriting instances
        self.file_tag = "CEAS"

    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
            if sample in self.sample_data["Controls"].values():      
                continue        # This is a control sample. Will not contain a "chip_seq" slot.
        
            if not "wig" in self.sample_data[sample]:
                raise AssertionExcept("chip_seq dict does not contain a wig file\n",sample)
            if not "peak_bed" in self.sample_data[sample]:
                raise AssertionExcept("chip_seq dict does not contain a peak bed file.\n",sample)
            
        
        
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
        for sample in self.sample_data["Controls"].keys():      # Getting list of samples out of Controls dict.

            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)

            # Name of specific script:
            self.spec_script_name = "_".join([self.step,self.name,sample])
            self.script = ""

            # Name of control sample:
            control = self.sample_data["Controls"][sample]

            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)

            # Defined full path to output filename
            output_filename = "%s.%s"               % (sample, self.file_tag)
            
                
                
            self.script += self.get_script_const()

            self.script += "--name %s%s \\\n\t"     % (use_dir,output_filename)
            self.script += "-b %s \\\n\t"           % self.sample_data[sample]["peak_bed"]
            self.script += "-w %s \n\n"             % self.sample_data[sample]["wig"]

            # Adding plotting with R: (Assumes R is defined in PATH!)
            self.script += "R --vanilla < %s%s.R\n\n"   % (use_dir,output_filename)
            # Store CEAS output!
            
            self.sample_data[sample]["CEAS.xls"]    = "%s%s.xls" % (sample_dir,output_filename)
            self.sample_data[sample]["CEAS.R"]      = "%s%s.R"   % (sample_dir,output_filename)
            self.sample_data[sample]["CEAS.plots"]  = "%s%s.pdf" % (sample_dir,output_filename)
            
            self.stamp_file(sample_dir)

            self.stamp_file(self.sample_data[sample]["CEAS.xls"])
            self.stamp_file(self.sample_data[sample]["CEAS.R"])
            self.stamp_file(self.sample_data[sample]["CEAS.plots"])
            
        
            # Wrapping up function. Leave these lines at the end of every iteration:
            self.local_finish(use_dir,sample_dir)       # Sees to copying local files to final destination (and other stuff)

            
            self.create_low_level_script()
                    
