# -*- coding: UTF-8 -*-
""" 
``UCSC_BW_wig`` (Included in main NeatSeq-Flow repo)
-----------------------------------------------------------------
:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A module for creating wig and bigwig files using UCSC tools:

The module creates bigwig and wig files from the current active BedGraph file.


Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* BedGraph file in the following slot:

    * ``sample_data[<sample>]["bdg"]``
    

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Puts output sam files in the following slots:

    * self.sample_data[<sample>]["bw"]
    * self.sample_data[<sample>]["wig"]
    
Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "bedGraphToBigWig_params", "*e.g.* -blockSize=10 -itemsPerSlot=20", "Parameters to pass to ``bedGraphToBigWig``"
    "bigWigToWig_params", "*e.g.* -chrom X1 -start X2 -end X3", "Parameters to pass to ``bigWigToWig``"
    "script_path", "", "Path to dir where UCSC tools are located."

    
.. note:: Set ``script_path`` to the path ro the UCSC tools, not to a specific tool!!!
    Both ``bedGraphToBigWig`` and ``bigWigToWig`` will be executed. To set specific params, use ``bedGraphToBigWig_params`` and ``bigWigToWig_params``, respectively.

Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    UCSCmap_bams:
        module: UCSC_BW_wig
        base: genCovBed
        script_path: /path/to/kentUtils/bin
        bedGraphToBigWig_params: -blockSize 10 -itemsPerSlot 20
        bigWigToWig_params: -chrom X1 -start X2 -end X3
        genome: /path/to/genome.chrom.sizes

References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Kent, W.J., Sugnet, C.W., Furey, T.S., Roskin, K.M., Pringle, T.H., Zahler, A.M. and Haussler, D., 2002. **The human genome browser at UCSC**. *Genome research*, 12(6), pp.996-1006.

"""


import os
import sys
from PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.1.0"


class Step_UCSC_BW_wig(Step):
    
    def step_specific_init(self):
        self.shell = "csh"      # Can be set to "bash" by inheriting instances

        if not "genome"  in self.params:
            raise AssertionExcept("You must pass a 'genome' parameter!")
        
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        

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
            # self.spec_script_name
            # self.script
        
        self.base_dir    
    
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash

            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)

            # Name of specific script:
            self.spec_script_name = "_".join([self.step,self.name,sample])
            self.script = ""
            
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)

            
            # Define input and output files
            input_file = self.sample_data[sample]["bdg"]
            
            output_file_bw  = "%s.bw"  % os.path.basename(input_file)
            output_file_wig = "%s.wig" % os.path.basename(input_file)
            
            
            # Creating bedGraphToBigWig script:
            self.script += "# Converting bdg to bigWig:\n\n"
            # Adding env, if it exists:
            if "env" in self.params.keys():         # Add optional environmental variables.
                script_const += "env %s \\\n\t" % self.params["env"]
            # Adding bedGraphToBigWig executable
            self.script += "%s%sbedGraphToBigWig \\\n\t" % (self.params["script_path"],os.sep)
            # Adding parameters, if the exist
            if "bedGraphToBigWig_params" in self.params:
                self.script += "%s \\\n\t" % self.params["bedGraphToBigWig_params"]
            # Adding input, genome and output files:
            self.script += "%s \\\n\t" % (input_file)
            self.script += "%s \\\n\t" % self.params["genome"]
            self.script += "%s\n\n" % (use_dir + output_file_bw)
            
            # Creating bigWigToWig script:
            self.script += "# Converting bigWig to wig:\n\n"
            self.script += "%s%sbigWigToWig \\\n\t" % (self.params["script_path"],os.sep)
            if "bigWigToWig_params" in self.params:
                self.script += "%s \\\n\t" % self.params["bigWigToWig_params"]
            self.script += "%s \\\n\t" % (use_dir + output_file_bw)
            self.script += "%s\n\n" % (use_dir + output_file_wig)
            
            
            self.sample_data[sample]["bw"]  = "%s%s" % (sample_dir, output_file_bw)
            self.sample_data[sample]["wig"] = "%s%s" % (sample_dir, output_file_wig)
    
            # Stamping output files:
            self.stamp_file(self.sample_data[sample]["bw"])
            self.stamp_file(self.sample_data[sample]["wig"])
        
            # Move all files from temporary local dir to permanent base_dir
            self.local_finish(use_dir,sample_dir)       # Sees to copying local files to final destination (and other stuff)
       
            
            
            self.create_low_level_script()
                    
        
