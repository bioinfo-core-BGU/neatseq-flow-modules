#!/fastspace/bioinfo_apps/python-2.7_SL6/bin/python
# -*- coding: UTF-8 -*-
""" 
Module ``Prokka``
-----------------------

:Authors: Liron Levin
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

.. Note:: This module was developed as part of a study led by Dr. Jacob Moran Gilad

SHORT DESCRIPTION
    Runs Prokka on all samples

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    For each Sample, a fasta.nucl file type [e.g. an assembly result] in:
        sample_data[sample]["fasta.nucl"]

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    For each Sample, puts the location of the Sample's GFF file in:
        sample_data[sample]["GFF"]
    For each Sample, puts the location of the Sample's identified genes file in:    
        sample_data[sample]["fasta.nucl"]
    For each Sample, puts the location of the Sample's identified genes [translated] file in:    
        sample_data[sample]["fasta.prot"]
    if generate_GFF_dir option exist, puts the directory location of all Samples GFFs in:
        sample_data["GFF_dir"]

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
        module: Prokka                          # Name of the module to use
        base:                                   # Name of the step [or list of names] to run after [must be after a fasta file generator step like an assembly program or start the analysis with fasta files]
        script_path:                            # Command for running Prokka 
        env:                                    # env parameters that needs to be in the PATH for running this module
        qsub_params:
            -pe:                                # Number of CPUs to reserve for this analysis
        generate_GFF_dir:                       # Create GFF directory
        redirects:
            --cpus:                             # parameters for running Prokka
            --force:                            # parameters for running Prokka
            --genus:                            # parameters for running Prokka
            --kingdom:                          # parameters for running Prokka
            --proteins:                         # Use the location of a protein DB [FASTA] for extra annotation or use "VFDB" to use the module VFDB built-in virulence/resistance DB  

    
References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"""
import os
import sys
import re
from PLC_step import Step,AssertionExcept


__author__ = "Liron Levin"

class Step_Prokka(Step):
    
    
    def step_specific_init(self):
        """ Called on intiation
            Good place for parameter testing.
            Wrong place for sample data testing
        """
        self.shell = "csh"      # Can be set to "bash" by inheriting instances
        self.file_tag = ".gff"
        import inspect
        self.module_location=os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))

    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        
            
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
             # Testing for existance of assembly data
            assert "fasta.nucl" in self.sample_data[sample].keys(), \
                "In %s:\tThere are no assembly results (fasta.nucl) for sample %s.\n" % (self.get_step_name(), sample)
        pass
        
        

    def build_scripts(self):
        """ This is the actual script building function
            Most, if not all, editing should be done here 
            HOWEVER, DON'T FORGET TO CHANGE THE CLASS NAME AND THE FILENAME!
        """

        
        # Each iteration must define the following class variables:
            # spec_script_name
            # script
        if "generate_GFF_dir" in self.params.keys():
            #Make a dir for the GFF files:
            GFF_dir = self.make_folder_for_sample("GFF")
            self.sample_data["GFF_dir"]=GFF_dir
            
            
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
            output_filename = "".join([use_dir , sample ])

            self.script += self.get_script_const()
            if "--proteins VFDB" in self.script:
                if "Virulence_Resistance.fasta" in os.listdir(self.module_location):
                    self.script=self.script.replace("--proteins VFDB","--proteins %s" % os.path.join(self.module_location,"Virulence_Resistance.fasta") )
                else:
                    raise AssertionExcept("The file %s is not found in the Prokka module directory" % "Virulence_Resistance.fasta" )


            self.script += "--outdir %s \\\n\t"   % use_dir
            self.script += "--locustag %s \\\n\t" % sample
            self.script += "--strain %s \\\n\t"   % sample
            self.script += "--prefix %s \\\n\t"   % sample
            self.script += "%s \n\n" % self.sample_data[sample]["fasta.nucl"]
            if "generate_GFF_dir" in self.params.keys():
                self.script += "cp %s  %%s \n\n" % os.path.join(sample_dir,sample+".gff") % GFF_dir
            
            
            # Store Prokka result files:
            #set_Sample_data(self,sample,["GFF"],os.path.join(sample_dir,sample+".gff"))
            self.sample_data[sample]["GFF"]=os.path.join(sample_dir,sample+".gff")
            self.sample_data[sample]["fasta.nucl"]=os.path.join(sample_dir,sample+".ffn")
            self.sample_data[sample]["fasta.prot"]=os.path.join(sample_dir,sample+".faa")
            # Wrapping up function. Leave these lines at the end of every iteration:
            self.local_finish(use_dir,sample_dir)       # Sees to copying local files to final destination (and other stuff)
            
            self.create_low_level_script()
