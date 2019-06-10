# -*- coding: UTF-8 -*-
""" 
``Merge_bins_contigs``
-----------------------

:Authors: Liron Levin
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.


Short Description
~~~~~~~~~~~~~~~~~~~~~
    Merge Bin's contigs using stretch of Ns as contigs separator 

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    * For each Sample, a bins directory with fasta format files in:
        ``self.sample_data[sample]["bins_dir"]``
    * For project level , a bins directory with fasta format files in:
        ``self.sample_data["project_data"]["bins_dir"]``
        

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    * For each Sample, puts the location of the Sample's bins directory in:
        ``self.sample_data[sample]["bins_dir"]``
    * For project level, puts the location of the bins directory in:    
        ``self.sample_data[project_data]["bins_dir"]``
    * if Merge_all_bins option exist, puts the Merged fasta file of all bins in:
        ``self.sample_data["project_data"]["fasta.nucl"]``
        ``self.sample_data[sample]["fasta.nucl"]``

Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "Merge_all_bins",  "", "If set will merge all bins to one fata file"

Comments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    Step_Name:                                  # Name of this step
        module: Merge_bins_contigs              # Name of the module to use
        base:                                   # Name of the step [or list of names] to run after [must be after a binning generating step]
        script_path:                            # 
        NumberofNs:                             # The number of Ns for contigs separator, default = 100
        Merge_all_bins:                         # If set will merge all bins to one fata file
        bins_dir:                               # In project scope it is possible to indicate the bins directory.
        qsub_params:
            -pe:                                # Number of CPUs to reserve for this analysis
 
References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"""
import os
import sys
import re
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Liron Levin"


class Step_Merge_bins_contigs(Step):
    
    
    def step_specific_init(self):
        """ Called on intiation
            Good place for parameter testing.
            Wrong place for sample data testing
        """
        self.shell = "bash"
        import inspect
        self.module_location=os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))
        if "scope" not in list(self.params.keys()):
            self.params["scope"] = 'sample'
            
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        
        if "project" in self.params["scope"]:   
            if  "bins_dir" not in list(self.params.keys()):
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
        if "project" in self.params["scope"]:
            self.build_scripts_byproject()
        else:
            self.build_scripts_bysample()
        
    def build_scripts_byproject(self):
        
        # Name of specific script:
        self.spec_script_name = self.set_spec_script_name()
        self.script = ""
        # This line should be left before every new script. It sees to local issues.
        # Use the dir it returns as the base_dir for this step.
        use_dir = self.local_start(self.base_dir)
        if  "bins_dir" in list(self.params.keys()):
            bins_dir = self.params["bins_dir"]
        else:
            bins_dir = self.sample_data["project_data"]["bins_dir"]
        
        self.script += "for fasta_file in $(ls -p " + bins_dir  + " | grep -v / ); do \n\n"
             
        # Define output filename 
        output_filename = "".join([self.sample_data["Title"],'_','${fasta_file%.*}' ,'.fa'])
        if 'NumberofNs' in list(self.params.keys()):
            NumberofNs = self.params['NumberofNs']
        else:
            NumberofNs=100
            
        self.script += ''' awk ' /^>/ { if (c++ == 0) {print ">'${fasta_file%.*}'" }
                                        else {printf "%s", "''' + 'N'*int(NumberofNs)+ '''"}
                                        next
                                      } 
                                 /^$/ {next} 
                                      {printf "%s", $0} 
                                 END {print ""}' '''
        
        self.script += "%s > %%s \n\n" % os.path.join( bins_dir , '$fasta_file' ) % os.path.join( use_dir , output_filename )
       
        self.script += "done \n\n"
        
        self.sample_data["project_data"]["bins_dir"] = self.base_dir
        
        if 'Merge_all_bins' in list(self.params.keys()):
            #Make a dir for the merged bins file:
            All_bins_fasta_dir = self.make_folder_for_sample("All_bins_fasta")
            self.script += "cat %s > %%s \n\n" % os.path.join( use_dir , '*.fa' ) % os.path.join( All_bins_fasta_dir , "".join([self.sample_data["Title"],'_','All_bins.fasta' ]) )
            self.sample_data["project_data"]["fasta.nucl"] = os.path.join( self.base_dir ,"All_bins_fasta" ,"".join([self.sample_data["Title"],'_','All_bins.fasta' ]) )
        
        # Wrapping up function. Leave these lines at the end of every iteration:
        self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
        
        
        
        self.create_low_level_script()
    
        
        
        
    def build_scripts_bysample(self):
        # Each iteration must define the following class variables:
            # spec_script_name
            # script
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
            
            # Name of specific script:
            self.spec_script_name = self.set_spec_script_name(sample)
            self.script = ""
            bins_dir = self.sample_data[sample]["bins_dir"]
            self.script +='for fasta_file in $(ls -p ' + bins_dir + ' | grep -v / ); do\n'

            
            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)
            use_dir = self.local_start(sample_dir)
            output_filename = "".join([sample,'_','${fasta_file%.*}' ,'.fa'])
            if 'NumberofNs' in list(self.params.keys()):
                NumberofNs = self.params['NumberofNs']
            else:
                NumberofNs=100
                
            self.script += ''' awk ' /^>/ { if (c++ == 0) {print ">'${fasta_file%.*}'" }
                                            else {printf "%s", "''' + 'N'*int(NumberofNs)+ '''"}
                                            next
                                          } 
                                     /^$/ {next} 
                                          {printf "%s", $0} 
                                     END {print ""}' '''
            
            self.script += "%s > %%s \n\n" % os.path.join( bins_dir , '$fasta_file' ) % os.path.join( use_dir , output_filename )
           
            self.script += "done \n\n"
            
            self.sample_data[sample]["bins_dir"] = sample_dir
            
            if 'Merge_all_bins' in list(self.params.keys()):
                #Make a dir for the merged bins file:
                All_bins_fasta_dir = self.make_folder_for_sample("All_bins_fasta")
                self.script += "cat %s > %%s \n\n" % os.path.join( use_dir , '*.fa' ) % os.path.join( All_bins_fasta_dir , "".join([self.sample_data["Title"],'_','All_bins.fasta' ]) )
                self.sample_data[sample]["fasta.nucl"] = os.path.join( sample_dir ,"All_bins_fasta" ,"".join([self.sample_data["Title"],'_','All_bins.fasta' ]) )

            self.create_low_level_script()
