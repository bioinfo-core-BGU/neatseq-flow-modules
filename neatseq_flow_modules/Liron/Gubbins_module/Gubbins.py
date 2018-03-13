#!/fastspace/bioinfo_apps/python-2.7_SL6/bin/python

# -*- coding: UTF-8 -*-
""" 
``Gubbins``
-----------------------

:Authors: Liron Levin
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

.. Note:: This module was developed as part of a study led by Dr. Jacob Moran Gilad

Short Description
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    A module for running Gubbins on a project level nucleotide Multi-FASTA alignment file.

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    * Project level nucleotide Multi-FASTA alignment file in the following slot:
        ``sample_data["fasta.nucl"]``

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    
    * puts result Tree file location of all analyzed samples in the slot:
        ``self.sample_data["newick"]``
    * update the Multi-FASTA alignment in the slot:
        ``self.sample_data["fasta.nucl"]``
    * puts the filtered vcf file in the slot:
        ``self.sample_data["vcf"]``
            
    if pars is set to run, puts phyloviz ready to use files in the slots:
        * Alleles:
            ``self.sample_data["phyloviz_Alleles"]``
        * MetaData:
            ``self.sample_data["phyloviz_MetaData"]``

Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "",  "", ""
    

Comments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    *  This module was tested on:
        ``gubbins v2.2.0``
    * For the pars analysis the following python packages are required:
        ``pandas``

Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    Step_Name:                                  # Name of this step
        module: Gubbins                         # Name of the module used
        base:                                   # Name of the step [or list of names] to run after [must be after a step that generates a Project level nucleotide Multi-FASTA alignment]
        script_path:                            # Command for running the gubbins script, if empty or this line dose not exist will not run gubbins
        env:                                    # env parameters that needs to be in the PATH for running this module
        qsub_params:
            -pe:                                # Number of CPUs to reserve for this analysis
        phyloviz:                                   # Generate phyloviz ready to use files
            -M:                                 # Location of a MetaData file 
            --Cut:                              # Use only Samples found in the metadata file
            --S_MetaData:                       # The name of the samples ID column
            -C:                                 # Use only Samples that has at least this fraction of identified alleles
        redirects:
            --threads:                          # Parameters for running Gubbins
            

References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    gubbins:
        Croucher N. J., Page A. J., Connor T. R., Delaney A. J., Keane J. A., Bentley S. D., Parkhill J., Harris S.R. "Rapid phylogenetic analysis of large samples of recombinant bacterial whole genome sequences using Gubbins". doi:10.1093/nar/gku1196, Nucleic Acids Research, 2014
"""



import os
import sys
import re
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Levin Liron"
__version__= "2.2.0"

class Step_Gubbins(Step):

    def step_specific_init(self):
        self.shell = "bash"
        self.file_tag = ""
        import inspect
        self.module_location=os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))


    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        
        # Assert that there is no project level nucleotide FASTA file :
        if (self.params["script_path"]!=None)&(self.params["script_path"]!=''):
            assert {"fasta.nucl"} & set(self.sample_data.keys()), "No project level nucleotide FASTA file in step %s\n" % ( self.name)
        
        pass
        

    def build_scripts(self):
        """ This is the actual script building function
            Most, if not all, editing should be done here 
            HOWEVER, DON'T FORGET TO CHANGE THE CLASS NAME AND THE FILENAME!
        """
        
        # Name of specific script:
        self.spec_script_name = "_".join([self.step,self.name])
        self.script = ""
        sample_dir=self.base_dir
        # This line should be left before every new script. It sees to local issues.
        # Use the dir it returns as the base_dir for this step.
        use_dir = self.local_start(sample_dir)
                
        if (self.params["script_path"]!=None)&(self.params["script_path"]!=''):
            Gubbins_results_dir = use_dir         
            self.script +="cd '%s' \n\n" % Gubbins_results_dir 
            # Get constant part of script:
            self.script += self.get_script_const()
            self.script +="%s  \n\n" % self.sample_data["fasta.nucl"]
            prefix=re.sub('\.\w+$','',os.path.basename(self.sample_data["fasta.nucl"]))
            self.sample_data["vcf"]=os.sep.join([Gubbins_results_dir.rstrip(os.sep),prefix+".summary_of_snp_distribution.vcf"])
            self.sample_data["fasta.nucl"]=os.sep.join([Gubbins_results_dir.rstrip(os.sep),prefix+".filtered_polymorphic_sites.fasta"])
            self.sample_data["newick"]=os.sep.join([Gubbins_results_dir.rstrip(os.sep),prefix+".final_tree.tre"])
            


            if "phyloviz" in self.params.keys(): 
                if "MLST_parser.py" in os.listdir(self.module_location):
                    # Make a dir for the parsed files:
                    pars_dir = self.make_folder_for_sample("Data_for_Phyloviz")
                    if "env" in self.params.keys():
                        if self.shell=="bash":
                            self.script +="export env %s  \\\n\t" % self.params["env"]
                        else:
                            self.script +="env %s  \\\n\t" % self.params["env"]
                    self.script +="python %s  \\\n\t" % os.path.join(self.module_location,"MLST_parser.py")
                    if is_it_dict(self.params["phyloviz"]):
                        for par in self.params["phyloviz"].keys():
                            if len(par)>0:
                                if self.params["phyloviz"][par]!=None:
                                    self.script +="%s  %%s \\\n\t" % par \
                                                                   % self.params["phyloviz"][par]
                                else:
                                    self.script +="%s  \\\n\t" % par
                    self.script += " -F %s \\\n\t" %  self.sample_data["fasta.nucl"]
                    self.script += " --FASTA  \\\n\t"
                    self.script += " -O %s \n\n" % pars_dir
                    self.sample_data["phyloviz_Alleles"]= os.path.join(pars_dir,"phyloviz_Alleles.tab")  
                    self.sample_data["phyloviz_MetaData"]= os.path.join(pars_dir,"phyloviz_MetaData.tab")  
                else:
                    raise AssertionExcept("The file %s is not found in the Snippy module directory" % "MLST_parser.py" )

        
        
        
        # Move all files from temporary local dir to permanent base_dir
        self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)


        if "spec_dir" in self.params.keys():
            self.script += "cd " + self.pipe_data["home_dir"] + "\n\n";
        
                    
        self.add_jid_to_jid_list()
        self.create_low_level_script()
                    

def get_dic_data(dic,category,default=None):
    try:
        res=reduce(dict.get, category, dic)
        if res==None:
            return default
        else:
            return res
    except:
        return default

def is_it_dict(dic):
    try:
        dic.keys()
        return True
    except:
        return False