#!/fastspace/bioinfo_apps/python-2.7_SL6/bin/python

# -*- coding: UTF-8 -*-
""" 
``Snippy``
-----------------------

:Authors: Liron Levin
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

.. Note:: This module was developed as part of a study led by Dr. Jacob Moran Gilad

Short Description
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    A module for running Snippy on fastq files

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    * fastq files in at least one of the following slots:
        ``sample_data[<sample>]["fastq.F"]``
        ``sample_data[<sample>]["fastq.R"]``
        ``sample_data[<sample>]["fastq.S"]``

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    * puts Results directory location in:
        ``self.sample_data[<sample>]["Snippy"]``
    * puts for each sample the vcf file location in:
        ``self.sample_data[<sample>]["vcf"]``
        
    if snippy_core is set to run:
        * puts the core Multi-FASTA alignment location in:
            ``self.sample_data["fasta.nucl"]``
        * puts core vcf file location of all analyzed samples in the following slot:
            ``self.sample_data["vcf"]``
            
    if Gubbins is set to run: 
        * puts result Tree file location of all analyzed samples in:
            ``self.sample_data["newick"]``
        * update the core Multi-FASTA alignment in:
            ``self.sample_data["fasta.nucl"]``
        * update the core vcf file in the slot:
            ``self.sample_data["vcf"]``
            
    if pars is set to run, puts phyloviz ready to use files in:
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
        ``Snippy v3.2``
        ``gubbins v2.2.0``
    * For the pars analysis the following python packages are required:
        ``pandas``

Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    Step_Name:                                  # Name of this step
        module: Snippy                          # Name of the module used
        base:                                   # Name of the step [or list of names] to run after [must be after a merge step]
        script_path:                            # Command for running the Snippy script
        env:                                    # env parameters that needs to be in the PATH for running this module
        qsub_params:
            -pe:                                # Number of CPUs to reserve for this analysis
        gubbins:
            script_path:                        # Command for running the gubbins script, if empty or this line dose not exist will not run gubbins
            --STR:                              # More redirects arguments for running gubbins
        pars:                                   # Generate phyloviz ready to use files
            -M:                                 # Location of a MetaData file 
            --Cut:                              # Use only Samples found in the metadata file
            --S_MetaData:                       # The name of the samples ID column
            -C:                                 # Use only Samples that has at least this fraction of identified alleles
        snippy_core:
            script_path:                        # Command for running the snippy-core script, if empty or this line dose not exist will not run snippy-core
            --noref:                            # Exclude reference 
        redirects:
            --cpus:                             # Parameters for running Snippy
            --force:                            # Force overwrite of existing output folder (default OFF)
            --mapqual:                          # Minimum mapping quality to allow
            --mincov:                           # Minimum coverage of variant site
            --minfrac:                          # Minumum proportion for variant evidence
            --reference:                        # Reference Genome location
            --cleanup                           # Remove all non-SNP files: BAMs, indices etc (default OFF)            

References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Snippy:
        https://github.com/tseemann/snippy
    gubbins:
        Croucher N. J., Page A. J., Connor T. R., Delaney A. J., Keane J. A., Bentley S. D., Parkhill J., Harris S.R. "Rapid phylogenetic analysis of large samples of recombinant bacterial whole genome sequences using Gubbins". doi:10.1093/nar/gku1196, Nucleic Acids Research, 2014
"""



import os
import sys
import re
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Levin Liron"
__version__= "1.2.0"

class Step_Snippy(Step):

    def step_specific_init(self):
        self.shell = "bash"
        self.file_tag = ""
        import inspect
        self.module_location=os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))

        assert "--outdir" not in self.params["redir_params"] or len({"--pe1", "--pe2", "--R1", "--R2","--left", "--right", "--se", "--single"}  & set(self.params["redir_params"]))==0, \
            "you should not give output directory in step %s\n" % self.get_step_name()

        assert "--reference" in self.params["redir_params"] , \
            "you should give reference file [Supports FASTA, GenBank, EMBL (not GFF)] in step %s\n" % self.get_step_name()


    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        
        # Assert that all samples have reads files:
        for sample in self.sample_data["samples"]:    
            assert {"fastq.F", "fastq.R", "fastq.S"} & set(self.sample_data[sample].keys()), "Sample %s does not have read files in step %s\n" % (sample, self.name)
        
        pass
        
    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        prefix="core"
        # Make a core snp fasta file of all results:
        if "snippy_core" in self.params.keys():
            if get_dic_data(self.params,["snippy_core","script_path"])!=None:
                if self.params["snippy_core"]["script_path"]!=None:
                    # Make a dir for the results file:
                    results_dir = self.make_folder_for_sample("Results")         
                    #Running the file merge script
                    self.script = ""
                    self.script +="cd '%s' \n\n" % results_dir
                    if "env" in self.params.keys():
                        if self.shell=="bash":
                            self.script +="export env %s  \\\n\t" % self.params["env"]
                        else:
                            self.script +="env %s  \\\n\t" % self.params["env"]
                    self.script +="%s  \\\n\t" % self.params["snippy_core"]["script_path"]
                    for par in self.params["snippy_core"].keys():
                        if par !="script_path":
                            if len(par)>0:
                                if self.params["snippy_core"][par]!=None:
                                    self.script +="%s  %%s \\\n\t" % par \
                                                                   % self.params["snippy_core"][par]
                                else:
                                    self.script +="%s  \\\n\t" % par 
                            if par=="--prefix":
                                prefix=self.params["snippy_core"][par]
                    for sample in self.sample_data["samples"]:
                        self.script +="%s \\\n\t" % self.sample_data[sample]["Snippy"]
                    self.script +="\n\n"
                    self.sample_data["fasta.nucl"]=os.sep.join([results_dir.rstrip(os.sep),prefix+".full.aln"])
                    self.sample_data["vcf"]=os.sep.join([results_dir.rstrip(os.sep),prefix+".vcf"])
                    

                    #Run gubbins if gubbins_script_path is given
                    if "gubbins" in self.params.keys():
                        if get_dic_data(self.params,["gubbins","script_path"])!=None:
                            if self.params["gubbins"]["script_path"]!=None:
                                # Make a dir for the results file:
                                Gubbins_results_dir = self.make_folder_for_sample("Gubbins")         
                                #Running the file merge script
                                self.script +="cd '%s' \n\n" % Gubbins_results_dir
                                if "env" in self.params.keys():
                                    if self.shell=="bash":
                                        self.script +="export env %s  \\\n\t" % self.params["env"]
                                    else:
                                        self.script +="env %s  \\\n\t" % self.params["env"]
                                self.script +="%s  \\\n\t" % self.params["gubbins"]["script_path"]
                                for par in self.params["gubbins"].keys():
                                    if par !="script_path":
                                        if len(par)>0:
                                            if self.params["gubbins"][par]!=None:
                                                self.script +="%s  %%s \\\n\t" % par \
                                                                               % self.params["gubbins"][par]
                                            else:
                                                self.script +="%s  \\\n\t" % par 
                                self.script +="%s  \n\n" % self.sample_data["fasta.nucl"]
                                self.sample_data["vcf"]=os.sep.join([Gubbins_results_dir.rstrip(os.sep),prefix+".full.summary_of_snp_distribution.vcf"])
                                self.sample_data["fasta.nucl"]=os.sep.join([Gubbins_results_dir.rstrip(os.sep),prefix+".full.filtered_polymorphic_sites.fasta"])
                                self.sample_data["newick"]=os.sep.join([Gubbins_results_dir.rstrip(os.sep),prefix+".full.final_tree.tre"])
                                


                    if "pars" in self.params.keys(): 
                        if "Snippy_parser.py" in os.listdir(self.module_location):
                            # Make a dir for the parsed files:
                            pars_dir = self.make_folder_for_sample("pars")
                            if "env" in self.params.keys():
                                if self.shell=="bash":
                                    self.script +="export env %s  \\\n\t" % self.params["env"]
                                else:
                                    self.script +="env %s  \\\n\t" % self.params["env"]
                            self.script +="python %s  \\\n\t" % os.path.join(self.module_location,"Snippy_parser.py")
                            if is_it_dict(self.params["pars"]):
                                for par in self.params["pars"].keys():
                                    if len(par)>0:
                                        if self.params["pars"][par]!=None:
                                            self.script +="%s  %%s \\\n\t" % par \
                                                                           % self.params["pars"][par]
                                        else:
                                            self.script +="%s  \\\n\t" % par
                            self.script += " -F %s \\\n\t" %  self.sample_data["fasta.nucl"]
                            self.script += " --FASTA  \\\n\t"
                            self.script += " -O %s \n\n" % pars_dir
                            self.sample_data["phyloviz_Alleles"]= os.path.join(pars_dir,"phyloviz_Alleles.tab")  
                            self.sample_data["phyloviz_MetaData"]= os.path.join(pars_dir,"phyloviz_MetaData.tab")  
                        else:
                            raise AssertionExcept("The file %s is not found in the Snippy module directory" % "Snippy_parser.py" )
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
            
            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)

            # Name of specific script:
            self.spec_script_name = "_".join([self.step,self.name,sample])
            self.script = ""
            
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)

            # Get constant part of script:
            self.script += self.get_script_const()
            #fastq check if it is a paired-end add the fastq files
            if len({"fastq.F", "fastq.R"} & set(self.sample_data[sample].keys()))==2:
                self.script +="--pe1 %s  \\\n\t" % self.sample_data[sample]["fastq.F"]
                self.script +="--pe2 %s  \\\n\t" % self.sample_data[sample]["fastq.R"]
            else:
                self.script +="--se  %s  \\\n\t" % self.sample_data[sample]["fastq.S"]
            self.script +="--prefix  %s  \\\n\t" % sample
            self.script +="--outdir  %s  \n\n" % use_dir
            self.sample_data[sample]["Snippy"]=sample_dir
            self.sample_data[sample]["vcf"]=os.path.join(sample_dir,sample+".vcf") 
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