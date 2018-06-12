# -*- coding: UTF-8 -*-
""" 
``trinity`` :sup:`*`
-----------------------------------------------------------------
:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A class that defines a module for RNA_seq assembly using the `Trinity assembler`_.


.. Attention:: This module was tested on release 2.5.x. It should also work with 2.4.x 
    
    For old versions of Trinity, you might need to use ``trinity`` module.
    
    The main difference between the modules is that ``trinity_new`` creates an output directory with the word `trinity` in it as required by the newer release of Trinity.
    
    In order to run on the cluster, you need to install `HpcGridRunner`_.     
    
.. _Trinity assembler: https://github.com/trinityrnaseq/trinityrnaseq/wiki
.. _HpcGridRunner: http://hpcgridrunner.github.io/


Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    * ``fastq`` files in at least one of the following slots:
        
        * ``sample_data[<sample>]["fastq.F"]``
        * ``sample_data[<sample>]["fastq.R"]``
        * ``sample_data[<sample>]["fastq.S"]``

    
Output:
~~~~~~~~~~~~~

    * puts ``fasta`` output files in the following slots:
        
        * for sample-wise assembly:
        
            * ``sample_data[<sample>]["fasta.nucl"]``
            * ``sample_data[<sample>]["Trinity.contigs"]``
        
        * for project-wise assembly:
        
            * ``sample_data["fasta.nucl"]``
            * ``sample_data["Trinity.contigs"]``

                
Parameters that can be set        
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"

    "scope", "sample|project", "Set if project-wide fasta slot should be used"
    
    
Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    trinity1:
        module:                 trinity
        base:                   trin_tags1
        script_path:            {Vars.paths.Trinity}
        qsub_params:
            node:               sge213
            -pe:                shared 20
        redirects:
            --grid_exec:        "{Vars.paths.hpc_cmds_GridRunner} --grid_conf {Vars.paths.SGE_Trinity_conf} -c" 
            --grid_node_CPU:    40 
            --grid_node_max_memory: 80G 
            --max_memory:        80G 
            --seqType:          fq
            --min_kmer_cov:     2
            --full_cleanup:

References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Grabherr, M.G., Haas, B.J., Yassour, M., Levin, J.Z., Thompson, D.A., Amit, I., Adiconis, X., Fan, L., Raychowdhury, R., Zeng, Q. and Chen, Z., 2011. **Trinity: reconstructing a full-length transcriptome without a genome from RNA-Seq data**. *Nature biotechnology*, 29(7), p.644.

"""



import os
import sys
import re
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.1.0"


class Step_trinity_new(Step):
    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        self.file_tag = ".Trinity.fasta"
        
        
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
            Here you should do testing for dependency output. These will NOT exist at initiation of this instance. They are set only following sample_data updating
        """
        
        
        # Assert that all samples have reads files:
        for sample in self.sample_data["samples"]:    
            if not {"fastq.F", "fastq.R", "fastq.S"} & set(self.sample_data[sample].keys()):
                raise AssertionExcept("No read files\n",sample)
         
        if "scope" in self.params:
          
            if self.params["scope"]=="project":
                pass

            elif self.params["scope"]=="sample":
                
                for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
                    pass
            else:
                raise AssertionExcept("'scope' must be either 'sample' or 'project'")
        else:
            raise AssertionExcept("No 'scope' specified.")
        
        
        ##########################
        

            
        pass     

        
    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        
        pass

    def build_scripts(self):
    
        if self.params["scope"] == "project":
            self.build_scripts_project()
        else:
            self.build_scripts_sample()
            
            
    def build_scripts_project(self):
        
        
        # Name of specific script:
        self.spec_script_name = "_".join([self.step,self.name,self.sample_data["Title"]])

        self.script = ""

        # This line should be left before every new script. It sees to local issues.
        # Use the dir it returns as the base_dir for this step.
        use_dir = self.local_start(self.base_dir)

        forward = list()   # List of all forward files
        reverse = list()   # List of all reverse files
        single = list()    # List of all single files
        
        # Loop over samples and concatenate read files to $forward and $reverse respectively
        # add cheack if paiered or single !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
            # If both F and R reads exist, adding them to forward and reverse
            # Assuming upstream input testing to check that if there are F reads then there are also R reads.
            if "fastq.F" in self.sample_data[sample].keys():
                forward.append(self.sample_data[sample]["fastq.F"])
                reverse.append(self.sample_data[sample]["fastq.R"])
            if "fastq.S" in self.sample_data[sample].keys():
                single.append(self.sample_data[sample]["fastq.S"])

        # Concatenate all filenames separated by commas:
        single  = ",".join(single)   if (len(single) > 0) else None
        forward = ",".join(forward)  if (len(forward) > 0) else None
        reverse = ",".join(reverse)  if (len(reverse) > 0) else None

        # Adding single reads to end of left (=forward) reads
        if single != None and forward != None:
            forward = ",".join([forward,single])

        self.script += self.get_script_const()

        # Adding 'trinity' to output dir:
        # "output directory must contain the word 'trinity' as a safety precaution, given that auto-deletion can take place."
        output_basename = "{title}_trinity".format(title = self.sample_data["Title"])
        # The results will be put in data/step_name/name/Title
        self.script += "--output %s \\\n\t" % os.path.join(use_dir, output_basename)
            
        #add if single or paired!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (forward): 
            self.script += "--left %s \\\n\t" % forward
            self.script += "--right %s \n\n" % reverse
        elif (single):
            self.script += "--single %s \n\n" % single
        
        # The assembly itself is stored in the prefix + .Trinity.fasta:
        output_basename += ".Trinity.fasta"

        # Store results to fasta and assembly slots:
        self.sample_data["fasta.nucl"] = os.path.join(self.base_dir, output_basename) 
        self.sample_data[self.get_step_step() + ".contigs"] = self.sample_data["fasta.nucl"]

        self.stamp_file(self.sample_data["fasta.nucl"])


        #######################################################################################
        ## 
        ## Step 2: Create gene_to_trans_map
        
        # cmd_text = self.get_script_env_path() 
        if "skip_gene_to_trans_map" not in self.params:
            cmd_text = """
    {script_path} \\
        {transcriptome} \\
        > {map} 
""".format(transcriptome = os.path.join(use_dir, output_basename) ,\
           script_path   = os.path.join(os.path.dirname(os.path.normpath(self.params["script_path"])) , \
                                          "util/support_scripts/get_Trinity_gene_to_trans_map.pl"),\
           map           = "{contigs}.gene_trans_map".format(contigs=os.path.join(use_dir, output_basename)))
           
            
            self.script += """
# Creating transcript to gene mapping file:

if [ -e {transcriptome} ]
then
	{cmd_text}
fi
""".format(transcriptome = os.path.join(use_dir, output_basename),
           cmd_text = cmd_text)
           
            
            self.sample_data["gene_trans_map"] = "{contigs}.gene_trans_map".format(contigs=os.path.join(self.base_dir, output_basename))
            self.stamp_file(self.sample_data["gene_trans_map"])
              
        # Move all files from temporary local dir to permanent base_dir
        self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
     
            
        
        
        self.create_low_level_script()
                    
#################################################
    def build_scripts_sample(self):
        
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash

        # Name of specific script:
            self.spec_script_name = "_".join([self.step,self.name,sample])
            self.script = ""


            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)

            self.script += self.get_script_const()
            # Adding 'trinity' to output dir:
            # "output directory must contain the word 'trinity' as a safety precaution, given that auto-deletion can take place."
            output_basename = "{title}_trinity".format(title = sample)

            self.script += "--output %s \\\n\t" % os.path.join(use_dir,output_basename)
            
            if "fastq.F" in self.sample_data[sample].keys():
                self.script += "--left %s \\\n\t" % self.sample_data[sample]["fastq.F"]
                self.script += "--right %s \\\n\t" % self.sample_data[sample]["fastq.R"]
            if "fastq.S" in self.sample_data[sample].keys():
                self.script += "--single %s \n\n" % self.sample_data[sample]["fastq.S"]

            # If there is an extra "\\\n\t" at the end of the script, remove it.
            self.script = self.script.rstrip("\\\n\t") + "\n\n"

            # The assembly itself is stored in the prefix + .Trinity.fasta:
            output_basename += ".Trinity.fasta"

            # Store results to fasta and assembly slots:
            self.sample_data[sample]["fasta.nucl"] = os.path.join(sample_dir, output_basename) 

            self.sample_data[sample][self.get_step_step() + ".contigs"] = self.sample_data[sample]["fasta.nucl"]

            self.stamp_file(self.sample_data[sample]["fasta.nucl"])


            #######################################################################################
            ## 
            ## Step 2: Create gene_to_trans_map
            
            # cmd_text = self.get_script_env_path() 
            if "skip_gene_to_trans_map" not in self.params:
                cmd_text = """
        {script_path} \\
            {transcriptome} \\
            > {map} 
""".format(transcriptome = os.path.join(use_dir,output_basename),\
           script_path   = os.sep.join([os.path.split(os.path.normpath(self.params["script_path"]))[0] , \
                                          "util/support_scripts/get_Trinity_gene_to_trans_map.pl"]),\
           map           = "{contigs}.gene_trans_map".format(contigs=os.path.join(use_dir, output_basename)))
           
            
                self.script += """
# Create summary of biom table for use in rarefaction later

if [ -e {transcriptome} ]
then
    {cmd_text}
fi
""".format(transcriptome = os.path.join(use_dir, output_basename),
           cmd_text = cmd_text)
           
                
                self.sample_data[sample]["gene_trans_map"] = "{contigs}.gene_trans_map".format(contigs=os.path.join(sample_dir, output_basename))
                self.stamp_file(self.sample_data[sample]["gene_trans_map"])
           
                                
            # Wrapping up function. Leave these lines at the end of every iteration:
            self.local_finish(use_dir,sample_dir)       # Sees to copying local files to final destination (and other stuff)

            self.create_low_level_script()
                        
            
            
                 
            
     