#!/fastspace/bioinfo_apps/python-2.7_SL6/bin/python
# -*- coding: UTF-8 -*-
""" 
``RSEM``
-----------------------

:Authors: Liron Levin
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.


Short Description
~~~~~~~~~~~~~~~~~~~~
    A module for running RSEM

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    * fastq file in 
        ``self.sample_data[sample]["fastq.F"]``
        ``self.sample_data[sample]["fastq.R"]``
        ``self.sample_data[sample]["fastq.S"]``
    * or bam file in 
        ``self.sample_data[sample]["bam"]``

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    * puts output bam files (if the input is fastq) in:
        ``self.sample_data[sample]["bam"]``
        ``self.sample_data[sample]["unsorted_bam"]``
    * puts the location of RSEM results in:
        ``self.sample_data[sample]["RSEM"]``

Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "mode",  "transcriptome/genome ", "Is the reference is a genome or a transcriptome?"
    "gff3","None","Use if the mode is genome and the annotation file is in gff3 format"
    "del_unsorted_bam","None","Delete unsorted bam files to save space"

Comments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    *  This module was tested on:
        ``RSEM v1.2.25``
        ``bowtie2 v2.2.6``

Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    Step_Name:                                                   # Name of this step
        module: RSEM                                             # Name of the module used
        base:                                                    # Name of the step [or list of names] to run after [must be after a bam file generator step or merge with fastq files]
        script_path:                                             # Command for running the RSEM script 
        qsub_params:
            -pe:                                                 # Number of CPUs to reserve for this analysis
        mode:                                                    # transcriptome or genome
        annotation:                                              # For Genome mode: the location of GTF file [the default] , for GFF3 use the gff3 flag. For Transcriptome mode: transcript-to-gene-map file.
                                                                 # If annotation is set to Trinity the transcript-to-gene-map file will be generated using the from_Trinity_to_gene_map script
        from_Trinity_to_gene_map_script_path:                    # If the mode is transcriptome and the reference was assembled using Trinity it is possible to generate the transcript-to-gene-map file automatically using this script
                                                                 # If annotation is set to Trinity and this line is empty or missing it will try using the module's associated script
        gff3:                                                    # Use if the mode is genome and the annotation file is in gff3 format
        mapper:                                                  # bowtie/bowtie2/star 
        mapper_path:                                             # Location of mapper script
        rsem_prepare_reference_script_path:                      # Location of preparing reference script
        plot_stat:                                               # Generate statistical plots
        plot_stat_script_path:                                   # Location of statistical plot generating script
        reference:                                               # The reference genome/transcriptome location [FASTA file]
        rsem_generate_data_matrix_script_path:                   # Location of the final matrix generating script
                                                                 # If this line is empty or missing it will try using the module's associated script
        del_unsorted_bam:                                        # Delete unsorted bam files to save space.
        redirects:
            --append-names:                                      # RSEM will append gene_name/transcript_name to the result files
            --estimate-rspd:                                     # Enables RSEM to learn from the data how the reads are distributed across a transcript
            -p:                                                  # Number of CPUs to use in this analysis
            --bam:                                               # Will use bam files and not fastq
            --no-bam-output:
            --output-genome-bam:                                 # Alignments in genomic coordinates (only if mode is genome)

References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Li, Bo, and Colin N. Dewey. "RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome." BMC bioinformatics 12.1 (2011): 323.‏

"""


import os
import sys
import re
from neatseq_flow.PLC_step import Step,AssertionExcept

__author__ = "Levin Levin"
__version__= "1.2.0"


class Step_RSEM(Step):
 
    def step_specific_init(self):
        """ Called on intiation
            Good place for parameter testing.
            Wrong place for sample data testing
        """
        self.shell = "bash"
        self.file_tag = ""
        assert "mode"  in self.params.keys() , \
            "you should  provide mode type [transcriptome or genome] in step %s\n" % self.get_step_name()
        assert "reference"  in self.params.keys() , \
            "you should  provide reference file in step %s\n" % self.get_step_name()
        assert "mapper"  in self.params.keys() , \
            "you should  provide mapper type [bowtie, bowtie2 or star] in step %s\n" % self.get_step_name()
        assert "mapper_path"  in self.params.keys() , \
            "you should  provide mapper script location in step %s\n" % self.get_step_name()
        assert not ("--output-genome-bam" in self.params["redir_params"].keys()) &("transcriptome" in self.params["mode"]) , \
            "you can't use '--output-genome-bam' option when the mode is 'transcriptome' in step  %s\n" % self.get_step_name()
        assert "rsem_prepare_reference_script_path"  in self.params.keys() , \
            "you should  provide rsem_prepare_reference script location in step %s\n" % self.get_step_name()
        if "plot_stat" in self.params.keys():
            assert "plot_stat_script_path"  in self.params.keys() , \
                "you should  provide plot_stat script location in step %s\n" % self.get_step_name()
                
        import inspect
        self.module_location=os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        
        if "--bam" not in self.params["redir_params"].keys():                
            # Assert that all samples have reads files:
            for sample in self.sample_data["samples"]:    
                assert {"fastq.F", "fastq.R", "fastq.S"} & set(self.sample_data[sample].keys()), "Sample %s does not have read files in step %s\n if you have bam files use --bam\n" % (sample, self.name)
            
        else:
            for sample in self.sample_data["samples"]:
                if "bam" not in self.sample_data[sample].keys():
                    sys.exit("No Mapping or bam file information!!! \n")
        pass
        
    def create_spec_preliminary_script(self):
        """ Add script to run BEFORE all other steps
        """
        #Preparing reference genome/transcriptome
        #Creating new folder for the reference files
        REF_dir = self.make_folder_for_sample("Reference")
        #initiating new script 
        self.script = ""
        if ("transcriptome" in self.params["mode"])&("Trinity" in self.params["annotation"]):
            if "from_Trinity_to_gene_map_script_path" not in self.params.keys():
                if "Create_map_from_Trinity.py" not in os.listdir(self.module_location):
                    sys.exit("you should provide from_Trinity_to_gene_map_script_path !!! \n")
                else:
                    self.params["from_Trinity_to_gene_map_script_path"]== "python  %s "  % os.path.join(self.module_location,"Create_map_from_Trinity.py")
                    
            if self.params["from_Trinity_to_gene_map_script_path"]!=None:
                #preparing a transcript_to_gene map file from the reference transcriptome file [if it was created by Trinity] 
                self.script +="%s  %%s  %%%%s \n\n" % self.params["from_Trinity_to_gene_map_script_path"] \
                                                       % self.params["reference"] \
                                                       % os.sep.join([REF_dir.rstrip(os.sep),"transcript_to_gene_map.map"])
                #update the annotation slot to the new  transcript_to_gene_map annotation file
                self.params["annotation"]=os.sep.join([REF_dir.rstrip(os.sep),"transcript_to_gene_map.map"])
                
        #The main part of generating the reference files        
        self.script +=self.params["rsem_prepare_reference_script_path"]+" \\\n\t"
        if ("transcriptome" in self.params["mode"]):
            #If the reference is a transcriptome use the transcript_to_gene_map annotation file
            self.script +="-transcript-to-gene-map %s \\\n\t" % self.params["annotation"]
        elif ("genome" in self.params["mode"]):
            if "gff3" not in self.params.keys(): 
                #If the reference is a genome use the gtf annotation file
                self.script +="--gtf %s \\\n\t" % self.params["annotation"]
            else:
                #If the reference is a genome and the --gff3 flag is set, use the gff3 annotation file
                self.script +="--gff3 %s \\\n\t" % self.params["annotation"]
        else:
            sys.exit("mode can only be transcriptome or genome !!! \n")
        if  "bowtie" == self.params["mapper"]:
            self.script +="--%s-path %%s  \\\n\t" % self.params["mapper"] \
                                                  % self.params["mapper_path"]

        else:
            self.script +="--%s --%%s-path %%%%s  \\\n\t" % self.params["mapper"] \
                                                          % self.params["mapper"] \
                                                          % self.params["mapper_path"]
        self.script +="%s  \\\n\t%%s \n\n" % self.params["reference"] \
                                           % os.sep.join([REF_dir.rstrip(os.sep),"REF"])
        #update the reference slot to the new reference folder location and the reference files prefix 
        self.params["reference"]= os.sep.join([REF_dir.rstrip(os.sep),"REF"])        
                
        pass
        
    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        # Make a merge file of all results:
        if "rsem_generate_data_matrix_script_path" not in self.params.keys():
            if "Merge_RSEM.py" in os.listdir(self.module_location):
                 self.params["rsem_generate_data_matrix_script_path"]="python  %s "  % os.path.join(self.module_location,"Merge_RSEM.py")
        elif self.params["rsem_generate_data_matrix_script_path"]==None:
            if "Merge_RSEM.py" in os.listdir(self.module_location):
                 self.params["rsem_generate_data_matrix_script_path"]="python  %s "  % os.path.join(self.module_location,"Merge_RSEM.py")
        
        if "rsem_generate_data_matrix_script_path" in self.params.keys():
            if self.params["rsem_generate_data_matrix_script_path"]!=None:
                # Make a dir for the results file:
                results_dir = self.make_folder_for_sample("Results")         
                #Running the file merge script
                self.script = ""
                for sample in self.sample_data["samples"]:
                    self.script +="cp  '%s'  '%%s' \n\n" % (self.sample_data[sample]["RSEM"]+'.genes.results') \
                                                       % results_dir
                for sample in self.sample_data["samples"]:
                    self.script +="cp  '%s'  '%%s' \n\n" % (self.sample_data[sample]["RSEM"]+'.isoforms.results') \
                                                         % results_dir
                self.script +="\n\n"
                self.script +="cd '%s' \n\n" % results_dir
                self.script +="%s  \\\n\t" % self.params["rsem_generate_data_matrix_script_path"]
                # for sample in self.sample_data["samples"]:
                    # self.script +="%s  \\\n\t" % (self.sample_data[sample]["RSEM"]+'.genes.results')
                self.script +="%s  \\\n\t" % '*.genes.results'
                self.script +="> %s \n\n" % os.sep.join([results_dir.rstrip(os.sep),"GeneMat.results"])
                self.script +="%s  \\\n\t" % self.params["rsem_generate_data_matrix_script_path"]
                # for sample in self.sample_data["samples"]:
                    # self.script +="%s  \\\n\t" % (self.sample_data[sample]["RSEM"]+'.isoforms.results')
                self.script +="%s  \\\n\t" % '*.isoforms.results'
                self.script +="> %s \n\n" % os.sep.join([results_dir.rstrip(os.sep),"IsoMat.results"])
        
        if "plot_stat" in self.params.keys():
            for sample in self.sample_data["samples"]:
                self.script +="%s '%%s'  '%%%%s' \n\n"  % (self.params["plot_stat_script_path"]) \
                                                        % (self.sample_data[sample]["RSEM"]) \
                                                        % (results_dir+sample+"_diagnostic.pdf")
        
        if "del_unsorted_bam" in self.params.keys():
            for sample in self.sample_data["samples"]:
                try:  # Does a unsorted_bam slot exist? 
                    self.sample_data[sample]["unsorted_bam"]
                except KeyError:  # If failed...
                        pass
                else:  #Delete unsorted bams                    
                    self.script +="rm -f %s  \n\n" % self.sample_data[sample]["unsorted_bam"]
                    
                    
                    
    def build_scripts(self):
        """ This is the actual script building function
            Most, if not all, editing should be done here 
            HOWEVER, DON'T FORGET TO CHANGE THE CLASS NAME AND THE FILENAME!
        """
        
        
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash

            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)

            # Name of specific script:
            self.spec_script_name = "_".join([self.step,self.name,sample])
            self.script = "("
            
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)

            # Define location and prefix for output files:
            output_prefix = sample + "_RSEM"
            
            # Get constant part of script:
            self.script += self.get_script_const()
            # Adding the mapper type and script location
            if  "bowtie" == self.params["mapper"]:
                self.script +="--%s-path %%s  \\\n\t" % self.params["mapper"] \
                                                      % self.params["mapper_path"]

            else:
                self.script +="--%s --%%s-path %%%%s  \\\n\t" % self.params["mapper"] \
                                                              % self.params["mapper"] \
                                                              % self.params["mapper_path"]
            #Check if to use bam or fastq files
            if "--bam" not in self.params["redir_params"].keys(): 
                #if fastq check if it is a paired-end
                if len({"fastq.F", "fastq.R"} & set(self.sample_data[sample].keys()))==2:
                    self.script +="--paired-end \\\n\t"
                #Add the fastq files
                for i in self.sample_data[sample].keys():
                    if i in ["fastq.F", "fastq.R", "fastq.S"]:
                        self.script +="%s \\\n\t" % self.sample_data[sample][i]
                #self.script +=" \\\n\t"
                #Append the new bam file location to the bam slot 
                if "--output-genome-bam" in self.params["redir_params"].keys():
                    #if the --output-genome-bam option is present use the genome sorted bam
                    self.sample_data[sample]["bam"]=os.sep.join([sample_dir.rstrip(os.sep),sample+".genome.sorted.bam"])
                    #remember the unsorted bam as well
                    self.sample_data[sample]["unsorted_bam"]=os.sep.join([sample_dir.rstrip(os.sep),sample+".genome.bam"])
                else:
                    # the default is the transcript sorted bam
                    self.sample_data[sample]["bam"]=os.sep.join([sample_dir.rstrip(os.sep),sample+".transcript.sorted.bam"])
                    #remember the unsorted bam as well
                    self.sample_data[sample]["unsorted_bam"]=os.sep.join([sample_dir.rstrip(os.sep),sample+".transcript.bam"])
            else:
                #Add the bam file
                self.script +="%s \\\n\t" % self.sample_data[sample]["bam"]
            #The output information at the end 
            self.script +="%s \\\n\t%%s \\\n\t " % self.params["reference"]  % os.sep.join([use_dir.rstrip(os.sep),sample])
            #Generate log file:
            self.script += "> %s.out ) >& %%s.log\n\n"  % os.sep.join([use_dir.rstrip(os.sep),output_prefix]) \
                                                        % os.sep.join([use_dir.rstrip(os.sep),output_prefix])
            #Append the location of RSEM results
            self.sample_data[sample]["RSEM"]=os.sep.join([sample_dir.rstrip(os.sep),sample])

            # Wrapping up function. Leave these lines at the end of every iteration:
            self.local_finish(use_dir,sample_dir)       # Sees to copying local files to final destination (and other stuff)         
            self.add_jid_to_jid_list()
            self.create_low_level_script()

def set_Sample_data_dir(self,category,info,data):
    if category not in self.keys():
        self[category] = {}
    if info not in self[category].keys():
        self[category][info] = {}
    self[category][info] = data 