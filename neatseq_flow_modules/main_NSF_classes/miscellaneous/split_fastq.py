# -*- coding: UTF-8 -*-
""" 
``split_fasta``
------------------------


:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.


A module for splitting `fasta` files into parts.

Convenient for parallelizing processes on the cluster. You can take a project wide fasta file (such as a transcriptome), split it into sub-fasta files, and run various processes on the sub-files.

The parts can then be combined with ``merge_table`` module, which can concatenate any type of file.

.. Attention:: This module is not defined on the sample scope, yet. It will only take a project wide fasta and split it into pieces, each one in a new ``subsample``. The original set of samples will be overridden for the rest of the branch. (However, you can get them back by using one of the upstream instances as first base.

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* A `fasta` file in one of the following slots (scope = "project"):

    * ``sample_data["project_data"]["fasta.nucl"]``
    * ``sample_data["project_data"]["fasta.prot"]``

    

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Puts output files in the following slots:
        
    * ``sample_data[<sample>]["fasta.nucl"]``
    * ``sample_data[<sample>]["fasta.prot"]``



Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: Parameters that can be set:
    :header: "Parameter", "Values", "Comments"

    "type", "nucl|prot", "The type of fasta file to split"
    "subsample_num", "", "Number of fragments"
    
Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    split_fasta1:
        module:         split_fasta
        base:           Trinity1
        script_path:    
        type:           nucl
        subsample_num:      4


"""



import os
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.1.0"

class Step_split_fastq(Step):
   
    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        # self.file_tag = "Bowtie_mapper"
        
        # if "type" not in self.params:
        #     raise AssertionExcept("Please supply a 'type' parameter: 'nucl' or 'prot'")
        # if self.params["type"] not in ["nucl","prot"]:
        #     raise AssertionExcept("'type' parameter must be 'nucl' or 'prot'")
        # if "scope" not in self.params:
        #     raise AssertionExcept("Please supply a 'scope' parameter: 'sample' or 'project'")
        # if self.params["scope"] not in ["sample","project"]:
        #         raise AssertionExcept("'scope' parameter must be 'sample' or 'project'")

        
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """

        # Check that scope exists and fasta.nucl or fasta.prot exist
        # Check that subsample_num exists
        # Check that type exists

        # self.params["type"] = "fasta.{type}".format(type=self.params["type"])
        if self.params["scope"]=="project":
            if not any(map(lambda x: x in self.sample_data["project_data"], ["fastq.F","fastq.R","fastq.S"])):
                raise AssertionExcept("Project does not contains fastq files")
        else:
            for sample in self.sample_data["samples"]:
                # if self.params["type"] not in self.sample_data[sample]:
                if not any(map(lambda x: x in self.sample_data[sample], ["fastq.F", "fastq.R", "fastq.S"])):
                    raise AssertionExcept("Sample does not contains fastq files", sample)

        
    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        pass
        
    
    def build_scripts(self):
        """ This is the actual script building function
            Most, if not all, editing should be done here 
            HOWEVER, DON'T FORGET TO CHANGE THE CLASS NAME AND THE FILENAME!
        """
        
        if self.params["scope"] == "project":

            # # Name of specific script:
            # self.spec_script_name = self.set_spec_script_name()
            # self.script = ""

            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(self.base_dir)

            direction_dict = {"fastq.F":"Forward","fastq.R":"Reverse","fastq.S":"Single"}
            for fastq_type in [fastq for fastq in ["fastq.F","fastq.R","fastq.S"] if fastq in self.sample_data["project_data"]]:
                # Name of specific script:
                self.spec_script_name = self.jid_name_sep.join([self.set_spec_script_name(),
                                                                direction_dict[fastq_type]])
                self.script = """
fastq={project_fastq}
numseqs=$(wc -l $fastq| cut -f1 -d" ")
awk -v nseqs="$numseqs" '
    BEGIN {{n_seq=1; 
            file_cnt=1;
            nseqs=nseqs/4;
            seqspfile=int(nseqs/{subsample_num}); 
            numbigfiles=nseqs%{subsample_num};
            file=sprintf("{use_dir}/subsample%04d.{fastq_type}.fastq",file_cnt); 
            seq_line_cnt=1;
            }}
    {{
    if('NR%4==1')  {{
            if((n_seq==seqspfile+2 && file_cnt<=numbigfiles) || 
                (n_seq==seqspfile+1 && file_cnt>numbigfiles)){{
                file_cnt++; 
                file=sprintf("{use_dir}/subsample%04d.{fastq_type}.fastq",file_cnt); 
                n_seq=1;
            }} 
            print > file; 
            n_seq++; 
            next;
            }} 
    {{ print > file; }}
    }}
' < $fastq
            
            
""".format(project_fastq=self.sample_data["project_data"][fastq_type],
            subsample_num = self.params["subsample_num"],
           fastq_type=direction_dict[fastq_type],
            use_dir = use_dir)

            try:
                sample_list = ["subsample{num:0=4}".format(num=num) for num in range(1,int(self.params["subsample_num"])+1)]
            except ValueError:
                raise AssertionExcept("'subsample_num' must be an integer")
            self.stash_sample_list(sample_list)

            # Creating data container for subsamples:
            for sample in self.sample_data["samples"]:
                self.sample_data[sample] = dict()
                self.sample_data[sample][self.params["type"]] = "{use_dir}{sample}.fa".format(use_dir=self.base_dir,
                                                                                              sample=sample)
                # Stamping the files takes a long time. Cancelling for the time being
                # self.stamp_file(self.sample_data[sample][self.params["type"]])

            # Wrapping up function. Leave these lines at the end of every iteration:
            self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)

            self.create_low_level_script()
                    
        else:  # self.params["scope"] == "sample"
        
            # Each iteration must define the following class variables:
                # spec_script_name
                # script
            new_sample_list = list()
            new_sample_dict = dict()
            
            for sample in self.sample_data["samples"]:


                # Make a dir for the current sample:
                sample_dir = self.make_folder_for_sample(sample)

                # This line should be left before every new script. It sees to local issues.
                # Use the dir it returns as the base_dir for this step.
                use_dir = self.local_start(sample_dir)

                direction_dict = {"fastq.F":"Forward", "fastq.R":"Reverse", "fastq.S":"Single"}
                for fastq_type in [fastq for fastq in ["fastq.F", "fastq.R", "fastq.S"]
                                   if fastq in self.sample_data[sample]]:
                    # Name of specific script:
                    # self.spec_script_name = self.jid_name_sep.join([self.step, self.name, sample, direction_dict[fastq_type]])
                    self.spec_script_name = self.jid_name_sep.join([self.set_spec_script_name(sample),
                                                                    direction_dict[fastq_type]])

                    self.script = """

fastq={sample_fastq}
numseqs=$(wc -l $fastq| cut -f1 -d" ")
awk -v nseqs="$numseqs" '
    BEGIN {{n_seq=1; 
            file_cnt=1;
            nseqs=nseqs/4;
            seqspfile=int(nseqs/{subsample_num}); 
            numbigfiles=nseqs%{subsample_num};
            file=sprintf("{use_dir}/{sample}.subsample%04d.{fastq_type}.fastq",file_cnt); 
            seq_line_cnt=1;
            }}
    {{
    if('NR%4==1')  {{
            if((n_seq==seqspfile+2 && file_cnt<=numbigfiles) || 
                (n_seq==seqspfile+1 && file_cnt>numbigfiles)){{
                file_cnt++; 
                file=sprintf("{use_dir}/{sample}.subsample%04d.{fastq_type}.fastq",file_cnt); 
                n_seq=1;
            }} 
            print > file; 
            n_seq++; 
            next;
            }} 
    {{ print > file; }}
    }}
' < $fastq
            
""".format(sample_fastq = self.sample_data[sample][fastq_type],
           sample = sample,
           subsample_num = self.params["subsample_num"],
           fastq_type=direction_dict[fastq_type],
           use_dir = use_dir)
    

                    # A list of this sample's subsamples
                    subsample_list = ["{sample}.subsample{num:0=4}".format(sample=sample, num=num)
                                      for num
                                      in range(1, self.params["subsample_num"] + 1)]

                    new_sample_list.extend(subsample_list)
                    # CONTINUE HERE
                    for subsample in subsample_list:
                        if subsample not in self.sample_data:
                            self.sample_data[subsample] = dict()
                        self.sample_data[subsample][fastq_type] = \
                            "{use_dir}{sample}.{subsample}.{fastq_type}.fastq".format(use_dir=sample_dir,
                                                                                      fastq_type=direction_dict[fastq_type],
                                                                                      sample=sample,
                                                                                      subsample=subsample)
                        # Storing origin of subsample in grouping dict:
                        self.sample_data[subsample]["grouping"] = dict()
                        self.sample_data[subsample]["grouping"]["source"] = sample

                        # Stamping file
                        self.stamp_file(self.sample_data[sample][fastq_type])

                    # Wrapping up function. Leave these lines at the end of every iteration:
                    self.local_finish(use_dir,sample_dir)
                    self.create_low_level_script()

            self.sample_data["samples"] = new_sample_list
            for sample in self.sample_data["samples"]:
                self.sample_data[sample]["type"] = self.determine_sample_types(sample,
                                                                                  self.sample_data[sample])