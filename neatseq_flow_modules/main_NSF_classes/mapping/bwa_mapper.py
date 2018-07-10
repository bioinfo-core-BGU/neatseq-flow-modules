# -*- coding: UTF-8 -*-
""" 
``bwa_mapper`` :sup:`*`
-----------------------------------------------------------------
:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A module for running bwa mapper:

The reads stored in each sample are aligned to one of the following bwa indices:

* An external index passed with the ``ref_index`` parameter.
* A bwa index on a project fasta files, such as an assembly from all samples. Specify with ``bwa_mapper:scope  project``
* A sample bwa index on a sample-specific fasta file, such as from a sample-wise assembly or from the sample fasta file. Specify with ``bwa_mapper:scope  sample``

The latter two options **must** come after a ``bwa_builder`` instance.

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* fastq files in one of the following slots:

    * ``sample_data[<sample>]["fastq.F"]``
    * ``sample_data[<sample>]["fastq.R"]``
    * ``sample_data[<sample>]["fastq.S"]``
* If ``mod`` is one of ``samse, sampe``, the sai files are required as well (created by a ``bwa aln`` step:
    * ``self.sample_data[<sample>]["saiF|saiR|saiS"]``
    

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Puts output sam files in the following slots:
    * If ``mod`` is one of ``mem, samse, sampe, bwasw``:
        * ``self.sample_data[<sample>]["sam"]``
    * If ``mod`` is ``aln``:
        * ``self.sample_data[<sample>]["saiF|saiR|saiS"]``

* Puts the name of the mapper in:
    * ``self.sample_data[<sample>]["mapper"]``

* puts fasta of reference genome (if one is given in param file) in:
    * ``self.sample_data[<sample>]["reference"]``

Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "ref_index", "path to bwa index", "If not given, will look for a project bwa index and then for a sample bwa index"
    "ref_genome", "path to genome fasta", "If ref_index is NOT given, will use the equivalent internal fasta. If ref_index is passed, and ref_genome is NOT passed, will leave the reference slot empty"
    "scope", "project | sample", "Indicates whether to use a project or sample bwa index."

Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**For external index:**

1. Using ``mem``:

::

    bwa_mem_1:
        module: bwa_mapper
        base: trim1
        script_path: /path/to/bwa
        mod: mem
        qsub_params:
            -pe: shared 20
        ref_genome: /path/to/ref_genome.fna
        ref_index: /path/to/bwa_index/ref_genome
        redirects:
            -t: 20

**2. Using ``aln - samse/sampe``:**

::

    bwa_aln_1:
        module: bwa_mapper
        base: trim1
        script_path: /path/to/bwa_mapper
        mod: aln
        qsub_params:
            -pe: shared 20
        ref_genome: /path/to/ref_genome.fna
        ref_index: /path/to/bwa_index/ref_genome
        redirects:
            -t: 20
    bwa_samse_1:
        module: bwa_mapper
        base: bwt2_1
        script_path: /path/to/bwa
        mod: samse
        ref_genome: /path/to/ref_genome.fna
        ref_index: /path/to/bwa_index/ref_genome
        

**For project bwa index:**

::

    bwa_1:
        module: bwa_mapper
        base: bwa_bld_ind
        script_path: /path/to/bwa
        mod: mem
        scope: project

References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Li, H. and Durbin, R., 2009. **Fast and accurate short read alignment with Burrows–Wheeler transform**. *Bioinformatics*, 25(14), pp.1754-1760.
"""


import os
import sys, re
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.1.0"


class Step_bwa_mapper(Step):
    
   
    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        # self.file_tag = "Bowtie_mapper"
        
        # Check if you can split the script_path (by space) into path + mod:
        try:
            # try splitting script_path by space and extracting the mod from the snd element:
            mod = re.split("\s+", self.params["script_path"])[1]
        except IndexError:
            mod = None
        
        # Reduce script_path to path only. The mod is treated separately.
        self.params["script_path"] = re.split("\s+", self.params["script_path"])[0]
        
        
        # Make sure mod is defined only once, and if passed through script_path, add to params.
        if "mod" in self.params:
            if(mod):
                raise AssertionExcept("You supplied mod as parameter as well as in script path.")
        else:
            if(mod):
                self.params["mod"] = mod
            else:
                raise AssertionExcept("You must supply a 'mod' parameter or add the mod to the end of the script path.\n\te.g. /path/to/bwa mem")
        

    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        
        
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash

            if "sam" in self.sample_data[sample]:
                self.write_warning("SAM file exists for sample. Double mapping steps?\n", sample)
           
           
            if self.params["mod"] in ["samse"]:
                try:
                    self.sample_data[sample]["saiS"]
                    self.sample_data[sample]["fastq.S"]
                except KeyError:
                    raise AssertionExcept("'samse' requires sai and single-end fatsq files for the sample. Make sure you have a bwa aln step before this step and 'Single' files in the sample file.", sample)
            if self.params["mod"] in ["sampe"]:
                try:
                    self.sample_data[sample]["saiF"]
                    self.sample_data[sample]["saiR"]
                    self.sample_data[sample]["fastq.F"]
                    self.sample_data[sample]["fastq.R"]

                except KeyError:
                    raise AssertionExcept("'sampe' requires sai and paired-end fatsq files for the sample. Make sure you have a bwa aln step before this step and 'Forward' and 'Reverse' files in the sample file.", sample)
                    
            
 
        
        # Require either 'scope' or 'ref_index':
        if "scope" in self.params:
            # If scope defined, comment if also ref_index exists.
            if "ref_index" in self.params:
                raise AssertionExcept("Both 'scope' and 'ref_index' specified!\n")

            try:
                # Loop over samples to set the reference genome:
                for sample in self.sample_data["samples"]:
                    if self.params["scope"] == "project":
                        # Set project wide reference:
                        self.sample_data[sample]["reference"] = self.sample_data["bwa_fasta"]
                    elif self.params["scope"] == "sample":
                        # Set per-sample reference:
                        self.sample_data[sample]["reference"] = self.sample_data[sample]["bwa_fasta"]
                    else:
                        raise AssertionExcept("Scope must be either 'sample' or 'project'")
                
            except KeyError:
                
                raise AssertionExcept("There is a mismatch between 'scope' and the existing bwa index\n", sample)
                
            if "ref_genome" in self.params.keys():
                raise AssertionExcept("ref_genome was passed, and 'scope' was defined. Ignoring ref_genome\n")
        else:
            # If scope is not defined, require '-x'
            if not "ref_index" in self.params:
                raise AssertionExcept("Neither 'scope' nor 'ref_index' specified.\n")
            # Storing reference genome for use by downstream steps:
            if "ref_genome" in self.params.keys():
                for sample in self.sample_data["samples"]:
                    # If reference already exists, ignore ref_genome
                    if "reference" in self.sample_data[sample]:
                        self.write_warning("ref_genome was passed, but a reference already exists. Setting reference to 'ref_genome'\n")
                        
                
                    self.sample_data[sample]["reference"] = self.params["ref_genome"]
            else:
                self.write_warning("No reference given. It is highly recommended to give one!\n")

        
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

            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)
 
            
            if self.params["mod"] in ["aln"]:
                
                for direction in filter(lambda x: x in ["fastq.F","fastq.R","fastq.S"], self.sample_data[sample].keys()):

                    self.script = ""
                    direction_tag = direction[-1] # Get last letter in direction
                    # Name of specific script:
                    self.spec_script_name = self.jid_name_sep.join([self.step,self.name,sample,direction_tag])
                    self.script = ""

                    

                    output_filename = "%s.%s.bwa.sai" % (sample, direction)

                    
                    # Get constant part of script:
                    self.script = self.get_script_env_path()
                    # Add mod:
                    self.script += "%s \\\n\t" % self.params["mod"]
                    # Add redir_params:
                    self.script += self.get_redir_parameters_script()

                    
                    # Add ref_index (depends on scope)
                    if "scope" in self.params:  # If scope was passed, include either project or sample bwa index
                        if self.params["scope"] == "project":
                            self.script += "%s \\\n\t" % self.sample_data["bwa_index"]
                        else:
                            self.script += "%s \\\n\t" % self.sample_data[sample]["bwa_index"]
                    else:  # Otherwise add ref_index
                        self.script += "%s \\\n\t" % self.params["ref_index"]

                    # Add reads
                    self.script += "%s \\\n\t" % self.sample_data[sample][direction]

                    # Add output:
                    self.script += "> %s\n\n" % (use_dir + output_filename)

                    file_code = "sai%s" % direction_tag
                    self.sample_data[sample][file_code] = (sample_dir + output_filename)
                    self.stamp_file(self.sample_data[sample][file_code])
                    
                    # Storing name of mapper. might be useful:
                    self.sample_data[sample]["mapper"] = self.get_step_step()  
                    
                    # Move all files from temporary local dir to permanent base_dir
                    self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
                    
                    self.create_low_level_script()
                           
            else:  # Not 'aln': one of the mods that create sam files.
                
                # Name of specific script:
                self.spec_script_name = self.set_spec_script_name(sample)
                self.script = ""
                
                
                # Define location of output file:
                output_filename = "%s.bwa.sam" % (sample)

                # Get constant part of script:
                self.script += self.get_script_env_path()
                # Add mod:
                self.script += "%s \\\n\t" % self.params["mod"]
                # Add redir_params:
                self.script += self.get_redir_parameters_script()

                # Deal with read group title:
                if self.params["mod"] in ["mem"]:
                    self.script += "-R '@RG\\tID:%(sample)s\\tSM:%(sample)s' \\\n\t" % {"sample":sample}
                elif self.params["mod"] in ["sampe","samse"]:
                    self.script += "-r '@RG\\tID:%(sample)s\\tSM:%(sample)s' \\\n\t" % {"sample":sample}
                else:
                    pass
                
                # Add ref_index (depends on scope)
                if "scope" in self.params:  # If scope was passed, include either project or sample bwa index
                    if self.params["scope"] == "project":
                        self.script += "%s \\\n\t" % self.sample_data["bwa_index"]
                    else:
                        self.script += "%s \\\n\t" % self.sample_data[sample]["bwa_index"]
                else:  # Otherwise add ref_index
                    self.script += "%s \\\n\t" % self.params["ref_index"]

                # Add reads
                if self.params["mod"] in ["mem"]:
                    if "fastq.F" in self.sample_data[sample].keys():
                        self.script += "%s \\\n\t%s\\\n\t" % \
                            (self.sample_data[sample]["fastq.F"],
                            self.sample_data[sample]["fastq.R"])
                    elif "fastq.S" in self.sample_data[sample].keys():
                        self.script += "%s \\\n\t" % self.sample_data[sample]["fastq.S"]
                    else:
                        pass
                    if "fastq.F" in self.sample_data[sample] and "fastq.S" in self.sample_data[sample]:
                        self.write_warning("Both paired- and single-end sequence files exists for sample. Using only paired data\n", sample)
                    
                    
                    
                if self.params["mod"] in ["samse"]:
                    self.script += "%s \\\n\t%s\\\n\t" % \
                        (self.sample_data[sample]["saiS"],\
                         self.sample_data[sample]["fastq.S"])
                if self.params["mod"] in ["sampe"]:
                    self.script += "%s \\\n\t%s\\\n\t%s \\\n\t%s\\\n\t" % \
                        (self.sample_data[sample]["saiF"], \
                        self.sample_data[sample]["saiR"], \
                        self.sample_data[sample]["fastq.F"],
                        self.sample_data[sample]["fastq.R"])

                # Add output:
                self.script += "> %s\n\n" % (use_dir + output_filename)

                
                
                self.sample_data[sample]["sam"] = (sample_dir + output_filename)
                self.stamp_file(self.sample_data[sample]["sam"])
                
                # Storing name of mapper. might be useful:
                self.sample_data[sample]["mapper"] = self.get_step_step()  
                

       
                # Move all files from temporary local dir to permanent base_dir
                self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
                
                self.create_low_level_script()
                        
