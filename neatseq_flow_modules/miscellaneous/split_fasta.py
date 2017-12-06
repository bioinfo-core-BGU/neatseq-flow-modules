
""" 
A class that defines a module for splitting a BAM file by chromosome.

.. attention:: The module was tested on samtools 1.3


**Requires**:

* A SAM file in the following location:

    * ``sample_data[<sample>]["fastq"]["mapping"]["sam"]``

* A comma-separated list of chromosomes to split by and extract. **Make sure the list matches the chromosome names in the reference.** 

**Output**:

* Depending on the parameters, will put files in the following locations:

    * ``sample_data[<sample>]["fastq"]["mapping"]["bam"]``

.. csv-table:: Parameters that can be set:
    :header: "Parameter", "Values", "Comments"

    "chr_list", "*e.g.*: chr1,chr2", "A list of chromosomes to extract."
    
 
"""



import os
import sys
from PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "0.2.0"
class Step_split_by_chrom(Step):
       

    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances

        if not "chr_list" in self.params.keys():
            raise AssertionExcept("You must supply a 'chr_list' parameter.\n")
            
        self.params["chr_list"] = self.params["chr_list"].split(",")
        
            
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        
        # Checking a "mapping" exists for each sample:
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash

            # Check that a mapping slot exists
            try:
                self.sample_data[sample]["fastq"]["mapping"]
            except KeyError:
                raise AssertionExcept("No mapping slot defined.\n", sample)
            # Check that a sam or bam exists
            if "bam" in self.sample_data[sample]["fastq"]["mapping"]:
                self.file2use = "bam"
            elif "sam" in self.sample_data[sample]["fastq"]["mapping"]:
                self.file2use = "sam"
            else:
                raise AssertionExcept("Neither BAM nor SAM file exist for sample.\n", sample)

        # Create new sample names and dict mapping new names to old names:
        new_sample_dict = {}
        # Add new sample name lists
        for sample in self.sample_data["samples"]:
            # initiate with sample names
            new_sample_dict[sample] = []
            for chr in self.params["chr_list"]:
                new_sample_dict[sample].append("%s_%s" % (sample,chr))
                
        # print new_sample_dict
        
        temp_list = []
        for sample,smp_list in new_sample_dict.iteritems():
            temp_list.extend(smp_list)
        
        self.sample_data["original_samples"] = self.sample_data["samples"]
        self.sample_data["samples"] = temp_list
        self.sample_data["new_sample_index"] = new_sample_dict
        
        
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
        for sample in self.sample_data["original_samples"]:      # Getting list of samples out of samples_hash

            """samtools view -b in.bam chr1 > in_chr1.bam"""
        
            for chr in self.params["chr_list"]:
            
                new_sample = "%s_%s" % (sample,chr)
                # Make a dir for the current sample:
                sample_dir = self.make_folder_for_sample(new_sample)

                # Name of specific script:
                self.spec_script_name = "_".join([self.step,self.name,new_sample])
                self.script = ""
                

                # This line should be left before every new script. It sees to local issues.
                # Use the dir it returns as the base_dir for this step.
                use_dir = self.local_start(sample_dir)
                
                
                output_filename = "%s.chr_%s.bam" % (os.path.basename(self.sample_data[sample]["fastq"]["mapping"]["sam"]),chr)
                
                
                self.script += "%s view \\\n\t" % self.get_script_env_path()
                self.script += "-b \\\n\t"
                self.script += "-o %s \\\n\t" % (use_dir + output_filename)
                self.script += "%s \\\n\t" % (self.sample_data[sample]["fastq"]["mapping"][self.file2use])
                self.script += "%s \n\n" % (chr)
                
                
                self.sample_data[new_sample] = dict()
                self.sample_data[new_sample]["fastq"] = dict()
                self.sample_data[new_sample]["fastq"]["mapping"] = dict()
                
                self.sample_data[new_sample]["fastq"]["mapping"]["bam"] = sample_dir + output_filename

                self.stamp_file(self.sample_data[new_sample]["fastq"]["mapping"]["bam"])

                

                # Move all files from temporary local dir to permanent base_dir
                self.local_finish(use_dir,sample_dir)       # Sees to copying local files to final destination (and other stuff)
                
                
                self.create_low_level_script()
                        
