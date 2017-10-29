


""" A module for running bowtie2 mapper:

The reads stored in each sample are aligned to one of the following bowtie2 indices:

* An external index passed with the ``-x`` parameter.
* A bowtie2 index on a project fasta files, such as an assembly from all samples. Specify with ``bowtie2_mapper:scope  project``
* A sample bowtie2 index on a sample-specific fasta file, such as from a sample-wise assembly or from the sample file. Specify with ``bowtie2_mapper:scope  sample``

The latter two options **must** come after a ``bowtie2_builder`` instance.

.. tip:: See the documentation for the ``bowtie2_builder`` module.

.. note:: fastq files are never defined project-wide
    
    The ``scope`` parameter controls the origin of the index files, *i.e.* wheather the fasta file to map to is an assembly of the sample reads (scope: sample) or an assembly of all reads in the project (scope: project). The reads to be mapped are always saple reads, as a 'fastq' slot is not defined at the project level.
    
Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* fastq files in one of the following slots:

    * ``sample_data[<sample>]["fastqc"]["fastq.F"]``
    * ``sample_data[<sample>]["fastqc"]["fastq.R"]``
    * ``sample_data[<sample>]["fastqc"]["fastq.S"]``
    

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Puts output sam files in the following slots:
    * ``self.sample_data[<sample>]["fastq"]["mapping"]["sam"]``

* Puts the name of the mapper in:
    * ``self.sample_data[<sample>]["fastq"]["mapping"]["mapper"]``

* puts fasta of reference genome (if one is given in param file) in:
    * ``self.sample_data[<sample>]["fastq"]["mapping"]["reference"]``

Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "-x", "path to bowtie2 index", "If not given, will look for a project bowtie2 index and then for a sample bowtie2 index"
    "ref_genome", "path to genome fasta", "If -x is NOT given, will use the equivalent internal fasta. If -x is passed, and ref_genome is NOT passed, will leave the reference slot empty"
    "get_map_log", "", "Store the log produced by bowtie2 (This is bowtie2 standard output)"
    "scope", "project | sample", "Indicates whether to use a project or sample bowtie2 index."

Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**For external index:**

::

    bwt2_1:
        module: bowtie2_mapper
        base: trim1
        script_path: /path/to/bowtie2
        qsub_params:
            -pe: shared 20
        get_map_log:
        ref_genome: /path/to/ref_genome.fna
        redirects:
            -p: 20
            -q: null
            -x: /path/to/bowtie2_index/ref_genome

**Using a bowtie2 index constructed from a project fasta:**

::

    bwt2_1:
        module: bowtie2_mapper
        base: bwt2_bld1
        script_path: /path/to/bowtie2
        qsub_params:
            -pe: shared 20
        get_map_log:
        scope: project
        redirects:
            -p: 20
            -q: null


"""


import os
import sys
from PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.1.0"


class Step_bowtie2_mapper(Step):
    
   
    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        # self.file_tag = "Bowtie_mapper"

    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        
        
        # # Initializing a "mapping" dict for each sample:
        # for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash

            # try:
                # self.sample_data[sample]["fastq"]["mapping"]
            # except KeyError:
                # self.sample_data[sample]["fastq"]["mapping"] = {}
            # else:
                # self.write_warning("mapping dict exists for sample %s. Double mapping steps?\n", sample)

        
        # Require either 'scope' or '-x':
        if "scope" in self.params:
            # If scope defined, comment if also -x exists.
            if "-x" in self.params["redir_params"]:
                raise AssertionExcept("Both 'scope' and '-x' specified!\n")

            try:
                # Loop over samples to set the reference genome:
                for sample in self.sample_data["samples"]:
                    if self.params["scope"] == "project":
                        # Set project wide reference:
                        self.sample_data[sample]["reference"] = self.sample_data["bowtie2_fasta"]
                    elif self.params["scope"] == "sample":
                        # Set per-sample reference:
                        self.sample_data[sample]["reference"] = self.sample_data[sample]["bowtie2_fasta"]
                    else:
                        raise AssertionExcept("Scope must be either 'sample' or 'project'")
                
            except KeyError:
                raise AssertionExcept("There is a mismatch between 'scope' and the existing bowtie2 index\n", sample)
                
            if "ref_genome" in self.params.keys():
                raise AssertionExcept("ref_genome was passed, and 'scope' was defined. Ignoring ref_genome\n")
        else:
            # If scope is not defined, require '-x'
            if not "-x" in self.params["redir_params"]:
                raise AssertionExcept("Neither 'scope' nor '-x' specified.\n")
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

            # Name of specific script:
            self.spec_script_name = "_".join([self.step,self.name,sample])
            self.script = ""
            
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)
 
            # Define location and prefix for output files:
            output_prefix = sample + "_bowtie2_map"
            output_prefix = use_dir + output_prefix
            
            # Get constant part of script:
            self.script += self.get_script_const()
            
            self.script += "--rg-id %s \\\n\t" % sample
            self.script += "--rg   SM:%s \\\n\t" % sample
            
            try:  # If scope was passed, include either project or sample bowtie2 index
                if self.params["scope"] == "project":
                    self.script += "-x %s \\\n\t" % self.sample_data["bowtie2_index"]
                else:
                    self.script += "-x %s \\\n\t" % self.sample_data[sample]["bowtie2_index"]
            except KeyError:  # Otherwise do nothing - '-x' is included through redirect params
                pass
            
                
            # assert set("fastq.F","fastq.R","fastq.S") & self.sample_data["sample"].keys(), "There are no reads for sample %s" % sample
            if "fastq.F" in self.sample_data[sample].keys():
                self.script += "-1 %s \\\n\t-2 %s \\\n\t" % (self.sample_data[sample]["fastq.F"],self.sample_data[sample]["fastq.R"])
            if "fastq.S" in self.sample_data[sample].keys():
                self.script += "-U %s \\\n\t" % self.sample_data[sample]["fastq.S"]

            
            
            
            self.script += "--met-file %s.stats \\\n\t" % output_prefix
            self.script += "-S %s.sam " % output_prefix
            if "get_map_log" in self.params.keys():
                self.script += "\\\n\t2> %s.log\n\n" % output_prefix
            else:
                self.script += " \n\n";


            self.sample_data[sample]["sam"] = "%s.sam" % output_prefix
            self.stamp_file(self.sample_data[sample]["sam"])
            
            # Storing name of mapper. might be useful:
            self.sample_data[sample]["mapper"] = self.get_step_step()  
            

   
            # Move all files from temporary local dir to permanent base_dir
            self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
       
            
            
            self.create_low_level_script()
                    
