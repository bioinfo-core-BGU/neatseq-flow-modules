# -*- coding: UTF-8 -*-
""" 
``makePICARDfiles`` 
-----------------------------------------------------------------
:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A class that prepares the various files required by PICARD and GATK:

* ``gtf`` -> ``bed`` conversion
* ``dict`` file from ``fasta`` file
* ``bed`` -> ``interval_list`` conversion
* ``bed`` -> ribosomal ``interval_list``
* ``gtf`` -> ``refFlat`` conversion


Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
* A nucleotide fasta file in 
    * self.sample_data["fasta.nucl"]
    * self.sample_data[sample]["fasta.nucl"]
* A ``gtf`` file in:
    * self.sample_data["gtf"]
    * self.sample_data[sample]["gtf"]

.. Attention:: If ``scope`` is set to ``sample``, all of the above files should be in the sample scope!

Also, requires the following programs:

* ``gtf2bed`` (from ``BEDOPS``)
* ``gtfToGenePred`` (kent utils. Can be installed from bioconda)


Output:
~~~~~~~~~~~~~

* puts Trinotate report file in:

    * ``sample_data["bed"]`` 
    * ``sample_data["dict"]`` 
    * ``sample_data["interval_list"]`` 
    * ``sample_data["rRNA.interval_list"]`` 
    * ``sample_data["refFlat"]`` 
        
If ``scope`` is set to ``sample``, all the above will be in the sample scope.
                
Parameters that can be set        
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"

    "scope", "sample|project", ""
    "gtfToGenePred_path","",""
    "gtf2bed_path","",""
    
Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    makePICARDfiles:
        module:             MakePICARDfiles
        base:               
        script_path:        {Vars.paths.PICARD}
        scope:              project
        gtfToGenePred_path: /path/to/gtfToGenePred
        gtf2bed_path:       /path/to/gtf2bed
        setenv:             PATH="/path/to/bedops/bin:$PATH"

  
        
References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


"""



import os
import sys
import re
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.1.0"

class Step_MakePICARDfiles(Step):
    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        
        self.arg_separator='='

    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
            Here you should do testing for dependency output. These will NOT exist at initiation of this instance. They are set only following sample_data updating
        """
        
        if "scope" not in self.params:
            raise AssertionExcept("No 'scope' specified.")
        elif self.params["scope"]=="project":
            if "fasta.nucl" not in self.sample_data or "gtf" not in self.sample_data:
                raise AssertionExcept("Project does not have fasta.nucl and gtf files.")

        elif self.params["scope"]=="sample":
            
            for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
                if "fasta.nucl" not in self.sample_data[sample] or "gtf" not in self.sample_data[sample]:
                    raise AssertionExcept("Sample does not have fasta.nucl and gtf files.", sample)
        else:
            raise AssertionExcept("'scope' must be either 'sample' or 'project'")

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

        bed_file = "{title}.bed".format(title = self.sample_data["Title"])
        dict_file = "{title}.dict".format(title = self.sample_data["Title"])
        interval_list_file = "{title}.interval_list".format(title = self.sample_data["Title"])
        rRNA_interval_list_file = "{title}.rRNA_interval_list".format(title = self.sample_data["Title"])
        refFlat_file = "{title}.refFlat".format(title = self.sample_data["Title"])

        # Get constant part of script:
        self.script += self.get_setenv_part()
        self.script += """

{gtf2bed} \\
    < {gtf} \\
    > {bed}

{PICARD} \\
    CreateSequenceDictionary \\
      R={fasta} \\
      O={dict}
      
{PICARD} \\
    BedToIntervalList \\
      I={bed} \\
      O={interval_list} \\
      SD={dict}
      
grep -i rrna {bed} |
    {PICARD} \\
        BedToIntervalList \\
          I=/dev/stdin \\
          O={rRNA_interval_list} \\
          SD={dict}

# Make refFlat.txt:
{gtfToGenePred} \\
    -genePredExt \\
    -geneNameAsName2 \\
    -ignoreGroupsWithoutExons \\
    {gtf} \\
    /dev/stdout | \\
    awk 'BEGIN {{ OFS="\\t"}} {{print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}}' \\
    > {refFlat}
""".format(PICARD=self.params["script_path"],
           gtf2bed=self.params["gtf2bed_path"],
           gtfToGenePred=self.params["gtfToGenePred_path"],
           gtf=self.sample_data["gtf"],
           fasta=self.sample_data["fasta.nucl"],
           bed=use_dir + bed_file,
           dict=use_dir + dict_file,
           interval_list=use_dir + interval_list_file,
           rRNA_interval_list=use_dir + rRNA_interval_list_file,
           refFlat=use_dir + refFlat_file)
    
        # Store results to fasta and assembly slots:
        self.sample_data["bed"]                = self.base_dir +  bed_file
        self.sample_data["sequence.dict"]               = self.base_dir +  dict_file
        self.sample_data["INTERVAL_LIST"]      = self.base_dir +  interval_list_file
        self.sample_data["RIBOSOMAL_INTERVALS"] = self.base_dir +  rRNA_interval_list_file
        self.sample_data["REF_FLAT"]            = self.base_dir +  refFlat_file
 
        self.stamp_file(self.sample_data["bed"])
        self.stamp_file(self.sample_data["sequence.dict"])
        self.stamp_file(self.sample_data["INTERVAL_LIST"])
        self.stamp_file(self.sample_data["RIBOSOMAL_INTERVALS"])
        self.stamp_file(self.sample_data["REF_FLAT"])
    
        # Move all files from temporary local dir to permanent base_dir
        self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
        self.create_low_level_script()
                    
#################################################
    def build_scripts_sample(self):
        
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash

        # Name of specific script:
            self.spec_script_name = self.set_spec_script_name(sample)
            self.script = ""


            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)
            
            bed_file = "{title}.bed".format(title = sample)
            dict_file = "{title}.dict".format(title = sample)
            interval_list_file = "{title}.interval_list".format(title = sample)
            rRNA_interval_list_file = "{title}.rRNA_interval_list".format(title = sample)
            refFlat_file = "{title}.refFlat".format(title = sample)
     
            
            # Get constant part of script:
            self.script += self.get_setenv_part()
            self.script += """

{gtf2bed} \\
    < {gtf} \\
    > {bed}

{PICARD} \\
    CreateSequenceDictionary \\
      R={fasta} \\
      O={dict}
      
{PICARD} \\
    BedToIntervalList \\
      I={bed} \\
      O={interval_list} \\
      SD={dict}
      
grep -i rrna {bed} |
    {PICARD} \\
        BedToIntervalList \\
          I=/dev/stdin \\
          O={rRNA_interval_list} \\
          SD={dict}

# Make refFlat.txt:
{gtfToGenePred} \\
    -genePredExt \\
    -geneNameAsName2 \\
    -ignoreGroupsWithoutExons \\
    {gtf} \\
    /dev/stdout | \\
    awk 'BEGIN {{ OFS="\\t"}} {{print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}}' \\
    > {refFlat}
""".format(PICARD=self.params["script_path"],
           gtf2bed=self.params["gtf2bed_path"],
           gtfToGenePred=self.params["gtfToGenePred_path"],
           gtf=self.sample_data["gtf"],
           fasta=self.sample_data["fasta.nucl"],
           bed=use_dir + bed_file,
           dict=use_dir + dict_file,
           interval_list=use_dir + interval_list_file,
           rRNA_interval_list=use_dir + rRNA_interval_list_file,
           refFlat=use_dir + refFlat_file)
    
            # Store results to fasta and assembly slots:
            self.sample_data[sample]["bed"]                = self.base_dir + bed_file
            self.sample_data[sample]["sequence.dict"]               = self.base_dir + dict_file
            self.sample_data[sample]["INTERVAL_LIST"]      = self.base_dir + interval_list_file
            self.sample_data[sample]["RIBOSOMAL_INTERVALS"] = self.base_dir + rRNA_interval_list_file
            self.sample_data[sample]["REF_FLAT"]            = self.base_dir + refFlat_file
     
            self.stamp_file(self.sample_data[sample]["bed"])
            self.stamp_file(self.sample_data[sample]["sequence.dict"])
            self.stamp_file(self.sample_data[sample]["INTERVAL_LIST"])
            self.stamp_file(self.sample_data[sample]["RIBOSOMAL_INTERVALS"])
            self.stamp_file(self.sample_data[sample]["REF_FLAT"])
     
            # Move all files from temporary local dir to permanent base_dir
            self.local_finish(use_dir,sample_dir)       # Sees to copying local files to final destination (and other stuff)
            self.create_low_level_script()
