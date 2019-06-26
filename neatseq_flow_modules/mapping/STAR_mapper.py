# -*- coding: UTF-8 -*-
""" 
``STAR_mapper``
-----------------------------------------------------------------
:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A module for running STAR mapper:


Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


* fastq files in one of the following slots:

    * ``sample_data[<sample>]["fastq.F"]``
    * ``sample_data[<sample>]["fastq.R"]``
    * ``sample_data[<sample>]["fastq.S"]``
    
* If ``scope`` is set (must come after ``STAR_builder`` module which populates the required slots):
    
    * STAR index directories in:

        * ``sample_data[<sample>]["STAR.index"]``  if ``scope`` = "sample"
        * ``sample_data["STAR.index"]``  if ``scope`` = "project"

    * Reference fasta files in:

        * ``sample_data[<sample>]["STAR.fasta"]``  if ``scope`` = "sample"
        * ``sample_data["STAR.fasta"]``  if ``scope`` = "project"

    
Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


* Puts output sam files in the following slots:

    * ``self.sample_data[<sample>]["sam"]``

* Alternatively, if ``--outSAMtype`` is set to ``BAM``, puts output BAM files in the following slots:

    * ``self.sample_data[<sample>]["bam"]``
    * ``self.sample_data[<sample>]["bam_unsorted"]``

* High confidence collapsed splice junctions (SJ.out.tab  file) will be stored in:

    * ``self.sample_data[<sample>]["SJ.out.tab"]``

* If ``--quantMode`` contains ``TranscriptomeSAM``, alignments BAM translated into transcript coordinates will be stored in:

    * ``self.sample_data[<sample>]["TranscriptomeSAM"]``

* If ``--quantMode`` contains ``GeneCounts``, the ``ReadsPerGene.out.tab`` file will be stored:

    * ``self.sample_data[<sample>]["GeneCounts"]``

* If ``--outWigType`` is set, will store outputs in:

    * if ``--outWigType`` is ``wiggle``

        * ``self.sample_data[<sample>]["wig2_UniqueMultiple"]``
        * ``self.sample_data[<sample>]["wig2_Unique"]``
        * ``self.sample_data[<sample>]["wig1_UniqueMultiple"]``
        * ``self.sample_data[<sample>]["wig1_Unique"]``
        * ``self.sample_data[<sample>]["wig"]``

    * if ``--outWigType`` is ``bedGraph``

        * ``self.sample_data[<sample>]["bdg2_UniqueMultiple"]``
        * ``self.sample_data[<sample>]["bdg2_Unique"]``
        * ``self.sample_data[<sample>]["bdg1_UniqueMultiple"]``
        * ``self.sample_data[<sample>]["bdg1_Unique"]``
        * ``self.sample_data[<sample>]["bdg"]``
    
                
* Puts the name of the mapper in:
    ``self.sample_data[<sample>]["mapper"]``

* Puts fasta of reference genome (if one is given in param file) in:
    ``self.sample_data[<sample>]["reference"]``

Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "ref_genome", "path to genome fasta", ""
    "scope", "project | sample", "The scope from which to take the genome directory"

.. Note::
    You can set the RG atrribute of the resulting SAM/BAM files with the redirected parameter ``--outSAMattrRGline``
    This will set the equivalent STAR parameter.

    By default, the parameter will be set to include ID and SM tags, both set to the sample name.
    You can set the SM tag, but any ID tags will be removed and replaced with the sample name.

Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**For external index:**

::

    STAR_map:
        module:             STAR_mapper
        base:               STAR_bld_ind
        script_path:        /path/to/STAR
        redirects:
            --readMapNumber:    1000
            --genomeDir:        /path/to/genome/STAR_index/
            
            
**For project STAR index:**

::

    STAR_map:
        module:             STAR_mapper
        base:               STAR_bld_ind
        script_path:        /path/to/STAR
        scope:              project
        redirects:
            --readMapNumber:    1000
    
References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Dobin, A., Davis, C.A., Schlesinger, F., Drenkow, J., Zaleski, C., Jha, S., Batut, P., Chaisson, M. and Gingeras, T.R., 2013. **STAR: ultrafast universal RNA-seq aligner**. *Bioinformatics*, 29(1), pp.15-21.

"""


import os, re
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept

__author__ = "Menachem Sklarz"
__version__ = "1.6.0"

class Step_STAR_mapper(Step):

    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances

        self.auto_redirs = ["--readFilesCommand", "--readFilesIn", "--outFileNamePrefix", "--outTmpDir", "--outStd"]

        if "ref_genome" not in list(self.params.keys()):
            self.write_warning("No reference given with 'ref_genome' (path to fasta file). It is highly recommended to give one!\n")
        
        if "--runDirPerm" not in self.params["redir_params"]:
            self.params["redir_params"]["--runDirPerm"] = "All_RWX"
            self.write_warning("No --runDirPerm specified. Using 'All_RWX'")

        if "--outSAMtype" in self.params["redir_params"]:
            outSAMtype = re.split("\s+", self.params["redir_params"]["--outSAMtype"])
            if outSAMtype[0] not in ["SAM","BAM","None"]:
                raise AssertionExcept("Bad value for --outSAMtype: Has to be 'BAM', 'SAM' or 'None'")
            self.output_type = outSAMtype[0] 
            if self.output_type == "BAM":
                self.bam_types = outSAMtype[1:]
                if "Unsorted" not in self.bam_types and "SortedByCoordinate" not in self.bam_types:
                    raise AssertionExcept("When --outSAMtype is BAM, you must supply a type: "
                                          "'Unsorted', 'SortedByCoordinate' or both.")
        else:
            self.output_type = "SAM"
            
        if "--outSAMattrRGline" in self.params["redir_params"]:
            if re.match("ID\:\S+",self.params["redir_params"]["--outSAMattrRGline"]):
                self.write_warning("Removing 'ID:' from --outSAMattrRGline line!")
            self.params["outSAMattrRGline"] = re.sub("ID\:\S+","",self.params["redir_params"]["--outSAMattrRGline"])

        if "--outWigType" in self.params["redir_params"]:
            outWigType = re.split("\s+", self.params["redir_params"]["--outWigType"])
            if outWigType[0] not in ['None', 'bedGraph', 'wiggle']:
                raise AssertionExcept("Bad value for --outWigType: Has to be 'None', 'bedGraph' or 'wiggle'")
            self.wig_type = outWigType[0] 
            if self.wig_type == "wiggle":
                # See in build_scripts below. STAR produces 4 different wig files. Storing one in wig slot and others in wig* slots.
                self.write_warning("Saving UniqueMultiple wig from strand 1 as main 'WIG' file. If you want something else, you have to move it to the right slot...")
            elif self.wig_type == "bedGraph":
                self.write_warning("Saving UniqueMultiple bedGraph from strand 1 as main 'bdg' file. If you want something else, you have to move it to the right slot...")

        else:
            self.wig_type = "None"

    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """

        # Require either 'scope' or '--genomeDir':
        if "scope" in self.params:
            # If scope defined, comment if also -x exists.
            if "--genomeDir" in self.params["redir_params"]:
                raise AssertionExcept("Both 'scope' and '--genomeDir' specified!\n")

            # Loop over samples to set the reference genome:
            for sample in self.sample_data["samples"]:
                if self.params["scope"] == "project":
                    # Set project wide reference:
                    try:
                        self.sample_data[sample]["reference"] = self.sample_data["project_data"]["STAR.fasta"]
                    except:
                        raise AssertionExcept("No reference exists at 'project' scope. Do you have a STAR_builder step defined?")
                elif self.params["scope"] == "sample":
                    # Set per-sample reference:
                    try:
                        self.sample_data[sample]["reference"] = self.sample_data[sample]["STAR.fasta"]
                    except:
                        raise AssertionExcept("No reference exists at 'sample' scope. Do you have a STAR_builder step defined?",sample)
                else:
                    raise AssertionExcept("Scope must be either 'sample' or 'project'")

            if "ref_genome" in list(self.params.keys()):
                raise AssertionExcept("ref_genome was passed, and 'scope' was defined. Resolve!\n")
        else:
            # If scope is not defined, require '--genomeDir'
            if not "--genomeDir" in self.params["redir_params"]:
                raise AssertionExcept("Neither 'scope' nor '--genomeDir' specified.\n")
            # Storing reference genome for use by downstream steps:
            if "ref_genome" in list(self.params.keys()):
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
        
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash

            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)

            # Name of specific script:
            self.spec_script_name = self.set_spec_script_name(sample)
            self.script = ""

            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)

            # Define location and prefix for output files:
            output_prefix = sample + "_STAR_map"
            # output_prefix = use_dir + output_prefix
            
            # Adding sample ID to ID: attribute:
            if "outSAMattrRGline" in self.params:
                self.params["redir_params"]["--outSAMattrRGline"] = "ID:{ID} SM:{ID} {rest}".format(ID=sample, \
                                                                                rest=self.params["outSAMattrRGline"])
            else:
                self.params["redir_params"]["--outSAMattrRGline"] = "ID:{ID} SM:{ID}".format(ID=sample)
            
            # If using internal index, define it here:
            if "scope" in self.params:
                if self.params["scope"] == "sample":
                    self.params["redir_params"]["--genomeDir"] = self.sample_data[sample]["STAR.index"]
                else:   
                    self.params["redir_params"]["--genomeDir"] = self.sample_data["project_data"]["STAR.index"]
            self.script += """
if [ -e {tmpdir} ]; then 
    rm -rf {tmpdir}; 
fi\n\n""".format(tmpdir=use_dir+"STAR_tmp")

            # Get constant part of script:
            self.script += self.get_script_const()
            # Setting location of temporary dir:
            self.script += "--outTmpDir {tmp_dir} \\\n\t".format(tmp_dir=use_dir+"STAR_tmp")

            if "fastq.F" in self.sample_data[sample]:
                self.script += "--readFilesIn {fastqF} {fastqR}\\\n\t".\
                    format(fastqF=self.sample_data[sample]["fastq.F"],
                           fastqR=self.sample_data[sample]["fastq.R"])
            elif "fastq.S" in self.sample_data[sample]:
                self.script += "--readFilesIn {fastqS} \\\n\t".format(fastqS=self.sample_data[sample]["fastq.S"])
            else:
                raise AssertionExcept("No fastq files exist for sample!!\n" , sample)
        
            self.script += "--outFileNamePrefix %s%s. \n\n" % (use_dir,output_prefix)

            if self.output_type == "SAM":
                self.sample_data[sample]["sam"] = "%s%s.Aligned.out.sam" % (sample_dir,output_prefix)
                self.stamp_file(self.sample_data[sample]["sam"])
            elif self.output_type == "BAM":
                if "Unsorted" in self.bam_types:
                    self.sample_data[sample]["bam"] = "%s%s.Aligned.out.bam" % (sample_dir,output_prefix)
                    self.sample_data[sample]["bam_unsorted"] = "%s%s.Aligned.out.bam" % (sample_dir,output_prefix)
                    self.stamp_file(self.sample_data[sample]["bam_unsorted"])
                if "SortedByCoordinate" in self.bam_types:
                    self.sample_data[sample]["bam"] = "%s%s.Aligned.sortedByCoord.out.bam" % (sample_dir,output_prefix)
                self.stamp_file(self.sample_data[sample]["bam"])
            else:  # None
                pass

            # Storing the SJ file:
            self.sample_data[sample]["SJ.out.tab"] = "%s%s.SJ.out.tab" % (sample_dir,output_prefix)
            self.stamp_file(self.sample_data[sample]["SJ.out.tab"])

            if self.wig_type == "bedGraph":
                if "--outWigStrand" not in self.params["redir_params"] or self.params["redir_params"]["--outWigStrand"] == "Stranded":
                    self.sample_data[sample]["bdg2_UniqueMultiple"] = "%s%s.Signal.UniqueMultiple.str2.out.bg" % (sample_dir,output_prefix)
                    self.sample_data[sample]["bdg2_Unique"] = "%s%s.Signal.Unique.str2.out.bg" % (sample_dir,output_prefix)
                    self.stamp_file(self.sample_data[sample]["bdg2_UniqueMultiple"])
                    self.stamp_file(self.sample_data[sample]["bdg2_Unique"])
                self.sample_data[sample]["bdg1_UniqueMultiple"] = "%s%s.Signal.UniqueMultiple.str1.out.bg" % (sample_dir,output_prefix)
                self.sample_data[sample]["bdg1_Unique"] = "%s%s.Signal.Unique.str1.out.bg" % (sample_dir,output_prefix)
                self.stamp_file(self.sample_data[sample]["bdg1_UniqueMultiple"])
                self.stamp_file(self.sample_data[sample]["bdg1_Unique"])
                self.sample_data[sample]["bdg"] = self.sample_data[sample]["bdg1_UniqueMultiple"]
            elif self.wig_type == "wiggle":
                if "--outWigStrand" not in self.params["redir_params"] or self.params["redir_params"]["--outWigStrand"] == "Stranded":
                    self.sample_data[sample]["wig2_UniqueMultiple"] = "%s%s.Signal.UniqueMultiple.str2.out.wig" % (sample_dir,output_prefix)
                    self.sample_data[sample]["wig2_Unique"] = "%s%s.Signal.Unique.str2.out.wig" % (sample_dir,output_prefix)
                    self.stamp_file(self.sample_data[sample]["wig2_UniqueMultiple"])
                    self.stamp_file(self.sample_data[sample]["wig2_Unique"])
                self.sample_data[sample]["wig1_UniqueMultiple"] = "%s%s.Signal.UniqueMultiple.str1.out.wig" % (sample_dir,output_prefix)
                self.sample_data[sample]["wig1_Unique"] = "%s%s.Signal.Unique.str1.out.wig" % (sample_dir,output_prefix)
                self.stamp_file(self.sample_data[sample]["wig1_UniqueMultiple"])
                self.stamp_file(self.sample_data[sample]["wig1_Unique"])
                self.sample_data[sample]["wig"] = self.sample_data[sample]["wig1_UniqueMultiple"]
            else:
                pass

            if "--quantMode" in self.params["redir_params"]:
                if re.search(string=self.params["redir_params"]["--quantMode"], pattern="GeneCounts"):
                    self.sample_data[sample]["GeneCounts"] = "%s%s.ReadsPerGene.out.tab" % (sample_dir,output_prefix)
                    self.stamp_file(self.sample_data[sample]["GeneCounts"])
                if re.search(string=self.params["redir_params"]["--quantMode"], pattern="TranscriptomeSAM"):
                    self.sample_data[sample]["bam_transcriptome"] = "%s%s.Aligned.toTranscriptome.out.bam" % (sample_dir,output_prefix)
                    self.stamp_file(self.sample_data[sample]["bam_transcriptome"])

            # Storing name of mapper. might be useful:
            self.sample_data[sample]["mapper"] = self.get_step_step()  
            
            # Storing reference genome for use by downstream steps:
            if "ref_genome" in list(self.params.keys()):
                self.sample_data[sample]["reference"] = self.params["ref_genome"]

            # Move all files from temporary local dir to permanent base_dir
            self.local_finish(use_dir,self.base_dir)
            self.create_low_level_script()
