# -*- coding: UTF-8 -*-
""" 
``macs2_bdgcmp``
---------------------------


:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A module for running macs2 bdgcmp:

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Files in the following slots:

    * ``self.sample_data[<sample>]["control_lambda"]`` - Control BedGraph
    * ``self.sample_data[<sample>]["treat_pileup"]`` - Treatment BedGraph


Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Puts output macs2 output files in the following slots:

    * ``self.sample_data[<sample>]["bdg"])`` - The comparison bedgraph!
    * ``self.sample_data[<sample>]["bigwig"])`` - if ``slop_path`` and ``ucscTools_path`` were passed
    * ``self.sample_data[<sample>]["wig"])`` - if ``slop_path`` and ``ucscTools_path`` were passed
    * ``self.sample_data[<sample>]["tdf"])`` - in TDF format (if ``slop_path``, ``ucscTools_path`` and ``toTDF_path`` were passed)

    
Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "slop_path", "path to bedtools slop", "Is part of the process for converting bdg files into bigwig and wig"
    "ucscTools_path", "path to ucscTools", "UCSCtools bedClip, bedGraphToBigWig and bigWigToWig are part of the process for converting bdg files into bigwig and wig"
    "toTDF_path", "path to toTDF", "Converts the wig file into TDF file."
    "genome", "path to chrom.sizes for reference genome", "If running bedToBigBed, you must supply the genome chrom.sizes file."
    




Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    bdgcmp:
        module: macs2_bdgcmp
        base: macs1
        script_path: /path/to/macs2 bdgcmp
        genome: /path/to/chrom.sizes.txt
        slop_path: /path/to/bin/bedtools slop
        ucscTools_path: /path/to/ucscTools/bin
        toTDF_path: /path/to/bin/java -Xmx1500m -jar /path/to/igvtools.jar toTDF
        redirects:
            --method: FE

References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Feng, J., Liu, T., Qin, B., Zhang, Y. and Liu, X.S., 2012. **Identifying ChIP-seq enrichment using MACS**. *Nature protocols*, 7(9), pp.1728-1740.

"""

import os
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.6.0"


class Step_macs2_bdgcmp(Step):
    """ A class that defines a pipeline step name (=instance).
    """
    


    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        self.file_tag = "macs2"

        if self.params["ucscTools_path"]:
            self.params["ucscTools_path"] = self.params["ucscTools_path"].rstrip(os.sep) + os.sep

    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """

        # Initializing a "mapping" dict for each sample:
        for sample in self.sample_data["Controls"]:      # Getting list of samples out of samples_hash

            # Make sure there are ChIP results
            try:
                self.sample_data[sample]["control_lambda"]
                self.sample_data[sample]["treat_pileup"]
            except KeyError:
                raise AssertionExcept("chip_seq results do not contain lambda and pileup files for sample %s.\n",sample)
        
        
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
        for sample in list(self.sample_data["Controls"].keys()):      # Getting list of samples out of Controls dict.

            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)

            # Name of specific script:
            self.spec_script_name = self.set_spec_script_name(sample)
            self.script = ""

            # Name of control sample:
            control = self.sample_data["Controls"][sample]

            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)

            if "-m" in self.params["redir_params"]:
                method2use = self.params["redir_params"]["-m"]
            elif "--method" in self.params["redir_params"]:
                method2use = self.params["redir_params"]["--method"]
            else:
                method2use = "ppois"

            # Defined full path to output filename
            output_filename = "%s.%s_%s.bdg" % (sample, self.file_tag, method2use)
            
                
                
            self.script += self.get_script_const()
                
            # Add lines for sample mapping files:
            self.script += "-t %s \\\n\t" % self.sample_data[sample]["treat_pileup"]
            if not "nocontrol" in list(self.params.keys()):
                self.script += "-c %s \\\n\t" % self.sample_data[sample]["control_lambda"]
            else:
                print("Running bdgcmp with no control defined. I don;t know how this will work...\n")

                
            # Add output directory
            self.script += "-o %s \n\n" % "".join([use_dir, output_filename])
            
            if "slop_path" in self.params and "ucscTools_path" in self.params:
                self.script += """\n\n
# Running slop | bedClip
{slop_path} \\
        -i {INPUT_BDG} \\
        -g {genome} \\
        -b 0 | \\
{ucscTools_path}bedClip \\
    stdin \\
    {genome} \\
    {OUT_CLIP}

# Running bedGraphToBigWig
{ucscTools_path}bedGraphToBigWig \\
    {OUT_CLIP} \\
    {genome} \\
    {OUT_BW}

# Running bigWigToWig 
{ucscTools_path}bigWigToWig  \\
    {OUT_BW} \\
    {OUT_W}


rm -f {OUT_CLIP}
""".format(slop_path= self.params["slop_path"] ,
           ucscTools_path= self.params["ucscTools_path"] ,
           genome= self.params["genome"] ,
           INPUT_BDG= "%s%s.%s_%s.bdg" % (use_dir,sample, self.file_tag, method2use),
           OUT_BW= "%s%s.%s_%s.bw" % (use_dir,sample, self.file_tag, method2use),
           OUT_W= "%s%s.%s_%s.wig" % (use_dir,sample, self.file_tag, method2use),
           OUT_CLIP= "%s%s.%s_%s.clip" % (use_dir,sample, self.file_tag, method2use))

                if "toTDF_path" in list(self.params.keys()):
                    self.script += """\n\n
# Converting to TDF:
{toTDF_path}  \\
    {OUT_W} \\
    {OUT_TDF} \\
    {genome} 
    
""".format(toTDF_path= self.params["toTDF_path"] ,
           genome= self.params["genome"] ,
           OUT_TDF= "%s%s.%s_%s.tdf" % (use_dir,sample, self.file_tag, method2use),
           OUT_W= "%s%s.%s_%s.wig" % (use_dir,sample, self.file_tag, method2use))

                    self.sample_data[sample]["tdf"] = "".join([sample_dir, "%s.%s_%s.tdf" % (sample, self.file_tag, method2use)])

                    # Stamping files produced in this step:
                    self.stamp_file(self.sample_data[sample]["tdf"])
                # End if "toTDF_path" in self.params.keys
        
                # Storing the output file in $samples_hash
                self.sample_data[sample]["bigwig"] = "".join([sample_dir, "%s.%s_%s.bw" % (sample, self.file_tag, method2use)])
                self.sample_data[sample]["wig"] = "".join([sample_dir, "%s.%s_%s.wig" % (sample, self.file_tag, method2use)])

                # Stamping files produced in this step:
                self.stamp_file(self.sample_data[sample]["bigwig"])
                self.stamp_file(self.sample_data[sample]["wig"])
            # End if "slop_path" in self.params and "ucscTools_path" in self.params
                
            # Storing the resulting bedgraph in mapping->bdg. 
            # This will occupy the bdg slot from macs2 and from mapping! (=will overwrite the slot, not the files)
            self.sample_data[sample]["bdg"] = "".join([sample_dir, output_filename])
            self.stamp_file(self.sample_data[sample]["bdg"])
            
            # Wrapping up function. Leave these lines at the end of every iteration:
            self.local_finish(use_dir,sample_dir)       # Sees to copying local files to final destination (and other stuff)

            
            self.create_low_level_script()
                    
