# -*- coding: UTF-8 -*-
""" 
``kaiju``
--------------------------------

:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

.. Note:: This module was developed as part of a study led by Dr. Jacob Moran Gilad


A module for running ``kaiju``:


Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* fastq files, either paired end or single:

    * ``sample_data[<sample>]["fastq.F"]``
    * ``sample_data[<sample>]["fastq.R"]``
    * ``sample_data[<sample>]["fastq.S"]``
    
    

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Puts the ``kaiju`` output files in:  

    * ``self.sample_data[<sample>]["raw_classification"]``
    
* If  'kaiju2krona' is set, puts the editted classification files in:

    * ``self.sample_data[<sample>]["classification"]``

* If ``ktImportText`` parameter was passed, puts the krona reports in

    * ``self.sample_data["project_data"]["krona"]``

    
Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "ktImportText", "", "A block with a ``path:`` element containing the path to ``ktImportText``, and a optionally a ``redirects`` element. If the path is left empty, will assume ``ktImportText`` is in the same dir as ``kaiju``"
    "kaiju2krona",  "", "A block with a ``path:`` element containing the path to ``kaiju2krona``, and a ``redirects`` element with ``-t`` and ``-n`` defined. If the path is left empty, will assume ``kaiju2krona`` is in the same dir as ``kaiju``"
    "use_fasta",  "", "Will use a ``fasta.nucl`` and NOT reads files"

.. Attention::
   Make sure to provide ``-t`` (nodes file) to kaiju via the (main) redirects block.

   Also, you must provide ``-n`` (names file) and ``-t`` (nodes file) to ``kaiju2krona`` via its redirects block. See example.

Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    kaiju1:
        module: kaiju
        base: trim1
        script_path: {Vars.paths.kaiju}
        redirects:
            -f: /path/to/kaijudb/kaiju_db.fmi
            -t: /path/to/kaijudb/nodes.dmp
            -z: 40            
        kaiju2krona:
            path:       '{Vars.Programs_path.kaiju.kaiju2krona}'
            redirects:
                -n:     '{Vars.databases.kaiju}/names.dmp'
                -t:     '{Vars.databases.kaiju}/nodes.dmp'
        ktImportText:
            path:       '{Vars.Programs_path.ktImportText}'

References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Menzel, P., Ng, K.L. and Krogh, A., 2016. **Fast and sensitive taxonomic classification for metagenomics with Kaiju**. Nature communications, 7.

"""




import os
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.6.0"



class Step_kaiju(Step):
    """ A class that defines a pipeline step name (=instance).
    """
    
    def step_specific_init(self):
        """ Called on intiation
            Good place for parameter testing.
            Wrong place for sample data testing
        """
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        self.file_tag = ".kaiju.out"
        if "scope" not in list(self.params.keys()):
            self.params["scope"]="sample"
            
        # Checking this once and then applying to each sample:
        try:
            self.params["redir_params"]
        except KeyError:
            raise AssertionExcept("You must specify -t and -f as redirected params.\n")
            
        if "-t" not in list(self.params["redir_params"].keys()) or "-f" not in list(self.params["redir_params"].keys()):
            raise AssertionExcept("You must specify -t and -f as redirected params.\n")


        if "kaiju2krona" in list(self.params.keys()):
            if self.params["kaiju2krona"] == None:
                raise AssertionExcept("Please provide -t and -n via the 'redirects' to 'kaiju2krona'")
            elif isinstance(self.params["kaiju2krona"], str):
                raise AssertionExcept("You must supply -t and -n via the 'redirects' to 'kaiju2krona'")
            else:
                if "-n" not in self.params["kaiju2krona"]["redirects"] or "-t" not in self.params["kaiju2krona"]["redirects"]:
                    raise AssertionExcept("Please provide -t and -n via the 'redirects' to 'kaiju2krona'")
                if self.params["kaiju2krona"]["path"] is None:
                    self.params["kaiju2krona"]["path"] = os.sep.join([os.path.basename(self.params["script_path"]),
                                                                        "kaiju2krona"])

        if "ktImportText" in list(self.params.keys()):
            if self.params["ktImportText"] == None:

                self.params["ktImportText"] = {"path": os.sep.join([os.path.basename(self.params["script_path"]),
                                                                       "ktImportText"])}
            elif isinstance(self.params["ktImportText"], str):
                self.params["ktImportText"] = {"path":self.params["ktImportText"]}
            else:
                if self.params["ktImportText"]["path"] is None:
                    self.params["ktImportText"]["path"] = os.sep.join([os.path.basename(self.params["script_path"]),
                                                                        "ktImportText"])




    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """

        if self.params["scope"] == "project":
            sample_list = ["project_data"]
        elif self.params["scope"] == "sample":
            sample_list = self.sample_data["samples"]
        else:
            raise AssertionExcept("'scope' must be either 'sample' or 'project'")
        #sample_list = self.sample_data["samples"]
        for sample in sample_list:  # Getting list of samples out of samples_hash
            if 'use_fasta' in list(self.params.keys()):
                if "fasta.nucl" not in list(self.sample_data[sample].keys()):
                    raise AssertionExcept("No Nucleotide FASTA in: \n", sample)
            else:
                if len({"fastq.F", "fastq.R"} & set(self.sample_data[sample].keys())) == 1:
                    raise AssertionExcept(
                        "Sample has only forward or reverse reads. It must have either pairs or single reads\n", sample)
                # if len({"fastq.F", "fastq.R", "fastq.S"} & set(self.sample_data[sample].keys())) ==3:
                    # raise AssertionExcept("Kaiju is not defined for mixed paired and single reads\n", sample)



    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        if self.params["scope"] == "project":
            sample_list = ["project_data"]
        elif self.params["scope"] == "sample":
            sample_list = self.sample_data["samples"]
            
        if "ktImportText" in self.params:
            if "redirects" in self.params["ktImportText"]:
                redirects = " \\\n\t".join(
                    [key + " " + (val if val else "")
                     for key, val
                     in self.params["ktImportText"]["redirects"].items()])
            else:
                redirects = ""

            krona_report_fn = self.base_dir + self.sample_data["Title"] + "_krona_report.html"
            self.script = "# Creating krona html reports\n"
            # Adding env and setenv lines to script
            self.script += self.get_setenv_part()
            # Main part of script:
            self.script += "%s \\\n\t" % self.params["ktImportText"]["path"]
            if redirects:
                self.script += "%s \\\n\t" % redirects
            self.script += "-o %s \\\n\t" % krona_report_fn
            for sample in sample_list:      # Getting list of samples out of samples_hash
                self.script += "%s,%s \\\n\t" % (self.sample_data[sample]["classification"],sample)
            # Removing extra \\
            self.script = self.script.rstrip("\\\n\t") 

            # Storing and stamping results:
            self.sample_data["project_data"]["krona"] = krona_report_fn
            self.stamp_file(self.sample_data["project_data"]["krona"])

        else:
            self.write_warning("You did not supply 'ktImportText'. Will not create krona reports...\n")
            self.script = ""

    def build_scripts(self):
        """ This is the actual script building function
            Most, if not all, editing should be done here 
            HOWEVER, DON'T FORGET TO CHANGE THE CLASS NAME AND THE FILENAME!
        """
        if self.params["scope"] == "project":
            sample_list = ["project_data"]
        elif self.params["scope"] == "sample":
            sample_list = self.sample_data["samples"]
            
        for sample in sample_list:      # Getting list of samples out of samples_hash
            
            # Name of specific script:
            self.spec_script_name = self.set_spec_script_name(sample)
            self.script = ""

            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)
                
            # Define output filename
            output_filename = "".join([use_dir , sample , self.file_tag])

            self.script += self.get_script_const()
            if 'use_fasta' in list(self.params.keys()):
                  self.script += "-i %s \\\n\t" % self.sample_data[sample]["fasta.nucl"]
            else:
                # Adding reads
                if ("fastq.F" in self.sample_data[sample]) and ("fastq.R" in self.sample_data[sample]):
                    self.script += "-i %s \\\n\t" % self.sample_data[sample]["fastq.F"]
                    self.script += "-j %s \\\n\t" % self.sample_data[sample]["fastq.R"]
                elif "fastq.S" in self.sample_data[sample]:
                    self.script += "-i %s \\\n\t" % self.sample_data[sample]["fastq.S"]
                else:
                    raise AssertionExcept("Weird combination of reads!\n")

            self.script += "-o %s\n\n" % output_filename

            # Storing the output file in $samples_hash
            self.sample_data[sample]["raw_classification"]        = "%s" % (sample_dir + os.path.basename(output_filename))
            self.stamp_file(self.sample_data[sample]["raw_classification"])

            if "kaiju2krona" in self.params:

                if "redirects" in self.params["kaiju2krona"]:
                    if "-u" not in self.params["kaiju2krona"]["redirects"].keys():
                        self.params["kaiju2krona"]["redirects"]["-u"] = ""
                    redirects = " \\\n\t".join(
                        [key + " " + (val if val else "")
                         for key, val
                         in self.params["kaiju2krona"]["redirects"].items()])
                    redirects = "\n\t{redirs}\\".format(redirs=redirects)
                else:
                    redirects = " \\\n\t".join("-u")
                self.script += """

# Creating text report for krona
{path} \\{redirects}
\t-i {input} \\
\t-o {input}.4krona.txt                 
                
                """.format(path=self.params["kaiju2krona"]["path"],
                           redirects=redirects,
                           input=self.sample_data[sample]["raw_classification"])

                self.sample_data[sample]["classification"] = "%s.4krona.txt" % self.sample_data[sample]["raw_classification"]
                self.stamp_file(self.sample_data[sample]["classification"])

            # Move all files from temporary local dir to permanent base_dir
            self.local_finish(use_dir,self.base_dir)
            self.create_low_level_script()

    def make_sample_file_index(self):
        """ Make file containing samples and target file names for use by kraken analysis R script
        """
        pass
