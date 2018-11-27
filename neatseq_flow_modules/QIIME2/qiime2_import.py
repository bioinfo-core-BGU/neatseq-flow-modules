# -*- coding: UTF-8 -*-
""" 
``qiime2_import`` 
-----------------------------------------------------------------
:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A module for running ``TransDecoder`` on a transcripts file.

.. Note:: Tested on TransDecoder version 5.5.0.. The main difference being that in this version an output directory can be specified in the command line.

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
``fasta`` files in at least one of the following slots:
    
    * ``sample_data[<sample>]["fasta.nucl"]``  (if ``scope = sample``)
    * ``sample_data["fasta.nucl"]``  (if ``scope = project``)

    
Output:
~~~~~~~~~~~~~

* If ``scope = project``:

    * Protein fasta in ``self.sample_data["project_data"]["fasta.prot"]``
    * Gene fasta in ``self.sample_data["project_data"]["fasta.nucl"]``
    * Original transcripts in ``self.sample_data["project_data"]["transcripts.fasta.nucl"]``
    * GFF file in ``self.sample_data["project_data"]["gff3"]``

* If ``scope = sample``:

    * Protein fasta in ``self.sample_data[<sample>]["fasta.prot"]``
    * Gene fasta in ``self.sample_data[<sample>]["fasta.nucl"]``
    * Original transcripts in ``self.sample_data[<sample>]["transcripts.fasta.nucl"]``
    * GFF file in ``self.sample_data[<sample>]["gff3"]``


                
Parameters that can be set        
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"

    "scope", "sample|project", "Determine weather to use sample or project transcripts file."
    
    
Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    trino_Transdecode_highExpr:
        module:             TransDecoder
        base:               Split_Fasta
        script_path:        {Vars.paths.TransDecoder}
        scope:              sample
        
        
References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"""



import os
import sys
import re
from neatseq_flow.PLC_step import Step,AssertionExcept
from qiime2_functions import *
import yaml
try:
    # from yaml import CLoader as Loader
    from yaml import CSafeLoader as Loader
except ImportError:
    from yaml import SafeLoader as Loader


__author__ = "Menachem Sklarz"
__version__ = "1.1.0"


class Step_qiime2_import(Step):
    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances


        # Read YAML of plugin arguments
        with open(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                               "importable_types_and_formats.yaml"),"r") as fileh:
            filelines = fileh.readlines()

        self.qiime_types_formats = yaml.load("".join(filelines),  Loader=yaml.Loader)

        if "scope" not in self.params:
            self.params["scope"] = "project"

        
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
            Here you should do testing for dependency output. These will NOT exist at initiation of this instance. They are set only following sample_data updating
        """

        if "--type" in self.params["redir_params"]:
            self.type = self.params["redir_params"]["--type"]

            if "--input-format" not in self.params["redir_params"]:
                raise AssertionExcept("Please supply an input-format to use!",sample)

            if self.type == "SampleData[SequencesWithQuality]":
                # Check single reads exist
                for sample in self.sample_data["samples"]:
                    # if "fastq.S" not in self.sample_data[sample]:
                    if not any(map(lambda x: x in self.sample_data[sample], ["fastq.F", "fastq.S"])):
                        raise AssertionExcept(
                            "For type {type}, you must have single or paired end files in samples".format(
                                type=self.type))
            elif self.type == "SampleData[PairedEndSequencesWithQuality]":
                # Check single reads exist
                for sample in self.sample_data["samples"]:
                    if not all(map(lambda x: x in self.sample_data[sample], ["fastq.F", "fastq.R"])):
                        raise AssertionExcept(
                            "For type {type}, you must have paired end files in samples".format(type=self.type))
            else:
                raise AssertionExcept("For type '{type}', please supply an '--input-path'.".format(type=self.type))

        else:
            self.type = None

            # if len(self.qiime_types) == 0:   # No QIIME types to import. Using sequences
            # if "type" in self.params:
            #     # raise AssertionExcept("If not importing QIIME2 types, you have to have read files and define a "
            #     #                       "'type' parameter")
            #     imp_types  = {"SampleData[SequencesWithQuality]":["fastq.F", "fastq.S"],
            #                  "SampleData[PairedEndSequencesWithQuality]":["fastq.F", "fastq.R"]}
            #     for imp_type in imp_types:
            #         if self.params["type"] == imp_type:
            #             if sample=="project_data":
            #                 for sample in self.sample_data["samples"]:
            #                     if not any(map(lambda x: x in self.sample_data[sample], imp_types[imp_type])):
            #                         raise AssertionExcept("For type {type}, you must have the following types in "
            #                                               "samples: {types}".format(type=self.params["type"],
            #                                                                         types=", ".join(imp_types[imp_type])))
            #             else:
            #                 if not any(map(lambda x: x in self.sample_data[sample], imp_types[imp_type])):
            #                     raise AssertionExcept("For type {type}, you must have the following types in "
            #                                           "sample: {types}".format(type=self.params["type"],
            #                                                                     types=", ".join(imp_types[imp_type])),
            #                                           sample)
            #
            #     if "--input-path" not in self.params["redir_params"]:
            #
            #         if self.params["type"] == "SampleData[SequencesWithQuality]":
            #             # Check single reads exist
            #             for sample in self.sample_data["samples"]:
            #                 # if "fastq.S" not in self.sample_data[sample]:
            #                 if not any(map(lambda x: x in self.sample_data[sample], ["fastq.F", "fastq.S"])):
            #                     raise AssertionExcept("For type {type}, you must have single or paired end files in samples".format(type=self.type))
            #         elif self.type == "SampleData[PairedEndSequencesWithQuality]":
            #             # Check single reads exist
            #             for sample in self.sample_data["samples"]:
            #                 if not all(map(lambda x: x in self.sample_data[sample], ["fastq.F","fastq.R"])):
            #                     raise AssertionExcept("For type {type}, you must have paired end files in samples".format(type=self.type))
            #         else:
            #             raise AssertionExcept("For type '{type}', please supply an '--input-path'.".format(type=self.type))



    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        
        pass

    def build_scripts(self):


        # if self.type in ["SampleData[PairedEndSequencesWithQuality]", "SampleData[SequencesWithQuality]"] and \
        #     self.sample_data["project_data"]:

        if self.type:

            # Name of specific script:
            self.spec_script_name = self.set_spec_script_name()
            self.script = ""

            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample="project_data")

            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)

            # Create manifest file
            manifest_fn = use_dir + self.sample_data["Title"] + ".manifest"
            manifest_file = open(manifest_fn, "w")
            manifest_file.write("# Created by NeatSeq-Flow\n")
            manifest_file.write("sample-id,absolute-filepath,direction\n")
            direction_dict = {"fastq.F":"forward","fastq.R":"reverse","fastq.S":"forward"}

            if self.type == "SampleData[PairedEndSequencesWithQuality]":
                directions_to_include = {"fastq.F", "fastq.R"}
            else:
                directions_to_include = {"fastq.F", "fastq.S"}
            for sample in self.sample_data["samples"]:
                for direction in list(directions_to_include & set(self.sample_data[sample])):
                    manifest_file.write("{sampleid},{absolutefilepath},{direction}\n".
                                        format(sampleid=sample,
                                               absolutefilepath=self.sample_data[sample][direction],
                                               direction=direction_dict[direction]))

            manifest_file.close()
            # Create script
            self.script = self.get_script_const()

            self.script += """\
\t--output-path {dir}{title}.import.qza \\
\t--input-path {manifest} 
""".format(dir=use_dir,
           title=self.sample_data["Title"],
           manifest=manifest_fn)

            self.sample_data["project_data"][self.type] = "{dir}{title}.import.qza".\
                format(dir=use_dir,
                       title=self.sample_data["Title"])
            self.stamp_file(self.sample_data["project_data"][self.type])
            self.local_finish(use_dir, sample_dir)

            self.create_low_level_script()

        self.qiime_types = list(set(self.qiime_types_formats["importable-types"]) & set(self.sample_data["project_data"].keys()))
        print self.qiime_types
        # sys.exit()

        for qtype in self.qiime_types:
            print self.sample_data["project_data"][qtype]
            # Name of specific script:
            self.spec_script_name = self.jid_name_sep.join([self.step,self.name,edit_qiime_types(qtype)])

            self.script = ""

            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample="project_data")

            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)

            # Create script
            self.script = self.get_script_const()

            self.script += """\
--type {qtype} \\
\t--input-path {inp_path} \\{inp_format}
\t--output-path {dir}{title}.{type}.import.qza  
""".format(dir=use_dir,
           qtype=qtype,
           type=edit_qiime_types(qtype),
           inp_path=self.sample_data["project_data"][qtype][0],
           inp_format="" if len(self.sample_data["project_data"][qtype])==1
               else "\n\t--input-format " + self.sample_data["project_data"][qtype][1] + " \\",
           title=self.sample_data["Title"])

            self.sample_data["project_data"][qtype] = "{dir}{title}.{type}.import.qza". \
                format(dir=use_dir,
                       type=edit_qiime_types(qtype),
                       title=self.sample_data["Title"])
            self.stamp_file(self.sample_data["project_data"][qtype])
            self.local_finish(use_dir, sample_dir)

            self.create_low_level_script()
