# -*- coding: UTF-8 -*-
""" 
``qiime2_import`` 
-----------------------------------------------------------------
:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A module for running ``qiime tools import`` on various importable types

.. Note:: Tested on ``qiime2`` version ``2018.11``

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
* If importing reads:

    * ``sample_data[<sample>]["fastq.F|R|S"]``

* If importing other types, requires that type in the sample file

  The file can be defined in the sample file either just as a path, or as a path, format pair, as follows:

  Only path::

    EMPSingleEndSequences	/path/to/emp-single-end-sequences

  Path, format pair::

    EMPSingleEndSequences	/path/to/emp-single-end-sequences
    EMPSingleEndSequences	EMPPairedEndDirFmt


Output:
~~~~~~~~~~~~~

* If importing reads, will create the imported artifact in one of:

    * ``sample_data["project_data"]["SampleData[SequencesWithQuality]"]``
    * ``sample_data["project_data"]["SampleData[PairedEndSequencesWithQuality]"]``

* If importing other types:

    * ``sample_data["project_data"]["<type imported>"]``



                
Parameters that can be set        
~~~~~~~~~~~~~~~~~~~~~~~~~~


    
Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Importing paired end reads::

    import_reads:
        module:                     qiime2_import
        base:                       trim1
        script_path:                qiime tools import
        redirects:
            --type:                 SampleData[PairedEndSequencesWithQuality]
            --input-format:         PairedEndFastqManifestPhred33

Importing internal types::

    merge_data:
        module:         merge
        src:            EMPSingleEndSequences
        trg:            EMPSingleEndSequences
        script_path:    ..import..
        scope:          project

    import:
        module:         qiime2_import
        base:           merge_data
        script_path:    qiime tools import

        
References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Bolyen, E., Rideout, J.R., Dillon, M.R., Bokulich, N.A., Abnet, C., Al-Ghalith, G.A., Alexander, H., Alm, E.J., Arumugam, M. and Asnicar, F., 2018. **QIIME 2: Reproducible, interactive, scalable, and extensible microbiome data science**. *PeerJ Preprints* 6:e27295v1 https://doi.org/10.7287/peerj.preprints.27295v1


"""



import os
import os.path
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
                               "importable_types_and_formats.yml"),"r") as fileh:
            filelines = fileh.readlines()

        self.qiime_types_formats = yaml.load("".join(filelines),  Loader=yaml.Loader)

        if "scope" not in self.params:
            self.params["scope"] = "project"

        
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
            Here you should do testing for dependency output. These will NOT exist at initiation of this instance. They are set only following sample_data updating
        """

        # Types that do not require formats.
        # Not used. At the moment, not checking at all for existance of --input-format
        self.types_no_formats = ["EMPSingleEndSequences"]

        if "--type" in self.params["redir_params"]:
            self.type = self.params["redir_params"]["--type"]


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
            elif self.type in self.sample_data["project_data"]:
                if isinstance(self.sample_data["project_data"][self.type],list) and len(self.sample_data["project_data"][self.type])>1:
                    self.params["redir_params"]["--input-format"] = self.sample_data["project_data"][self.type][1]
                # elif "--input-format" not in self.params["redir_params"]:
                #     raise AssertionExcept("Please define a --input-format for type {type}".format(type=self.type))

            else:
                raise AssertionExcept("For type '{type}', please supply an '--input-path'.".format(type=self.type))

            # print self.params["redir_params"]["--input-format"]
            # print self.type
            # sys.exit()
            # if "--input-format" not in self.params["redir_params"]:
            #     raise AssertionExcept("Please supply an --input-format to use!")

            if "--input-format" in self.params["redir_params"]:
                if self.params["redir_params"]["--input-format"] not in self.qiime_types_formats["importable-formats"]:
                    raise AssertionExcept("Unrecognised format {format} for type {type}".
                                          format(format=self.params["redir_params"]["--input-format"],
                                                 type=self.type))

        else:
            self.type = None

        # Create list of types to import. In two steps:
        # 1. Get list of types in project_data, which are legitimate qiime2 types:
        self.params["qiime_types"] = list(
            set(self.qiime_types_formats["importable-types"]) & set(self.sample_data["project_data"].keys()))
        # 2. Keep only types which are either
        #   - strings (so they were imported by a merge step) or
        #   - a list of length>1 where the 2nd element is a legitimate qiime2 format or
        #   - a list of length=1
        self.params["qiime_types"] = [type for type
                                      in self.params["qiime_types"]
                                      if isinstance(self.sample_data["project_data"][type],str) or
                                      (isinstance(self.sample_data["project_data"][type],list) and
                                       (len(self.sample_data["project_data"][type]) > 1 and
                                        self.sample_data["project_data"][type][1] in
                                        self.qiime_types_formats["importable-formats"]) or
                                       (len(self.sample_data["project_data"][type]) == 1))]
        # 3. Remove those with qza or qzv extensions. They are already imported. This will be done in the script
        # building below
        # 4. If type requested in in the list, keep only it:
        if self.type in self.params["qiime_types"]:
            self.params["qiime_types"] = [self.type]
            self.type = None

        # print self.get_step_name()
        # print self.params["qiime_types"]
        # print self.type
        # print self.params["redir_params"]
        # sys.exit()

    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        
        pass

    def build_scripts(self):


        # If type is defined, import only the requested type

        if self.type:


            # Name of specific script:
            self.spec_script_name = self.set_spec_script_name()
            self.script = ""

            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample="project_data")

            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)

            if self.type in ["SampleData[PairedEndSequencesWithQuality]","SampleData[SequencesWithQuality]"]:

                # Create manifest file
                manifest_fn = use_dir + self.sample_data["Title"] + ".manifest"
                manifest_file = open(manifest_fn, "w")
                manifest_file.write("# Created by NeatSeq-Flow\n")
                manifest_file.write("sample-id,absolute-filepath,direction\n")
                direction_dict = {"fastq.F": "forward", "fastq.R": "reverse", "fastq.S": "forward"}

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
                input_path = manifest_fn
                output_path = "{dir}{title}.import.qza".format(dir=use_dir,title=self.sample_data["Title"])

            else:
                print self.params
                sys.exit()
            # Create script
            self.script = self.get_script_const()

            self.script += """\
--output-path {outp} \\
\t--input-path {inp} 
""".format(inp=input_path,
           outp=output_path)

            self.sample_data["project_data"][self.type] = "{dir}{title}.import.qza". \
                format(dir=sample_dir,
                       title=self.sample_data["Title"])
            self.stamp_file(self.sample_data["project_data"][self.type])
            self.local_finish(use_dir, sample_dir)

            self.create_low_level_script()

        # If type is not defined, import all available types:
        else:
            for qtype in self.params["qiime_types"]:

                # Skipping qza's and qzv's:
                if isinstance(self.sample_data["project_data"][qtype], str) and \
                        os.path.splitext(self.sample_data["project_data"][qtype])[1] in [".qza",".qzv"]:
                    self.write_warning("Skipping type {qtype}. Extension is a qiime extension. Already imported!".format(qtype=qtype))
                    # print qtype
                    # print self.sample_data["project_data"][qtype]
                    # print os.path.splitext(self.sample_data["project_data"][qtype])
                    continue

                self.spec_script_name = self.jid_name_sep.join([self.step, self.name, edit_qiime_types(qtype)])

                self.script = ""

                # Make a dir for the current sample:
                sample_dir = self.make_folder_for_sample(sample="project_data")

                # This line should be left before every new script. It sees to local issues.
                # Use the dir it returns as the base_dir for this step.
                use_dir = self.local_start(sample_dir)
                input_format = ""
                if isinstance(self.sample_data["project_data"][qtype], list):
                    input_path = self.sample_data["project_data"][qtype][0]
                    if len(self.sample_data["project_data"][qtype]) == 1:
                        try:
                            input_format = self.params["redir_params"]["--input-format"]
                        except KeyError:
                            # Here you can try fitting an input-format to the type, or test it there is a default
                            pass
                            # raise AssertionExcept("--input-format not defined anywhere")
                    else:
                        input_format = self.sample_data["project_data"][qtype][1]
                else:
                    input_path = self.sample_data["project_data"][qtype]
                    try:
                        input_format = self.params["redir_params"]["--input-format"]
                    except KeyError:
                        pass
                        # raise AssertionExcept("--input-format not defined anywhere")

                if input_format:
                    input_format = "\n\t--input-format {inform} \\".format(inform=input_format)

                # Create script
                self.script = """\
{script_path} \\
\t--type {qtype} \\
\t--input-path {inp_path} \\{inp_format}
\t--output-path {dir}{title}.{type}.import.qza  
""".format(script_path=self.params["script_path"],
           dir=use_dir,
           qtype=qtype,
           type=edit_qiime_types(qtype),
           inp_path=input_path,
           inp_format=input_format,
           title=self.sample_data["Title"])

                self.sample_data["project_data"][qtype] = "{dir}{title}.{type}.import.qza". \
                    format(dir=sample_dir,
                           type=edit_qiime_types(qtype),
                           title=self.sample_data["Title"])
                self.stamp_file(self.sample_data["project_data"][qtype])
                self.local_finish(use_dir, sample_dir)

                self.create_low_level_script()

