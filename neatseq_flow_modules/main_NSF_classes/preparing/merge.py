# -*- coding: UTF-8 -*-
"""
``merge`` :sup:`*`
-------------------

:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A module for merging and staging files from the sample file into **NeatSeq-Flow**. 

Can be used in two modes: 

**The Basic mode**
    **NeatSeq-Flow** will attempt to guess all the parameters it requires. Files will be concatenated and stored in the file type index according to the table below. File types not included in the table will be stored in the file type index by the type specified in the sample file, in lower case.
    
    You have to make sure that all files of each file type have the same extension for **NeatSeq-Flow** to guess the ``script_path`` and ``pipe`` parameters. 
    
    
**The Advanced mode**
    is used when more control on data merging is required. It enables full control over which file types are merged, how they are copied and in which slots they are placed in the file type index. It also enables merging file types not recognized by **NeatSeq-Flow** (see list below).
    
    In this mode, you have to define the following lists: ``src``, ``trg``, ``script_path``, ``scope`` and ``ext``. For each file type in the sample file, you should have an entry in the ``src`` list. The other lists should apply to the equivalent entry in ``src``. ``trg`` is the target file type (in the file type index) for the merged files, ``script_path`` is the shell command to use to merge the source type, ``scope`` is the scope for which the source type is defined and ``ext`` is the suffix to append to the merged filenames. Strings are expanded to the length of ``src`` list, so if ``script_path`` is the same for all source types, it is enough to specify it once. 
    
    When using the Advanced mode, by passing the ``src`` list, you must also define the other lists, *i.e.* ``trg``, ``ext``, ``scope`` and ``script_path``. However, **NeatSeq-Flow** will try guessing the lists based on the lists of recognized file types and extensions. 
    
    If some of the file types in ``src`` are recognized and some are not, you can pass the lists mentioned above with values for the unrecognized types, leaving *null* in the positions of the recognized types. These *null* values will be guessed by **NeatSeq-Flow**.
    
    The advanced mode is experimental, and documentation will hopefully improve as we gain experience with it.

.. Note:: **Definition of ``script_path`` in the ``merge`` module**: ``script_path`` should be a shell program that receives a list of files and produces one single output file **to the standard error**. Examples of such programs are ``cat`` for text files and ``gzip -cd`` for gzipped files. Other types of compressed files should have such a command as well. 


.. Tip:: **NeatSeq-Flow** attempts to guess the ``script_path`` and ``pipe`` values based on the input file extensions. For this to work, leave the ``script_path`` and ``pipe`` lists empty and make sure all files from the same source have the same extensions (*e.g.* all gzipped files should have *.gz* as file extension). 

    If you want **NeatSeq-Flow** to guess only some of the ``script_path`` values, set them to `null` or to ``..guess..``, *e.g.* if ``src`` is ``[Single,TYP1]`` and ``script_path`` is ``[null,cat]``, then the ``script_path`` for *Single* will be guessed and the ``script_path`` for *TYP1* will be set to *cat*.

    Two more options are available for ``script path``: ``..skip..`` will skip the type entirely, while ``..import..`` will import the values from the sample file into the relevant slots without actually producing any scripts (This is useful for including entities which are not files in the sample file. `e.g.` in the qiime2 pipeline you might want to include a semantic type in the sample file).


    The following extensions are recognized:

    .. csv-table:: File extensions recognized by **NeatSeq-Flow**
        :header: "Extension", "``script_path``", "``pipe``"
        :widths: 20,20,60
        
        ".fasta", "cat",""
        ".faa", "cat",""
        ".fna", "cat",""
        ".txt", "cat",""
        ".csv", "cat",""
        ".fastq", "cat",""
        ".fa", "cat",""
        ".fq", "cat",""
        ".gz", "gzip -cd", ""
        ".zip", "echo", 'xargs -d " " -I % sh -c "unzip -p %"'
        ".bz2","bzip -cd",""
        ".dsrc2","echo", 'xargs -d " " -I % sh -c "dsrc2 d -s %"'
        ".dsrc","echo", 'xargs -d " " -I % sh -c "dsrc d -s %"'
        

    
.. _sample definition file: http://neatseq-flow.readthedocs.io/en/latest/02.build_WF.html#sample-file-definition

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* For the basic mode:
    * A list of files of the following types, either in ``[<sample>]`` or in ``[project_data]``:

.. csv-table:: File types recognized by **NeatSeq-Flow**
    :header: "Source", "Target"

    "Forward",      "fastq.F"
    "Reverse",      "fastq.R"
    "Single",       "fastq.S"
    "Nucleotide",   "fasta.nucl"
    "Protein",      "fasta.prot"
    "SAM",          "sam"
    "BAM",          "bam"
    "REFERENCE",    "reference"
    "VCF",          "vcf"
    "G.VCF",        "g.vcf"

* For the Advanced mode:
    * Lists of files in any file type, either in ``[<sample>]`` or in ``[project_data]``.
    
Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Merged files of the types in the table above are placed in slots according to the types in the 2nd column of the table. 

.. attention:: If you want to do something more complex with the combined files, you can use the ``pipe`` parameter to send extra commands to be piped on the files after the main command. **This is an experimental feature and should be used with care**. 

    e.g.: You can get files from a remote location by setting ``script_path`` to ``curl`` and ``pipe`` to ``gzip -cd``. This will download the files with curl, unzip them and concatenate them into the target file.  In the sample file, specify remote URLs instead of local pathes. **This will work only for one file per sample**.

    As of version 1.3.0, ``pipe`` can be a list of the same length as ``src`` and it we be treated like the other lists describe above.

Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 10, 20, 70

    "script_path", "", "The shell command to use for merging the source files."
    "src", "", "A list of source file types as the appear in the sample file."
    "trg", "", "A list of target file type for the merged files."
    "scope", "sample | project", "The scope at which each of the sources can be found."
    "ext", "", "The suffix to append to the merged filename."
    "pipe", "", "Additional commands to be piped on the files before writing to file."

Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Basic mode, gzipped files::

    merge1:
        module: merge
        script_path: gzip -cd

Basic mode, remote files::

    merge1:
        module: merge
        script_path: curl
        pipe:  gzip -cd

Advanced mode, mixture of types and scopes::

    merge1:
        module:         merge
        src:            [UR1,       UR2]
        script_path:    [gzip -cd,  cat]
        scope:          [sample,    project]
        trg:            [unrecog1,  unrecog2]
        ext:            [ur1,       ur2]

        
Advanced mode, both recognized and unrecognized file types::

    merge1:
        module:         merge
        src:            [UR1,       Forward,    Reverse]
        script_path:    [gzip -cd,  null,       null]
        scope:          # Guess!
        trg:            [unrecog1,  null,       null]
        ext:            [ur1,       null,       null]
        
Advanced mode, same types in samples and project::

    merge1:
        module:         merge
        src:            [Nucleotide,    Nucleotide]
        script_path:    [cat,           cat]
        scope:          [sample,        project]
        trg:            
        ext:            
        

"""

import os
import sys
import re
from neatseq_flow.PLC_step import Step,AssertionExcept

import yaml

from pprint import pprint as pp

__author__ = "Menachem Sklarz"
__version__ = "1.1.0"


class Step_merge(Step):
    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        self.file_tag = "merge"
        

        # Load YAML of file type stored in merge_file_types.yml
        with open(os.path.join(os.path.dirname(os.path.realpath(__file__)),"merge_file_types.yml"),"r") as fileh:
            try:
                self.default_src_trg_map = yaml.load("".join(fileh.readlines()),  Loader=yaml.SafeLoader)
            except yaml.YAMLError, exc:
                if hasattr(exc, 'problem_mark'):
                    mark = exc.problem_mark
                    print "Error position: (%s:%s)" % (mark.line+1, mark.column+1)
                    print mark.get_snippet()
                raise AssertionExcept("Error loading file types index 'merge_file_types.yml'")
            except:
                raise AssertionExcept("Error loading file types index 'merge_file_types.yml'")
            
        # Load YAML of script_paths stored in merge_script_path_types.yml
        with open(os.path.join(os.path.dirname(os.path.realpath(__file__)),"merge_script_path_types.yml"),"r") as fileh:
            try:
                self.script_path_map = yaml.load("".join(fileh.readlines()),  Loader=yaml.SafeLoader)
            except yaml.YAMLError, exc:
                if hasattr(exc, 'problem_mark'):
                    mark = exc.problem_mark
                    print "Error position: (%s:%s)" % (mark.line+1, mark.column+1)
                    print mark.get_snippet()
                raise AssertionExcept("Error loading script_path index 'merge_script_path_types.yml'")
            except:
                raise AssertionExcept("Error loading script_path index 'merge_script_path_types.yml'")
            
    
    def step_sample_initiation(self):
        """ Two operations are performed in this function:
            1. Creating good lists for both Basic and Advanced modes.
            2. Testing the values passed by the user and guessing values where necessary.
        """

        # Getting list of possible src values.
        # This is done for sample and project scope separately.
        # Used for creating src or testing src list passed by user
        src = []
        scope = []
        # Get list of existing file types in samples file:
        for sample in self.sample_data["samples"]:
            sample_src = list(set(self.sample_data[sample].keys()))
            if "type" in sample_src: 
                sample_src.remove("type")   # 'type' is the type of sample (PE, etc.)
            sample_scope = ["sample"] * len(sample_src)
            # If sample data exists, store in 'src' and 'scope'
            # Do this only if one of the following:
                # 1. 'scope' is not defined, or 
                # 2. 'scope' is a string equal to 'sample' or
                # 3. 'scope' is a list containing 'sample'.
            if "scope" not in self.params \
                or (isinstance(self.params["scope"],str) and self.params["scope"] == "sample") \
                or (isinstance(self.params["scope"],list) and "sample" in self.params["scope"]):
                src += sample_src
                scope += sample_scope
        if "project_data" in self.sample_data:
            project_src = list(set(self.sample_data["project_data"].keys()))
            project_scope = ["project"] * len(project_src)
            # If project data exists, add to 'src' and 'scope'
            # Do this only if one of the following:
                # 1. 'scope' is not defined, or 
                # 2. 'scope' is a string equal to 'project' or
                # 3. 'scope' is a list containing 'project'.
            if "scope" not in self.params \
                or (isinstance(self.params["scope"],str) and self.params["scope"] == "project") \
                or (isinstance(self.params["scope"],list) and "project" in self.params["scope"]):
                src = src + project_src
                scope = scope + project_scope

        # Getting unique pairs of src and scope:
        uniq_src_scope = list(set(zip(src,scope)))
        src, scope = zip(*uniq_src_scope)
        src = list(src)
        scope = list(scope)
        # src = sample_src + project_src
        # scope = sample_scope + project_scope
        # If 'src' is NOT user-defined: (Basic mode)
        if "src" not in self.params or not self.params["src"]:
            self.params["src"] = src
            # For each of the other required lists
            for param in ["trg","script_path","ext","pipe","scope"]:
                if param in self.params and self.params[param]:  # list is defined!
                    if param in ['trg','ext']:
                        raise AssertionExcept("'src' not specified. Please do not specify '%s'." %param)
                    else:  # Param is script_path or pipe or scope
                        if isinstance(self.params[param], str):
                            # If param is a str, convert it into a list:
                            self.params[param] = re.split(pattern="\s*,\s*", string=self.params[param])
                            
                        if isinstance(self.params[param], list):
                            if len(self.params[param])>1:
                                raise AssertionExcept("'src' not specified. Please do not specify '%s' as list "
                                                      "(only string values are accepted)." %param)
                            else:
                                self.params[param] = self.params[param] * len(self.params["src"])
                        else:
                            raise AssertionExcept("Unrecognized value for '%s'" % param)
                           
                else: # 4. If not passed by user, creating list with None. This is so that pipe can be populated by automatic file extension recognition
                    self.params[param] = [None] * len(self.params["src"])
            # Defining 'scope' (basic mode. user should not pass this. Setting it here)
            self.params["scope"] = scope

            
        else: # 'src' is user-defined (Advanced mode)
            # Converting 'src' into list
            if isinstance(self.params['src'], str):
                self.params['src'] = re.split(pattern="\s*,\s*", string=self.params['src'])
            if isinstance(self.params['src'], list):
                pass
            else:  # src is not a string or a list. What is it???
                raise AssertionExcept("Unrecognized format in 'src'")
            for param in ["trg","script_path","ext","pipe","scope"]:
                if param in self.params and self.params[param]:  # list is defined!

                    if isinstance(self.params[param], str):
                        self.params[param] = re.split(pattern="\s*,\s*", string=self.params[param])
                    if isinstance(self.params[param], list):
                        if len(self.params[param]) == 1:
                            self.params[param] = self.params[param] * len(self.params["src"])
                        elif len(self.params[param]) == len(self.params["src"]):
                            pass
                        else:
                            raise AssertionExcept("Parameter '%s' must be a single value or a list the length of "
                                                  "'src'. Set null in places to be determined automatically.")
                    else: # param is not a string or a list. What is it???
                        raise AssertionExcept("Unrecognized format in '%s'" % param)
                else:
                    self.params[param] = [None] * len(self.params["src"])

        #
        # #---------------------------------------
        # for param in ["script_path","src","trg","ext","pipe","scope"]:
        #     print param
        #     pp(self.params[param])
        # #---------------------------------------
        # sys.exit()
        #
                    
        # For each src in the list of sources:
        for src_ind in range(len(self.params["src"])):  
            src = self.params["src"][src_ind]
            # print "---> ",src
            scope = self.params["scope"][src_ind]
            script_path = self.params["script_path"][src_ind]
            trg = self.params["trg"][src_ind]
            ext = self.params["ext"][src_ind]
            # A list of srcs to remove. These are sources that do not exist in samples or projects
            bad_srcs = []
            # Guessing 'trg'
            if not trg:
                if src not in self.default_src_trg_map.keys():
                    self.write_warning("The following 'src' is  not recognized: {src}. "
                                       "Setting 'trg' to {trg}".format(src=src,trg=src.lower()))
                    self.params["trg"][src_ind] = src.lower()
                else:
                    self.params["trg"][src_ind] = self.default_src_trg_map[src][0] 
                # Guessing 'ext'
            if ext == None:
                if src not in self.default_src_trg_map.keys():
                    self.write_warning("The following 'src' is  not recognized: {src}. "
                                       "Setting 'ext' to {ext}".format(src=src,ext=src.lower()))
                    self.params["ext"][src_ind] = src.lower()
                else:
                    self.params["ext"][src_ind] = self.default_src_trg_map[src][1] 

            # Guessing scope if None
            if not scope:
                if all([src in self.sample_data[x].keys() for x in self.sample_data["samples"]]):
                    self.params["scope"][src_ind] = "sample"
                elif src in self.sample_data["project_data"]:
                    self.params["scope"][src_ind] = "project"
                else:
                    raise AssertionExcept("{src} does not exist in all samples or in project. Make sure sample "
                                          "file is correct, or pass 'scope' explicitly.".format(src=src))
            scope = self.params["scope"][src_ind]
            # Testing scope and guessing 'script_path':
            if scope=="sample":
                for sample in self.sample_data["samples"]:
                    if src not in self.sample_data[sample]:
                        self.write_warning("Type '{src}' does not exist for sample '{smp}'!".format(src=src,smp=sample))
                        bad_srcs += [src_ind]  # Adding bad source to bad_srcs
                # print "==>", bad_srcs
                # Guessing script_path:
                # Get file extensions:
                if not script_path or script_path == "..guess..":
                    # Is none or ..guess.. - try guessing
                    # src_exts is defined as follows: For each sample in samples list,
                    # get the list of file extensions. Creates a list of lists.

                    src_exts = [[os.path.splitext(filename)[1]
                                 for filename
                                 in self.sample_data[sample][src]]
                                for sample
                                in self.sample_data["samples"]
                                if src in self.sample_data[sample]]

                    # Flatten the list of lists, and uniqify:
                    src_exts = list(set([item for sublist in src_exts for item in sublist]))

                    if len(src_exts)>1:
                        pass
                    else:
                        # Convert set to string:
                        src_exts = src_exts[0]
                        if src_exts not in self.script_path_map.keys():
                            raise AssertionExcept("Unidentified extension in source '{src}' ({ext}). Can't guess "
                                                  "'script_path'".format(src=src, ext=src_exts))
                        else:
                            if isinstance(self.script_path_map[src_exts],list):
                                self.params["script_path"][src_ind] = self.script_path_map[src_exts][0]
                                self.params["pipe"][src_ind] = self.script_path_map[src_exts][1]
                            else:
                                self.params["script_path"][src_ind] = self.script_path_map[src_exts]
                        # print "===> ",src_exts

            elif scope=="project":
                if src not in self.sample_data["project_data"]:
                    self.write_warning("Type '{src}' does not exist in project data!".format(src=src))
                    bad_srcs += [src_ind]  # Adding bad source to bad_srcs

                # Guessing script_path:
                # Get file extensions:
                if not script_path or script_path == "..guess..":
                    # Is none or ..guess.. - try guessing
                    src_exts = list(set([os.path.splitext(filename)[1] for filename in self.sample_data["project_data"][src]]))
                    if len(src_exts)>1:
                        raise AssertionExcept("More than one file extension in source '{src}' for project ({ext}). "
                                              "Can't guess 'script_path'".format(src=src, ext=", ".join(src_exts)))
                    # Convert set to string:
                    src_exts = src_exts[0]
                    if src_exts not in self.script_path_map.keys():
                        raise AssertionExcept("Unidentified extension in source '{src}' for project ({ext}). "
                                              "Can't guess 'script_path'".format(src=src, ext=src_exts))
                    else:
                        if isinstance(self.script_path_map[src_exts],list):
                            self.params["script_path"][src_ind] = self.script_path_map[src_exts][0]
                            self.params["pipe"][src_ind] = self.script_path_map[src_exts][1]
                        else:
                            self.params["script_path"][src_ind] = self.script_path_map[src_exts]
                                    
            else:
                raise AssertionExcept("'scope' must be either 'sample' or 'project'")

        # #---------------------------------------
        # for param in ["script_path","src","trg","ext","pipe","scope"]:
            # print param
            # # self.params[param] = [i for j, i in enumerate(self.params[param]) if j not in bad_srcs]
            # self.params[param] = [i for j, i in enumerate(self.params[param]) ]
            # pp(self.params[param])
        # print bad_srcs
        # # ---------------------------------------
        # sys.exit()

    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        
        pass
    
    def build_scripts(self):

        for scope_ind in range(len(self.params["scope"])):
            src = self.params["src"][scope_ind]
            scope = self.params["scope"][scope_ind]
            trg = self.params["trg"][scope_ind]
            ext = self.params["ext"][scope_ind]
            script_path = self.params["script_path"][scope_ind]
            pipe = self.params["pipe"][scope_ind]

            if scope == "sample":
                # Each iteration must define the following class variables:
                    # spec_script_name
                    # script
                for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
                    # General comment: If there is a parallel routine for each direction (forward, reverse), add this loop	
                    # if  in self.sample_data[sample].keys():

                    # The following two may be modified per sample. Therefore, reading them again for each sample
                    script_path = self.params["script_path"][scope_ind]
                    pipe = self.params["pipe"][scope_ind]

                    # src_type not defined for this sample. Move on.
                    if src not in self.sample_data[sample]:
                        continue

                    if script_path == "..import..":
                        self.sample_data[sample][trg] = self.sample_data[sample][src]
                        return
                    if script_path == "..skip..":
                        return

                    self.spec_script_name = self.jid_name_sep.join([self.step,self.name,sample,src])
                    
                    # This line should be left before every new script. It sees to local issues.
                    # Use the dir it returns as the base_dir for this step.
                    use_dir = self.local_start(self.base_dir)
                    
                    fq_fn = ".".join([sample, src, self.file_tag,ext])          #The filename containing the end result. Used both in script and to set reads in $sample_params

                    if not script_path or script_path == "..guess..":
                        # Not all samples have the same file types. Sample-specific guessing...

                        src_exts = list(set([os.path.splitext(filename)[1]
                                             for filename
                                             in self.sample_data[sample][src]]))
                        if len(src_exts)>1:
                            raise AssertionExcept("More than one file extension in source '{src}' for sample "
                                                  "'{sample}' ({ext}). Can't guess 'script_path'".
                                                  format(src=src, sample = sample, ext=", ".join(src_exts)))
                        # Convert set to string:
                        src_exts = src_exts[0]
                        if src_exts not in self.script_path_map.keys():
                            raise AssertionExcept("Unidentified extension in source '{src}' for sample {sample} "
                                                  "({ext}). Can't guess 'script_path'".
                                                  format(src=src, sample = sample, ext=src_exts))
                        else:
                            if isinstance(self.script_path_map[src_exts],list):
                                script_path = self.script_path_map[src_exts][0]
                                pipe = self.script_path_map[src_exts][1]
                            else:
                                script_path = self.script_path_map[src_exts]
                    # Composing script:
                    self.script = ""
                    self.script += script_path + " \\\n\t"
                    # The following line concatenates all the files in the direction separated by a " "
                    self.script += " ".join(self.sample_data[sample][src]) 
                    self.script += " \\\n\t"
                    if pipe:  # pipe is not 'None'
                        self.script += "| {pipe} \\\n\t".format(pipe = pipe)
                    self.script += "> %s%s \n\n"  % (use_dir, fq_fn)

                    # Move all files from temporary local dir to permanent base_dir
                    self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)

                    # Store file in active file for sample:
                    self.sample_data[sample][trg] = self.base_dir + fq_fn
                    
                    self.stamp_file(self.sample_data[sample][trg])

                    self.create_low_level_script()
                    
            elif scope == "project":

                # src_type not defined for this sample. Move on.
                if src not in self.sample_data["project_data"]:
                    continue

                if script_path == "..import..":
                    self.sample_data[trg] = self.sample_data["project_data"][src]
                    return
                if script_path == "..skip..":
                    return

                self.spec_script_name = self.jid_name_sep.join([self.step,self.name,self.sample_data["Title"],src])
                
                # This line should be left before every new script. It sees to local issues.
                # Use the dir it returns as the base_dir for this step.
                use_dir = self.local_start(self.base_dir)
                
                fq_fn = ".".join([self.sample_data["Title"], src, self.file_tag,ext])          #The filename containing the end result. Used both in script and to set reads in $sample_params

                # Defining script:
                self.script = ""
                self.script += script_path + " \\\n\t"
                # The following line concatenates all the files in the direction separated by a " "
                self.script += " ".join(self.sample_data["project_data"][src]) 
                self.script += " \\\n\t"
                if self.params["pipe"][scope_ind]:  #"pipe" in self.params:
                    self.script += "| {pipe} \\\n\t".format(pipe = self.params["pipe"][scope_ind])
                self.script += "> %s%s \n\n"  % (use_dir, fq_fn)

                # Move all files from temporary local dir to permanent base_dir
                self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)

                # Store file in active file for sample:
                self.sample_data[trg] = self.base_dir + fq_fn
                
                self.stamp_file(self.sample_data[trg])

                self.create_low_level_script()
