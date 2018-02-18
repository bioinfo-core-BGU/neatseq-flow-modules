# -*- coding: UTF-8 -*-
"""
``merge`` :sup:`*`
-------------------

:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A module for merging and staging files from the sample file into **NeatSeq-Flow**. 

Can be used in two modes: 

**The basic mode**
    is used when the sample file includes only sample-scope entries or only project-scope entries, and the file types are a subset of the list of recognized types. In this case, all you need to supply is the ``script_path`` (note below on legitimate ``script_path`` values). If the files are project-scope, you also need to set ``scope`` to ``project``. 
    
    You have to make sure that all your files are zipped or unzipped, so that they can all be concatenated with the same shell command. *i.e.* with the basic mode **you cannot** have both zipped and unzipped files, or files zipped with ``zip`` and ``gzip``. 
    
**The advanced mode**
    is used when the basic mode can't be used. It enables merging files from sample and project scope, including files not in the recognized list. It also enables mixing zipped and unzipped files. 
    
    In this mode, you have to define 4 lists: ``script_path``, ``src``, ``trg``, ``scope`` and ``ext``. For each file type in the sample file, you should have an entry in the ``src`` list. The other lists should apply to the equivalent entry in ``src``. ``trg`` is the target file type for the merged files, ``script_path`` is the shell command to use to merge the source type, ``scope`` is the scope for which the source type is defined and ``ext`` is the suffix to append to the merged filenames. Strings are expanded to the length of ``src`` list, so if ``script_path`` is the same for all source types, it is enough to specify it once. 
    
    If ``src`` is not passed, it is extracted from the list of types in the sample file. If an unrecognised source type is used, **NeatSeq-Flow** will demand you define the ``trg`` list as well. 
    
    The advanced mode is experimental, and documentation will hopefully improve as we gain experience with it.

.. Note:: **Definition of ``script_path`` in the ``merge`` module**: ``script_path`` should be a shell program that receives a list of files and produces one single output file **to the standard error**. Examples of such programs are ``cat`` for text files and ``gzip -cd`` for gzipped files. Other types of compressed files should have such a command as well. 


Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For the basic mode:

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
        src:            [UR1, UR2]
        script_path:    [gzip -cd, cat]
        scope:          [sample,project]
        trg:            [unrecognised1, unrecognised2]
        ext:            [ur1,ur2]

        
Advanced mode, single source type that exists both in sample and project::

    merge1:
        module:         merge
        script_path:    cat 
        scope:          [sample,project]
        

"""

import os
import sys
import re
from neatseq_flow.PLC_step import Step,AssertionExcept

from  neatseq_flow.modules.global_defs import ZIPPED_EXTENSIONS, ARCHIVE_EXTENSIONS, KNOWN_FILE_EXTENSIONS

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
            
    
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        
######################################



        for param_list in ["script_path","src","trg","ext","pipe","scope"]:
            if param_list in self.params:
                if not isinstance(self.params[param_list], list):
                    self.params[param_list] = re.split(pattern="\s*,\s*", string=self.params[param_list])
        if "scope" in self.params:
            if any(map(lambda x: x not in ["sample", "project"], self.params["scope"])):
                raise AssertionExcept("'scope' param must be 'sample' or 'project'")
        
######################################
                

        src = set()
        # Get list of existing file types in samples file:
        for sample in self.sample_data["samples"]:
            src = src | set(self.sample_data[sample].keys())
        if "project_data" in self.sample_data:
            src = src | set(self.sample_data["project_data"].keys())
        if "src" not in self.params:
            src = list(src)
            if "type" in src: 
                src.remove("type")   # 'type' is the type of sample (PE, etc.)
            self.params["src"] = src
        else: # Check that all src's exist somewhere in sample_data:
            if not set(self.params["src"]) <= src:
                raise AssertionExcept("The following types in 'src' do not exist in samples file: %s. " % " ".join(list(set(self.params["src"]) - src)))
        

        
        # Create trg list
        if "trg" not in self.params:
            if not set(self.params["src"]) < set(self.default_src_trg_map.keys()):
                raise AssertionExcept("The following types in 'src' are not recognized: %s. Please use explicit 'trg' list in merge" % " ".join(list(set(self.params["src"]) - set(self.default_src_trg_map.keys()))))
            self.params["trg"] = [self.default_src_trg_map[src][0] for src in self.params["src"]]
            
        if "ext" not in self.params:
            if not set(self.params["src"]) < set(self.default_src_trg_map.keys()):
                raise AssertionExcept("The following types in 'src' are not recognized: %s. Please use explicit 'ext' list." % " ".join(list(set(self.params["src"]) - set(self.default_src_trg_map.keys()))))
            self.params["ext"] = [self.default_src_trg_map[src][1] if len(self.default_src_trg_map[src])>1 else src for src in self.params["src"]]

        if "scope" not in self.params:
            self.params["scope"] = ["sample" for src in self.params["src"]]
        
            
        # Check all lists have len 1 or same length
        active_params = set(["script_path","src","trg","ext","pipe","scope"]) & set(self.params.keys())
        
        # Get those params with len>1 (i.e., lists)
        list_params = filter(lambda x: len(self.params[x])>1, active_params)
        str_params =  filter(lambda x: len(self.params[x])==1, active_params)
        
        
        if len(set(map(lambda x: len(self.params[x]), list_params))) > 1:
            raise AssertionExcept("More than one list with len>1 specified! (%s)" % ", ".join(list_params))

        if list_params:
            required_len = len(self.params[list_params[0]])
            for i in str_params:
                self.params[i] = self.params[i] * required_len

                
        # Check 'src's exist in 'scope's:
        for src_ind in range(len(self.params["src"])):
            src = self.params["src"][src_ind]
            scope = self.params["scope"][src_ind]
            if scope=="sample":
                for sample in self.sample_data["samples"]:
                    if src not in self.sample_data[sample]:
                        raise AssertionExcept("Type '{src}' does not exist for sample '{smp}'!".format(src=src,smp=sample))
            elif scope=="project":
                if src not in self.sample_data["project_data"]:
                    raise AssertionExcept("Type '{src}' does not exist in project data!".format(src=src))
            else:
                pass
        
        # pp(self.params["src"])
        # pp(self.params["trg"])
        # pp(self.params["ext"])
        # pp(self.params["scope"])
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

            if scope == "sample":
                # Each iteration must define the following class variables:
                    # spec_script_name
                    # script
                for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
                    # General comment: If there is a parallel routine for each direction (forward, reverse), add this loop	
                    # if  in self.sample_data[sample].keys():

                    # for type_i in range(len(self.params["src"])):
                
                    self.script = ""
                    
                    # Get index of src. Will be used to extract equivalent trg, script_path and ext.
                    # type_i = self.params["src"].index(src_type)
                    # src_type = self.params["src"][type_i]
                    
                    # src_type not defined for this sample. Move on.
                    if src not in self.sample_data[sample]:
                        continue
                        
                    self.spec_script_name = "_".join([self.step,self.name,sample,src]) 
                    
                    # This line should be left before every new script. It sees to local issues.
                    # Use the dir it returns as the base_dir for this step.
                    use_dir = self.local_start(self.base_dir)
                    
                    fq_fn = ".".join([sample, src, self.file_tag,ext])          #The filename containing the end result. Used both in script and to set reads in $sample_params


                    self.script += script_path + " \\\n\t"
                    # The following line concatenates all the files in the direction separated by a " "
                    self.script += " ".join(self.sample_data[sample][src]) 
                    self.script += " \\\n\t"
                    if "pipe" in self.params:
                        self.script += "| {pipe} \\\n\t".format(pipe = self.params["pipe"][scope_ind])
                    self.script += "> %s%s \n\n"  % (use_dir, fq_fn)

                    # Move all files from temporary local dir to permanent base_dir
                    self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)

                    
                    # Store file in active file for sample:
                    self.sample_data[sample][trg] = self.base_dir + fq_fn
                    
                    self.stamp_file(self.sample_data[sample][trg])
                    
                    
                    self.create_low_level_script()
                    
            elif scope == "project":

                self.script = ""
                
                # Get index of src. Will be used to extract equivalent trg, script_path and ext.
                # type_i = self.params["src"].index(src_type)
                # src_type = self.params["src"][type_i]
                
                # src_type not defined for this sample. Move on.
                if src not in self.sample_data["project_data"]:
                    continue
                    
                self.spec_script_name = "_".join([self.step,self.name,self.sample_data["Title"],src])

                # This line should be left before every new script. It sees to local issues.
                # Use the dir it returns as the base_dir for this step.
                use_dir = self.local_start(self.base_dir)
                
                fq_fn = ".".join([self.sample_data["Title"], src, self.file_tag,ext])          #The filename containing the end result. Used both in script and to set reads in $sample_params


                self.script += script_path + " \\\n\t"
                # The following line concatenates all the files in the direction separated by a " "
                self.script += " ".join(self.sample_data["project_data"][src]) 
                self.script += " \\\n\t"
                if "pipe" in self.params:
                    self.script += "| {pipe} \\\n\t".format(pipe = self.params["pipe"][scope_ind])
                self.script += "> %s%s \n\n"  % (use_dir, fq_fn)

                # Move all files from temporary local dir to permanent base_dir
                self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)

                
                # Store file in active file for sample:
                self.sample_data[trg] = self.base_dir + fq_fn
                
                self.stamp_file(self.sample_data[trg])
                
                
                self.create_low_level_script()
                