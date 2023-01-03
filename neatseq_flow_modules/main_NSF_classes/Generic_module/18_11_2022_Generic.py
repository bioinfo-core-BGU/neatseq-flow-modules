# -*- coding: UTF-8 -*-
""" 
``Generic``
-----------------------

:Authors: Liron Levin
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

Short Description
~~~~~~~~~~~~~~~~~~~~~

A generic module that enables the user to design new modules that can handle most cases.

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* In this module the users define the required file types in the inputs section 

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* In this module the users define the output file types in the outputs section
* The scope of the output file types is determinant by the module scope

Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "scope",  "sample/project", "The scope of this module could be sample/project, the default is by sample"
    "shell","csh/bash","Type of shell [csh OR bash]. bash is the default, **only bash can be used in conda environment**"
    "inputs_last","","The inputs arguments will be at the end of the command"

Comments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    The order of the input/output arguments in the final command 
    will be according to the order of their appearance in the parameter file.
    The redirect arguments are always first. 



Example of usage and implementation of the generic module:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. attention:: 

    
    .. figure:: Generic_Module_Example.png
        :align: center
        :alt: Generic Module Example
        :figclass: align-left
    
    **A generic module is used to generate a BLAST database for each sample and a subsequent generic step queries each database with sequences from an external FASTA file. This example is a typical use of BLAST in many biological scenarios such as searching for virulence/resistance genes (whose sequences are in the external FASTA file) in bacterial genomes**


    **A.** Calling a generic module to generate a BLAST database (using makeblastdb) from each sample. This step can be used after (base:) any step that creates a nucleotide FASTA file (File_Type:  ``fasta.nucl``), e.g. after merge (if the raw files are in nucleotide FASTA format) or after a de novo assembly step. The location of the BLAST database for each sample is saved as a blast_db file type (File_Type: blast_db) for downstream use. **B.** Calling a generic module which performs a BLAST search (tblastn) of an external query protein fasta file (``-query`` :  path to query protein fasta file) against the previously generated BLAST data base per sample. This step can be used after the Make_BLAST_DB step (base: Make_BLAST_DB). The user can pass additional parameters directly to the used program in the redirects section (e.g. ``–dbtype``, ``–evalue``, ``-num_descriptions`` etc.). 


Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    Step_Name:                      # Name of this step
        module: Generic             # Name of module
        base:                       # Name of the step [or list of names] to run after [must be after steps that generates the inputs File_Types] 
        script_path:                # Main command for this module
        scope:                      # The scope of this module could be sample/project, the default is by sample
        shell:                      # Type of shell [csh OR bash]. bash is the default. only bash can be used in conda environment  
        arg_separator:              # The separator between the arguments and values [The default is space].
        inputs_last:                # The inputs arguments will be at the end of the command. [The default is inputs arguments at the beginning of the command]
        use_base_dir:               # Use the base step directory as the output for this step, it is possible to specify the base to use.
        cd:                         # Change current working directory to the output location.
        no_sample_dir:              # In Sample Scope: will NOT create a dedicated folder for each sample and the location of the base folder will be stored
                                    # in a project level 'base_dir' File_Type
        remove_subsamples:          # Will remove subsamples created by previous steps (split_fasta for example)
        subsamples_string:          # A string to identify a subsample, all subsample will start with this string. [default: 'subsample']
        inputs:                     # The inputs for this module
            STR:                    # Input argument, e.g. -i, --input [could be also 'empty1', 'empty2'.. for no input argument string]
                scope:              # The scope of this input argument could be sample/project
                                    # If the module scope is project and the argument scope is sample:
                                    # all the samples inputs File_Types of this argument will be listed as: [input argument] [File_Type(sample#)] e.g. -i sample1.bam -i sample2.bam ... 
                File_Type:          # The input File_Type could be any File_Type available from previous (in this branch) steps
                                    # It is possible to indicate more then one File_Type separated by comma 'fastq.F,fastq.R'
                base:               # From which previous step to take the input File_Type. The default is the current step.
                sep:                # If the module scope is project and the argument scope is sample:  
                                    #       All the samples inputs File_Types of this argument will be listed delimited by sep. e.g. [sep=,] -i sample1.bam,sample2.bam ... 
                                    # If more then one File_Type was specify the inputs File_Types of this argument will be listed delimited by sep.
                prefix:             # A prefix for this input argument file name
                suffix:             # A suffix for this input argument file name
                use_dirname:        # Use only the input Directory and add suffix for constant file name and prefix to add a string before the input Directory
                del:                # Delete the files in the input File_Type after the step ends [use to save space for large files you don't need downstream]
                                    # Will generate empty file with the same name and a suffix of _DELETED
                constant_value:     # use a constant value instead of "File_Type".
                                    # it is the same as the "redirects".
                                    # use when the order of inputs are important!!
                                    # use '{{sample_name}}' to be replace with the sample name (or project name in project scope)
                                    # using the constant_value option will override all other input arguments!!!!!!
        outputs:                    # The outputs for this module
            STR:                    # Output argument, e.g. -o, --out , the scope of the output arguments is determinant by the module scope
                                    # could be also 'empty1', 'empty2'.. for no output argument string OR 'No_run1', 'No_run2'.. for only entering the file information to output File_Type
                File_Type:          # The output File_Type could be any File_Type name for the current branch downstream work 
                                    # If the File_Type exists its content will be override for the current branch downstream work 
                prefix:             # A prefix for this output argument file name
                suffix:             # A suffix for this output argument file name
                                    # between prefix and suffix will be the sample name [in sample scope] or the project title [in project scope] 
                constant_file_name: # Use constant file name for this output argument [ignore prefix and suffix]
                                    # If empty [''] will enter the output directory location
                use_base_name:      # use only the base name of the output file [ignored if constant_file_name is used]
        copy_File_Types:            # Transferring information between File_Types
            STR:                    # Unique name for the transfer
                source:
                    File_Type:      # Copy the content of source File_Type to the target File_Type [copy from here]
                    scope:          # Copy the source File_Type From this scope [if not specified the default is sample]
                    base:           # The source step to copy the File_Type from (from previous steps). The default it the current step.
                    constant_value: # Use to transfer information outside of the 'File_Type' system to a File Type, will always be considered as project scope
                                    # Using the constant_value option will override all other source arguments!!!!!!
                target:
                    File_Type:      # Copy the content of source File_Type to the target File_Type [copy to here]
                    scope:          # Copy to the target File_Type in this scope [if not specified the default is sample]
        qsub_params:                # Parameters for qsub [number of cpus or memory to reserve etc ]
            STR: 
        redirects:                  # Parameters to pass directly to the command
            STR: 


"""


import os
import sys
import re
from functools import reduce
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Liron Levin"

class Step_Generic(Step):
    
    def step_specific_init(self):
        """ Called on intiation
            Good place for parameter testing.
            Wrong place for sample data testing
        """
        self.shell = get_File_Type_data(self.params,["shell"],"bash")
        self.project_del_script=[]
        if 'arg_separator' in self.params.keys():
            if self.params['arg_separator']==None:
                self.params['arg_separator']=""
        self.params['arg_separator'] = get_File_Type_data(self.params,["arg_separator"]," ")
        
        if 'subsamples_string' in self.params.keys():
            if self.params['subsamples_string']==None:
                self.params['subsamples_string']=""
        self.params['subsamples_string'] = get_File_Type_data(self.params,["subsamples_string"],"subsample")
        pass
        
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        if 'remove_subsamples' in self.params.keys():
            New_Sample_List = list()
            for sample in self.sample_data["samples"]:
                if sample.startswith(self.params['subsamples_string']):
                    self.sample_data.pop(sample)
                else:
                    New_Sample_List.append(sample)
            self.sample_data["samples"] = New_Sample_List
        
        if len(get_File_Type_data(self.params,["copy_File_Types"]))>0:
            for transfer in self.params["copy_File_Types"]:
                dif=set(["source","target",]).difference(list(self.params["copy_File_Types"][transfer].keys()))
                if len(dif)==0:
                    constant_value = get_File_Type_data(self.params["copy_File_Types"],[transfer,"source","constant_value"],None)
                    if constant_value == None:
                        scope_in      = get_File_Type_data(self.params["copy_File_Types"],[transfer,"source","scope"],"sample")
                        scope_out     = get_File_Type_data(self.params["copy_File_Types"],[transfer,"target","scope"],"sample")
                        File_Type_in  = get_File_Type_data(self.params["copy_File_Types"],[transfer,"source","File_Type"],None)
                        File_Type_out = get_File_Type_data(self.params["copy_File_Types"],[transfer,"target","File_Type"],None)
                        base          = get_File_Type_data(self.params["copy_File_Types"],[transfer,"source","base"],None)
                    else:
                        scope_in      = 'project'
                        scope_out     = get_File_Type_data(self.params["copy_File_Types"],[transfer,"target","scope"],"sample")
                        File_Type_in  = 'Not_Used'
                        File_Type_out = get_File_Type_data(self.params["copy_File_Types"],[transfer,"target","File_Type"],None)
                        base          = None
                    source_sample_data = self.sample_data
                    if base!=None:
                        if base in list(self.get_base_sample_data().keys()):
                            source_sample_data=self.get_base_sample_data()[base]
                        else:
                            raise AssertionExcept("The step name %s is not one of the previous steps of the %%s step" % base  % self.step )
                    else:
                        base=self.step
                        
                    if (File_Type_in and File_Type_out)!=None:
                        if "sample" in [scope_in , scope_out]:
                            for sample in self.sample_data["samples"]:
                                if scope_in == "sample":
                                    if File_Type_in in list(source_sample_data[sample].keys()):
                                        if scope_out=="sample":
                                            if not isinstance(source_sample_data[sample][File_Type_in],str):
                                                raise AssertionExcept("The File_Type %s in SAMPLE %%s is empty !!  in step %%%%s" % File_Type_in % sample % base)
                                            else:
                                                self.sample_data[sample][File_Type_out]=source_sample_data[sample][File_Type_in]
                                        else:
                                            #self.sample_data["project_data"][File_Type_out]=source_sample_data[sample][File_Type_in]
                                            raise AssertionExcept("It is not possible to transfer data from SAMPLE level to PROJECT level in transfer name: '%s' " % transfer)
                                    else:
                                        raise AssertionExcept("The File_Type %s is not found in the SAMPLE level  File_Types available are : %%s in step %%%%s" % File_Type_in % list(source_sample_data[sample].keys()) % base)
                                else:
                                    if constant_value == None:
                                        if File_Type_in in list(source_sample_data["project_data"].keys()):
                                            if not isinstance(source_sample_data["project_data"][File_Type_in],str):
                                                raise AssertionExcept("The File_Type %s in the PROJECT level is empty !!  in step %%s" % File_Type_in  % base)
                                            else:
                                                self.sample_data[sample][File_Type_out]=source_sample_data["project_data"][File_Type_in]
                                        else:
                                            raise AssertionExcept("The File_Type %s is not found in the PROJECT level \n\t "
                                                                  "File_Types available are : %s in step %s (base %s)" % \
                                                                    (File_Type_in,
                                                                     list(source_sample_data["project_data"].keys()),
                                                                     base))
                                    else:
                                        self.sample_data[sample][File_Type_out]=constant_value
                        else:
                            if constant_value == None:
                                if File_Type_in in list(source_sample_data["project_data"].keys()):
                                    if not isinstance(source_sample_data["project_data"][File_Type_in],str):
                                        raise AssertionExcept("The File_Type %s in the PROJECT level is empty !!  in step %%s" % File_Type_in  % base)
                                    else:
                                        self.sample_data["project_data"][File_Type_out]=source_sample_data["project_data"][File_Type_in]
                                else:
                                    raise AssertionExcept("The File_Type %s is not found in the PROJECT level \n\t File_Types available are : %%s in step %%%%s" % File_Type_in %  list(source_sample_data["project_data"].keys()) % base )
                            else:
                                self.sample_data["project_data"][File_Type_out]=constant_value
                    else:
                        if File_Type_in==None:
                            raise AssertionExcept("The following argument/s are missing or empty in the copy_File_Types section: %s" % "source File_Type")
                        if File_Type_out==None:
                            raise AssertionExcept("The following argument/s are missing or empty in the copy_File_Types section: %s" % "target File_Type")
                else:
                    raise AssertionExcept("The following argument/s are missing in the copy_File_Types section: %s" % list(dif))
        
        if "scope" in list(self.params.keys()):
            if "project" in self.params["scope"]:
                self.step_sample_initiation_byproject()
            else:
                self.step_sample_initiation_bysample()
        else:
            self.step_sample_initiation_bysample()
        
        
        
        if "use_base_dir" in list(self.params.keys()):
            if self.params['use_base_dir']==None:
                self.base_step_to_use = self.get_base_step_list()[0]
            else:
                depend_dir_dic = get_base_dirs_dict(self)
                print((list(depend_dir_dic.keys())))
                if self.params['use_base_dir'] in list(depend_dir_dic.keys()):
                    self.base_step_to_use = depend_dir_dic[self.params['use_base_dir']]
                else:
                    raise AssertionExcept("The step name %s is not one of the previous steps of the %%s step" % self.params['use_base_dir'] % self.step )

        
    pass
        
    def step_sample_initiation_bysample(self):
        """ A place to do initiation stages following setting of sample_data
            This set of tests is performed for sample-level
        """
        if len(get_File_Type_data(self.params,["inputs"]))>0:
            # Test if the input File_Types exists 
            for inputs in list(self.params["inputs"].keys()):
                
                if "constant_value" not in list(self.params["inputs"][inputs].keys()):
                    
                    base=get_File_Type_data(self.params["inputs"],[inputs,"base"],None)                
                    inputs_sample_data=self.sample_data
                    if base!=None:
                        if base in list(self.get_base_sample_data().keys()):
                            inputs_sample_data=self.get_base_sample_data()[base]
                        else:
                            raise AssertionExcept("The step name %s is not one of the previous steps of the %%s step" % base  % self.step )
                    else:
                        base=self.step
                    
                    
                    if get_File_Type_data(self.params["inputs"],[inputs,"File_Type"],None)==None: #Test if the user specify a File_Type for the input argument
                        raise AssertionExcept("You must specify a File_Type argument in the input parameter: %s " % inputs)
                    else:  #Test if the File_Type for the input argument exists
                        if get_File_Type_data(self.params["inputs"],[inputs,"scope"])=="project":
                            for File_Type in str(get_File_Type_data(self.params["inputs"],[inputs,"File_Type"])).replace("'",'').replace(" ",'').strip('[').strip(']').strip('"').split(','):
                                if File_Type not in list(inputs_sample_data["project_data"].keys()):
                                    raise AssertionExcept("The File_Type %s is not found in the PROJECT level \n\t File_Types available are : %%s in step %%%%s" % File_Type % list(inputs_sample_data["project_data"].keys()) % base )
                                elif not isinstance(inputs_sample_data["project_data"][File_Type],str):
                                    raise AssertionExcept("The File_Type %s in the PROJECT level is empty !!  in step %%s" % File_Type  % base)
                        else:
                            if "control" in self.params["scope"]:
                                try:
                                    self.sample_data["Controls"]
                                except KeyError:
                                    raise AssertionExcept("You must define sample:control pairs\n")
                                    
                                for samples in list(self.sample_data["Controls"].keys()):
                                    if get_File_Type_data(self.params["inputs"],[inputs,"scope"])=="control":
                                        sample = self.sample_data["Controls"][samples]
                                    else:
                                        sample = samples
                                    for File_Type in str(get_File_Type_data(self.params["inputs"],[inputs,"File_Type"])).replace("'",'').replace(" ",'').strip('[').strip(']').strip('"').split(','):
                                        if File_Type not in list(inputs_sample_data[sample].keys()):
                                            raise AssertionExcept("The File_Type %s is not found in the SAMPLE level [in sample name %%s] \n\t File_Types available are : %%%%s in step %%%%%%%%s" % File_Type % sample % list(inputs_sample_data[sample].keys()) % base )
                                        elif not isinstance(inputs_sample_data[sample][File_Type],str):
                                            raise AssertionExcept("The File_Type %s in SAMPLE %%s is empty !!  in step %%%%s" % File_Type % sample % base)
                                    
                            else:
                                if get_File_Type_data(self.params["inputs"],[inputs,"scope"])=="control":
                                    raise AssertionExcept("You can not use a control scope File_Type if the Step %s scope is not set to control " % base)
                                for sample in self.sample_data["samples"]:
                                    for File_Type in str(get_File_Type_data(self.params["inputs"],[inputs,"File_Type"])).replace("'",'').replace(" ",'').strip('[').strip(']').strip('"').split(','):
                                        if File_Type not in list(inputs_sample_data[sample].keys()):
                                            raise AssertionExcept("The File_Type %s is not found in the SAMPLE level [in sample name %%s] \n\t File_Types available are : %%%%s in step %%%%%%%%s" % File_Type % sample % list(inputs_sample_data[sample].keys()) % base )
                                        elif not isinstance(inputs_sample_data[sample][File_Type],str):
                                            raise AssertionExcept("The File_Type %s in SAMPLE %%s is empty !!  in step %%%%s" % File_Type % sample % base)
                                        
                    if "del" in list(self.params["inputs"][inputs].keys()):
                        self.write_warning("!!! The file/directory in the input File_Type %s in step %%s will be DELETED at the end of this step!!! " % get_File_Type_data(self.params["inputs"],[inputs,"File_Type"]) % base)
            
        if len(get_File_Type_data(self.params,["outputs"]))>0:
            # Test if the output File_Types
            for outputs in list(self.params["outputs"].keys()):
                if get_File_Type_data(self.params["outputs"],[outputs,"File_Type"],None)!=None: #Test if the user specify a File_Type for the output argument
                    if "control" in self.params["scope"]:
                        try:
                            self.sample_data["Controls"]
                        except KeyError:
                            raise AssertionExcept("You must define sample:control pairs\n")
                        for sample in list(self.sample_data["Controls"].keys()):
                            if get_File_Type_data(self.params["outputs"],[outputs,"File_Type"]) in list(self.sample_data[sample].keys()): #Test if the File_Type for the output argument exists
                                if self.sample_data[sample][get_File_Type_data(self.params["outputs"],[outputs,"File_Type"])]==None: #Test if the File_Type was already defined in the output arguments
                                    raise AssertionExcept("The output File_Type %s in the SAMPLE level was defined more the once !!! " % get_File_Type_data(self.params["outputs"],[outputs,"File_Type"]) )
                                else:
                                    self.write_warning("The output File_Type %s already exists in the SAMPLE level, its content will be override !!! " % get_File_Type_data(self.params["outputs"],[outputs,"File_Type"]) )
                            else: # If the File_Type dose not exists, will generate empty File_Type 
                                self.sample_data[sample][get_File_Type_data(self.params["outputs"],[outputs,"File_Type"])]=None
                    else:
                        for sample in self.sample_data["samples"]:
                            if get_File_Type_data(self.params["outputs"],[outputs,"File_Type"]) in list(self.sample_data[sample].keys()): #Test if the File_Type for the output argument exists
                                if self.sample_data[sample][get_File_Type_data(self.params["outputs"],[outputs,"File_Type"])]==None: #Test if the File_Type was already defined in the output arguments
                                    raise AssertionExcept("The output File_Type %s in the SAMPLE level was defined more the once !!! " % get_File_Type_data(self.params["outputs"],[outputs,"File_Type"]) )
                                else:
                                    self.write_warning("The output File_Type %s already exists in the SAMPLE level, its content will be override !!! " % get_File_Type_data(self.params["outputs"],[outputs,"File_Type"]) )
                            else: # If the File_Type dose not exists, will generate empty File_Type 
                                self.sample_data[sample][get_File_Type_data(self.params["outputs"],[outputs,"File_Type"])]=None
                if "del" in list(self.params["outputs"][outputs].keys()):
                    raise AssertionExcept("Output File_Types cannot be deleted")
        pass
        
    def step_sample_initiation_byproject(self):
        """ A place to do initiation stages following setting of sample_data
            This set of tests is performed for project-level 
        """
        if len(get_File_Type_data(self.params,["inputs"]))>0:
            # Test if the input File_Types exists 
            for inputs in list(self.params["inputs"].keys()):
                
                if "constant_value" not in list(self.params["inputs"][inputs].keys()):
                    base=get_File_Type_data(self.params["inputs"],[inputs,"base"],None)                
                    inputs_sample_data=self.sample_data
                    if base!=None:
                        if base in list(self.get_base_sample_data().keys()):
                            inputs_sample_data=self.get_base_sample_data()[base]
                        else:
                            raise AssertionExcept("The step name %s is not one of the previous steps of the %%s step" % base  % self.step )
                    else:
                        base=self.step    
                    
                    if get_File_Type_data(self.params["inputs"],[inputs,"File_Type"],None)==None: #Test if the user specify a File_Type for the input argument
                        raise AssertionExcept("You must specify a File_Type argument in the input parameter: %s " % inputs)
                    else:  #Test if the File_Type for the input argument exists
                        if get_File_Type_data(self.params["inputs"],[inputs,"scope"])=="project":
                            for File_Type in str(get_File_Type_data(self.params["inputs"],[inputs,"File_Type"])).replace("'",'').replace(" ",'').strip('[').strip(']').strip('"').split(','):
                                if File_Type not in list(inputs_sample_data["project_data"].keys()):
                                    raise AssertionExcept("The File_Type %s is not found in the PROJECT level \n\t File_Types available are : %%s in step %%%%s" % File_Type % list(inputs_sample_data["project_data"].keys()) % base)
                                elif not isinstance(inputs_sample_data["project_data"][File_Type],str):
                                    raise AssertionExcept("The File_Type %s in the PROJECT level is empty !!  in step %%s" % File_Type  % base)
                        else:
                            for sample in self.sample_data["samples"]:
                                for File_Type in str(get_File_Type_data(self.params["inputs"],[inputs,"File_Type"])).replace("'",'').replace(" ",'').strip('[').strip(']').strip('"').split(','):
                                    if File_Type not in list(inputs_sample_data[sample].keys()):
                                        raise AssertionExcept("The File_Type %s is not found in the SAMPLE level [in sample name %%s] \n\t File_Types available are : %%%%s in step %%%%%%%%s" % File_Type % sample % list(inputs_sample_data[sample].keys()) % base)
                                    elif not isinstance(inputs_sample_data[sample][File_Type],str):
                                        raise AssertionExcept("The File_Type %s in SAMPLE %%s is empty !!  in step %%%%s" % File_Type % sample % base)
                    if "del" in list(self.params["inputs"][inputs].keys()):
                        self.write_warning("!!! The file/directory in the input File_Type %s in step %%s will be DELETED at the end of this step!!! " % get_File_Type_data(self.params["inputs"],[inputs,"File_Type"]),base)
             
        if len(get_File_Type_data(self.params,["outputs"]))>0:
            # Test if the output File_Types
            for outputs in list(self.params["outputs"].keys()):
                if get_File_Type_data(self.params["outputs"],[outputs,"File_Type"],None)!=None: #Test if the user specify a File_Type for the output argument
                    if get_File_Type_data(self.params["outputs"],[outputs,"File_Type"]) in list(self.sample_data["project_data"].keys()): #Test if the File_Type for the output argument exists
                        if self.sample_data["project_data"][get_File_Type_data(self.params["outputs"],[outputs,"File_Type"])]==None: #Test if the File_Type was already defined in the output arguments
                            raise AssertionExcept("The output File_Type %s in the PROJECT level was defined more the once !!! " % get_File_Type_data(self.params["outputs"],[outputs,"File_Type"]) )
                        else:
                            self.write_warning("The output File_Type %s already exists in the PROJECT level, it's content will be override !!! " % get_File_Type_data(self.params["outputs"],[outputs,"File_Type"]) )
                    else: # If the File_Type dose not exists, will generate empty File_Type 
                        self.sample_data["project_data"][get_File_Type_data(self.params["outputs"],[outputs,"File_Type"])]=None
                if "del" in list(self.params["outputs"][outputs].keys()):
                    raise AssertionExcept("Output File_Types cannot be deleted")
        
        pass
        
        
    def create_spec_preliminary_script(self):
        """ Add script to run BEFORE all other steps
        """
        pass
    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        self.project_del_script=set(self.project_del_script)
        if len(self.project_del_script)>0:
            self.script=""
            for line in self.project_del_script:
                self.script+=line
        pass
        
    
    
    def build_scripts(self):
        """ This is the actual script building function
            
        """
        if "scope" in list(self.params.keys()):
            if "project" in self.params["scope"]:
                self.build_scripts_byproject()
            else:
                self.build_scripts_bysample()
        else:
            self.build_scripts_bysample()
            
        pass
    def build_scripts_bysample(self):
        """ Script building function for sample-level"""
        del_script=""
        # Each iteration must define the following class variables:
            # spec_script_name
            # script
            
        if "control" in self.params["scope"]:
            samples_list = list(self.sample_data["Controls"].keys())
        else:
            samples_list = self.sample_data["samples"]
            
        for sample in samples_list:      # Getting list of samples out of samples_hash
            if "control" in self.params["scope"]:
                Control_Sample = self.sample_data["Controls"][sample]
            # Name of specific script:
            self.spec_script_name = self.set_spec_script_name(sample)
            self.script = ""
            
            inputs_script = ""
            outputs_script = ""
            
            
            
            
            if 'use_base_dir' in list(self.params.keys()):
                if "no_sample_dir" in list(self.params.keys()):
                    sample_dir = self.base_dir
                    self.sample_data["project_data"]['base_dir'] = sample_dir
                else:
                    sample_dir = self.base_step_to_use.make_folder_for_sample(sample)
                self.base_dir = self.base_step_to_use.base_dir
            else:
                if "no_sample_dir" in list(self.params.keys()):
                    sample_dir = self.base_dir
                    self.sample_data["project_data"]['base_dir'] = sample_dir
                else:
                    # Make a dir for the current sample:
                    sample_dir = self.make_folder_for_sample(sample)
                
                
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)
            
            if 'cd' in list(self.params.keys()):
                self.script += "cd %s \n\n" % use_dir
            
            # Add the script constant args 
            self.script += self.get_script_const()
            if len(get_File_Type_data(self.params,["inputs"]))>0:
                del_script=""
                # Adds inputs files
                for inputs in list(self.params["inputs"].keys()): 
                    
                    if "constant_value" in list(self.params["inputs"][inputs].keys()):
                        value =get_File_Type_data(self.params["inputs"],[inputs,"constant_value"],'')
                        value = value.replace('{{sample_name}}',sample)
                        value = value.replace('{{project_name}}',self.sample_data["Title"])
                        if "control" in self.params["scope"]:
                            value = value.replace('{{control_name}}',Control_Sample)
                        if 'type' in list(self.sample_data[sample].keys()):
                            if len(get_File_Type_data(self.params["inputs"],[inputs,"sep"]))>0:
                                sep = get_File_Type_data(self.params["inputs"],[inputs,"sep"])
                            else:
                                sep = self.params['arg_separator']
                            value = value.replace('{{type}}',sep.join(map(str,self.sample_data[sample]['type'])))
                            if "control" in self.params["scope"]:
                                value = value.replace('{{type_control}}',sep.join(map(str,self.sample_data[Control_Sample]['type'])))
                        else:
                            value = value.replace('{{type}}','')
                            value = value.replace('{{type_control}}','')
                        if value!='':
                            if inputs.startswith("Empty".lower()):
                                inputs_script +="%s   \\\n\t"    % value
                            else:
                                inputs_script +="%s%%s%%%%s \\\n\t" % inputs \
                                                                    % self.params['arg_separator'] \
                                                                    % value
                    else:
                        base=get_File_Type_data(self.params["inputs"],[inputs,"base"],None)                
                        inputs_sample_data=self.sample_data
                        if base!=None:
                            if base in list(self.get_base_sample_data().keys()):
                                inputs_sample_data=self.get_base_sample_data()[base]
                            else:
                                raise AssertionExcept("The step name %s is not one of the previous steps of the %%s step" % base  % self.step )
                        else:
                            base=self.step    
                        
                        prefix       = get_File_Type_data(self.params["inputs"],[inputs,"prefix"])
                        suffix       = get_File_Type_data(self.params["inputs"],[inputs,"suffix"])
                        
                        prefix = prefix.replace('{{sample_name}}',sample)
                        prefix = prefix.replace('{{project_name}}',self.sample_data["Title"])
                        
                        if "control" in self.params["scope"]:
                            prefix = prefix.replace('{{control_name}}',Control_Sample)
                        
                        suffix = suffix.replace('{{sample_name}}',sample)
                        suffix = suffix.replace('{{project_name}}',self.sample_data["Title"])
                        
                        if "control" in self.params["scope"]:
                            suffix = suffix.replace('{{control_name}}',Control_Sample)
                        
                        if len(get_File_Type_data(self.params["inputs"],[inputs,"sep"]))>0:
                            sep=get_File_Type_data(self.params["inputs"],[inputs,"sep"])
                        else:
                            if inputs.startswith("Empty".lower()):
                                sep=" "
                            else:
                                sep=" \\\n\t"+inputs + self.params['arg_separator']
                        File_Type = ''
                        for File_Type_slot in str(get_File_Type_data(self.params["inputs"],[inputs,"File_Type"])).replace("'",'').replace(" ",'').strip('[').strip(']').strip('"').split(','):
                            if len(File_Type)>0:
                                File_Type+=sep
                            if get_File_Type_data(self.params["inputs"],[inputs,"scope"])=="project":
                                if 'use_dirname' in self.params["inputs"][inputs].keys():
                                    File_Type += prefix + os.path.join(os.path.dirname(inputs_sample_data["project_data"][File_Type_slot]) , (suffix).lstrip(os.sep) )
                                else:
                                    File_Type +=          os.path.join(os.path.dirname(inputs_sample_data["project_data"][File_Type_slot]) , (prefix + os.path.basename(inputs_sample_data["project_data"][File_Type_slot]) + suffix).lstrip(os.sep) )
                                
                                if ("del" in list(self.params["inputs"][inputs].keys()) ):
                                    if File_Type.endswith(os.sep):
                                        self.project_del_script.append("rm -rf %s*   \n\n"      % File_Type )
                                    else:
                                        self.project_del_script.append("rm -rf %s   \n\n"       % File_Type )
                                        self.project_del_script.append("echo > %s_DELETED \n\n" % File_Type.rstrip(os.sep) )
                            else:
                                if get_File_Type_data(self.params["inputs"],[inputs,"scope"])=="control":
                                    if 'use_dirname' in self.params["inputs"][inputs].keys():
                                        File_Type += prefix + os.path.join(os.path.dirname(     inputs_sample_data[Control_Sample][File_Type_slot]   ) , (suffix).lstrip(os.sep) )
                                    else:
                                        File_Type +=          os.path.join(os.path.dirname(     inputs_sample_data[Control_Sample][File_Type_slot]   ) , (prefix + os.path.basename( inputs_sample_data[Control_Sample][File_Type_slot]       ) + suffix).lstrip(os.sep) )
                                else:
                                    if 'use_dirname' in self.params["inputs"][inputs].keys():
                                        File_Type += prefix + os.path.join(os.path.dirname(     inputs_sample_data[sample][File_Type_slot]   ) , (suffix).lstrip(os.sep) )
                                    else:
                                        File_Type +=          os.path.join(os.path.dirname(     inputs_sample_data[sample][File_Type_slot]   ) , (prefix + os.path.basename( inputs_sample_data[sample][File_Type_slot]       ) + suffix).lstrip(os.sep) )
                                    
                                if ("del" in list(self.params["inputs"][inputs].keys()) ):
                                    if File_Type.endswith(os.sep):
                                        self.project_del_script.append("rm -rf %s*   \n\n"       % File_Type )
                                    else:
                                        del_script +="rm -rf %s   \n\n"        % File_Type
                                        del_script +="echo > %s_DELETED  \n\n" % File_Type.rstrip(os.sep)
                                        
                        if inputs.startswith("Empty".lower()):
                            inputs_script +="%s   \\\n\t"    % File_Type
                        else:
                            inputs_script +="%s%%s%%%%s \\\n\t" % inputs \
                                                                % self.params['arg_separator'] \
                                                                % File_Type
                
                    
            if len(get_File_Type_data(self.params,["outputs"]))>0:
                # Add output files
                for outputs in list(self.params["outputs"].keys()): 
                    # Define output filename  
                    if get_File_Type_data(self.params["outputs"],[outputs,"constant_file_name"])=="":
                        prefix=get_File_Type_data(self.params["outputs"],[outputs,"prefix"])
                        suffix=get_File_Type_data(self.params["outputs"],[outputs,"suffix"])
                        if "use_base_name" not in list(self.params["outputs"][outputs].keys()):                            
                            output_filename = "".join([use_dir ,prefix, sample , suffix])    
                            real_filename   = "".join([sample_dir ,prefix, sample , suffix]) 
                        else:
                            output_filename = "".join([prefix, sample , suffix])
                            real_filename   = "".join([prefix, sample , suffix])
                    else:
                        output_filename = "".join([use_dir ,get_File_Type_data(self.params["outputs"],[outputs,"constant_file_name"])])
                        real_filename   = "".join([sample_dir ,get_File_Type_data(self.params["outputs"],[outputs,"constant_file_name"])])
                    
                    real_filename   = real_filename.replace('{{project_name}}',self.sample_data["Title"])
                    output_filename = output_filename.replace('{{project_name}}',self.sample_data["Title"])
                    real_filename   = real_filename.replace('{{sample_name}}',sample)
                    output_filename = output_filename.replace('{{sample_name}}',sample)
                    
                    if "control" in self.params["scope"]:
                        real_filename   = real_filename.replace('{{control_name}}',Control_Sample)
                        output_filename = output_filename.replace('{{control_name}}',Control_Sample)
                    
                    File_Type=""
                    if outputs.startswith("No_run")==False:
                        if outputs.startswith("Empty".lower()):
                            outputs_script +="%s   \\\n\t"    % output_filename
                        else:
                            outputs_script +="%s%%s%%%%s \\\n\t" % outputs \
                                                            % self.params['arg_separator'] \
                                                            % output_filename
                    
                    #updating the output File_Types
                    if get_File_Type_data(self.params["outputs"],[outputs,"File_Type"],None)!=None:
                        File_Type=get_File_Type_data(self.params["outputs"],[outputs,"File_Type"])
                        # Save output file location in File_Type
                        self.sample_data[sample][File_Type]=(real_filename)
                        # Stamp the output file
                        self.stamp_file(self.sample_data[sample][File_Type])
            
            if "inputs_last" in list(self.params.keys()):
                self.script+=outputs_script
                self.script+=inputs_script
            else:
                self.script+=inputs_script
                self.script+=outputs_script
            
            self.script=self.script.rstrip("\\")
            self.script +="\n\n"
            # Delete the input file/directory if specified, before updating the output File_Types!!! 
            self.script += del_script
                
            # Wrapping up function. Leave these lines at the end of every iteration:
            self.local_finish(use_dir,sample_dir)       # Sees to copying local files to final destination (and other stuff)
            
            
            self.create_low_level_script()
                    
    def build_scripts_byproject(self):
        """ Script building function for project-level """
        
        del_script=""
        # Each iteration must define the following class variables:
        # spec_script_name
        # script
        
        # Name of specific script:
        self.spec_script_name = self.set_spec_script_name()
        self.script = ""
        
        inputs_script = ""
        outputs_script = ""
        
        
        if 'use_base_dir' in list(self.params.keys()):
            self.base_dir = self.base_step_to_use.base_dir
        
        # This line should be left before every new script. It sees to local issues.
        # Use the dir it returns as the base_dir for this step.
        use_dir = self.local_start(self.base_dir)
        
        if 'cd' in list(self.params.keys()):
            self.script += "cd %s \n\n" % use_dir
        
        # Add the script constant args 
        self.script += self.get_script_const()
        # Adds inputs files
        if len(get_File_Type_data(self.params,["inputs"]))>0:
            for inputs in list(self.params["inputs"].keys()):
                if "constant_value" in list(self.params["inputs"][inputs].keys()):
                    value =get_File_Type_data(self.params["inputs"],[inputs,"constant_value"],'')
                    if get_File_Type_data(self.params["inputs"],[inputs,"scope"])=="sample":
                        if len(get_File_Type_data(self.params["inputs"],[inputs,"sep"]))>0:
                            sep=get_File_Type_data(self.params["inputs"],[inputs,"sep"])
                        else:
                            sep = self.params['arg_separator']
                            
                        value =  value.replace('{{sample_name}}', sep.join(self.sample_data["samples"]))
                        
                        if len(get_File_Type_data(self.params["inputs"],[inputs,"sep"]))>0:
                            sep = get_File_Type_data(self.params["inputs"],[inputs,"sep"])
                        else:
                            sep = self.params['arg_separator']
                        type = []
                        for sample in self.sample_data["samples"]:
                            if 'type' in list(self.sample_data[sample].keys()):
                                type.append(self.sample_data[sample]['type'])
                        value = value.replace('{{type}}',sep.join(map(str,type)))
                        
                    value = value.replace('{{sample_name}}',self.sample_data["Title"])
                    value = value.replace('{{project_name}}',self.sample_data["Title"])
                    if 'type' in list(self.sample_data["project_data"].keys()):
                        value = value.replace('{{type}}',self.sample_data["project_data"]['type'])
                    if value!='':
                        if inputs.startswith("Empty".lower()):
                            inputs_script +="%s   \\\n\t"    % value
                        else:
                            inputs_script +="%s%%s%%%%s \\\n\t" % inputs \
                                                                % self.params['arg_separator'] \
                                                                % value
                else:
                    base=get_File_Type_data(self.params["inputs"],[inputs,"base"],None)                
                    inputs_sample_data=self.sample_data
                    if base!=None:
                        if base in list(self.get_base_sample_data().keys()):
                            inputs_sample_data=self.get_base_sample_data()[base]
                        else:
                            raise AssertionExcept("The step name %s is not one of the previous steps of the %%s step" % base  % self.step )
                    else:
                        base=self.step
                    
                    # prefix    = get_File_Type_data(self.params["inputs"],[inputs,"prefix"])
                    # suffix    = get_File_Type_data(self.params["inputs"],[inputs,"suffix"])
                    
                    File_Type = ""
                    
                    if len(get_File_Type_data(self.params["inputs"],[inputs,"sep"]))>0:
                        sep=get_File_Type_data(self.params["inputs"],[inputs,"sep"])
                    else:
                        if inputs.startswith("Empty".lower()):
                            sep=" "
                        else:
                            sep=" \\\n\t"+inputs + self.params['arg_separator']
                    
                    if get_File_Type_data(self.params["inputs"],[inputs,"scope"])=="project":
                        
                        prefix    = get_File_Type_data(self.params["inputs"],[inputs,"prefix"])
                        suffix    = get_File_Type_data(self.params["inputs"],[inputs,"suffix"])
                        
                        prefix = prefix.replace('{{sample_name}}',self.sample_data["Title"])
                        prefix = prefix.replace('{{project_name}}',self.sample_data["Title"])
                        
                        suffix = suffix.replace('{{sample_name}}',self.sample_data["Title"])
                        suffix = suffix.replace('{{project_name}}',self.sample_data["Title"])
                        
                        for File_Type_slot in str(get_File_Type_data(self.params["inputs"],[inputs,"File_Type"])).replace("'",'').replace(" ",'').strip('[').strip(']').strip('"').split(','):
                            if len(File_Type)>0:
                                File_Type+=sep
                            
                            if 'use_dirname' in self.params["inputs"][inputs].keys():
                                File_Type+=     prefix + os.path.join(os.path.dirname(inputs_sample_data["project_data"][File_Type_slot])  , (suffix).lstrip(os.sep) )
                            else:
                                File_Type+=              os.path.join(os.path.dirname(inputs_sample_data["project_data"][File_Type_slot])  , ( prefix + os.path.basename(inputs_sample_data["project_data"][File_Type_slot]) + suffix).lstrip(os.sep) )
                    else:
                        for sample in self.sample_data["samples"]:
                            
                            prefix    = get_File_Type_data(self.params["inputs"],[inputs,"prefix"])
                            suffix    = get_File_Type_data(self.params["inputs"],[inputs,"suffix"])
                            
                            prefix = prefix.replace('{{sample_name}}',sample)
                            prefix = prefix.replace('{{project_name}}',self.sample_data["Title"])
                            
                            suffix = suffix.replace('{{sample_name}}',sample)
                            suffix = suffix.replace('{{project_name}}',self.sample_data["Title"])
                            
                            for File_Type_slot in str(get_File_Type_data(self.params["inputs"],[inputs,"File_Type"])).replace("'",'').replace(" ",'').strip('[').strip(']').strip('"').split(','):
                                if len(File_Type)>0:
                                    File_Type+=sep
                                if 'use_dirname' in self.params["inputs"][inputs].keys():
                                    File_Type+=   prefix + os.path.join(os.path.dirname(        inputs_sample_data[sample][File_Type_slot])  , (suffix).lstrip(os.sep) )
                                else:
                                    File_Type+=            os.path.join(os.path.dirname(        inputs_sample_data[sample][File_Type_slot])  , ( prefix + os.path.basename(        inputs_sample_data[sample][File_Type_slot]) + suffix).lstrip(os.sep) )

                    if inputs.startswith("Empty".lower()):
                        inputs_script +="%s   \\\n\t"    % File_Type
                    else:
                        inputs_script +="%s%%s%%%%s \\\n\t" % inputs \
                                                                % self.params['arg_separator'] \
                                                                % File_Type
            
            
        
        if len(get_File_Type_data(self.params,["inputs"]))>0:
            # Generating delete script for input File_Types if specified 
            del_script=""
            for inputs in list(self.params["inputs"].keys()): 
                if "constant_value" not in list(self.params["inputs"][inputs].keys()):
                    if "del" in list(self.params["inputs"][inputs].keys()):
                        
                        base=get_File_Type_data(self.params["inputs"],[inputs,"base"],None)                
                        inputs_sample_data=self.sample_data
                        if base!=None:
                            if base in list(self.get_base_sample_data().keys()):
                                inputs_sample_data=self.get_base_sample_data()[base]
                            else:
                                raise AssertionExcept("The step name %s is not one of the previous steps of the %%s step" % base  % self.step )
                        else:
                            base=self.step
                        
                        prefix   = get_File_Type_data(self.params["inputs"],[inputs,"prefix"])
                        suffix   = get_File_Type_data(self.params["inputs"],[inputs,"suffix"])
                        
                        for File_Type_slot in str(get_File_Type_data(self.params["inputs"],[inputs,"File_Type"])).replace("'",'').replace(" ",'').strip('[').strip(']').strip('"').split(','):
                            if get_File_Type_data(self.params["inputs"],[inputs,"scope"])=="project":
                                if 'use_dirname' in self.params["inputs"][inputs].keys():
                                    file2delete = prefix + os.path.join(os.path.dirname(inputs_sample_data["project_data"][File_Type_slot])  , (suffix).rstrip(os.sep) )
                                else:
                                    file2delete =          os.path.join(os.path.dirname(inputs_sample_data["project_data"][File_Type_slot])  , ( prefix + os.path.basename(inputs_sample_data["project_data"][File_Type_slot]) + suffix).rstrip(os.sep) )
                                
                                if file2delete.endswith(os.sep):
                                    self.project_del_script.append("rm -rf %s*   \n\n"      % file2delete )
                                else:
                                    self.project_del_script.append("rm -rf %s   \n\n"       % file2delete )
                                    self.project_del_script.append("echo > %s_DELETED \n\n" % file2delete.rstrip(os.path.sep))
                            else:
                                for sample in self.sample_data["samples"]:
                                    if 'use_dirname' in self.params["inputs"][inputs].keys():
                                        file2delete = prefix + os.path.join(os.path.dirname(     inputs_sample_data[sample][File_Type_slot]    )  , (suffix).rstrip(os.sep)  )
                                    else:
                                        file2delete =          os.path.join(os.path.dirname(     inputs_sample_data[sample][File_Type_slot]    )  , ( prefix + os.path.basename(     inputs_sample_data[sample][File_Type_slot]     ) + suffix).rstrip(os.sep)  )
                                    
                                    if file2delete.endswith(os.sep):
                                        del_script +="rm -rf %s*   \n\n"       % file2delete
                                    else:
                                        del_script +="rm -rf %s   \n\n"        % file2delete
                                        del_script +="echo > %s_DELETED  \n\n" % file2delete.rstrip(os.path.sep)

        if len(get_File_Type_data(self.params,["outputs"]))>0:
            # Add output files
            for outputs in list(self.params["outputs"].keys()):
                # Define output filename  
                if get_File_Type_data(self.params["outputs"],[outputs,"constant_file_name"])=="":
                    prefix=get_File_Type_data(self.params["outputs"],[outputs,"prefix"])
                    suffix=get_File_Type_data(self.params["outputs"],[outputs,"suffix"])
                    if "use_base_name" not in list(self.params["outputs"][outputs].keys()):                       
                        output_filename = "".join([use_dir ,prefix, self.sample_data["Title"] , suffix])
                        real_filename   = "".join([self.base_dir ,prefix, self.sample_data["Title"] , suffix])
                    else:
                        output_filename = "".join([prefix, self.sample_data["Title"] , suffix])
                        real_filename   = "".join([prefix, self.sample_data["Title"] , suffix])
                else:
                    output_filename = "".join([use_dir ,get_File_Type_data(self.params["outputs"],[outputs,"constant_file_name"])])            
                    real_filename   = "".join([self.base_dir ,get_File_Type_data(self.params["outputs"],[outputs,"constant_file_name"])])    
                
                real_filename   = real_filename.replace('{{project_name}}',self.sample_data["Title"])
                output_filename = output_filename.replace('{{project_name}}',self.sample_data["Title"])
                if outputs.startswith("No_run")==False:
                    if outputs.startswith("Empty".lower()):
                        outputs_script +="%s   \\\n\t"    % output_filename
                    else:
                        outputs_script +="%s%%s%%%%s \\\n\t" % outputs \
                                                            % self.params['arg_separator'] \
                                                            % output_filename

                #updating the output File_Types
                if get_File_Type_data(self.params["outputs"],[outputs,"File_Type"],None)!=None:
                    File_Type=get_File_Type_data(self.params["outputs"],[outputs,"File_Type"])
                    # Save output file location in File_Type
                    self.sample_data["project_data"][File_Type] = (real_filename)
                    # Stamp the output file
                    self.stamp_file(self.sample_data["project_data"][File_Type])
                    
        if "inputs_last" in list(self.params.keys()):
            self.script+=outputs_script
            self.script+=inputs_script
        else:
            self.script+=inputs_script
            self.script+=outputs_script
            
            
        self.script=self.script.rstrip("\\")
        self.script +="\n\n"
        # Delete the input file/directory if specified, before updating the output File_Types!!! 
        self.script += del_script
        # Wrapping up function. Leave these lines at the end of every iteration:
        self.local_finish(use_dir,self.base_dir)       # Sees to copying local files to final destination (and other stuff)
        
        
        self.create_low_level_script()




def get_File_Type_data(dic,category,default=""):
    try:
        res=reduce(dict.get, category, dic)
        if res==None:
            return default
        else:
            return res
    except:
        return default
        
        
def get_base_dirs_dict(step,depend_dir_dic={}):
    """ Recursion. Beware
    Returns list of base dirs of the base steps of 'step'
    """

    
    if step.base_step_list:
        for base_step in step.base_step_list:
            depend_dir_dic[base_step.name]=base_step
            get_base_dirs_dict(base_step,depend_dir_dic)
            
    return depend_dir_dic