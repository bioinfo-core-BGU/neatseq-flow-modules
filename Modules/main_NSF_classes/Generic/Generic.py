""" 
``Generic`` (Included in main NeatSeq-Flow repo)
-----------------------------------------------------------------

A class that defines a module for a generic module.
It enables the user to design new modules that can handle most cases.

Parameters file arguments:
~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    Generic_step:                   # Name of this step
        module: Generic             # Name of module
        base:                       # Name of the step [or list of names] to run after [mast be after steps that generates the inputs File_Types] 
        script_path:                # Main command for this module
        scope:                      # The scope of this module could be sample/project, the default is by sample
        shell:                      # Type of shell [csh OR bash]. csh is the default
        inputs:                     # The inputs for this module
            STR:                    # Input argument, could be also 'empty#'
                scope:              # The scope of this input argument could be sample/project
                                    # If the module scope is project and the argument scope is sample:
                                    # all the samples inputs File_Types of this argument will be listed as: argument File_Type(of sample #) 
                File_Type:          # The input File_Type could be any File_Type available from previous (in this branch) steps
                sep:                # If the module scope is project and the argument scope is sample:  
                                    # all the samples inputs File_Types of this argument will be listed delimited by sep
                del:                # Delete this File_Type after the step ends
                                    # Will generate empty file with the same name and a suffix of _DELETED
        outputs:                    # The outputs for this module
            STR:                    # Output argument, the scope of the output arguments is determinant by the module scope
                                    # could be also 'empty#' OR "No_run" for only entering the file information to output File_Type 
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
                    scope:          # From the scope of scope_in [if scope_in is not specified the default is sample]
                target:
                    File_Type:      # Copy the content of source File_Type to the target File_Type [copy to here]
                    scope:          # From the scope of scope_out [if scope_out is not specified the default is sample]
        qsub_params:                # Parameters for qsub [number of cpus or memory to reserve etc ]
            STR: 
        redirects:                  # Parameters to pass directly to the command
            STR: 





"""

import os
import sys
import re
from PLC_step import Step,AssertionExcept


__author__ = "Liron Levin"

class Step_Generic(Step):
    
    def step_specific_init(self):
        """ Called on intiation
            Good place for parameter testing.
            Wrong place for sample data testing
        """
        self.shell = get_File_Type_data(self.params,["shell"],"csh")
        self.project_del_script=[]
        pass
        
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        if len(get_File_Type_data(self.params,["copy_File_Types"]))>0:
            for transfer in self.params["copy_File_Types"]:
                dif=set(["source","target",]).difference(self.params["copy_File_Types"][transfer].keys())
                if len(dif)==0:
                    scope_in=get_File_Type_data(self.params["copy_File_Types"],[transfer,"source","scope"],"sample")
                    scope_out=get_File_Type_data(self.params["copy_File_Types"],[transfer,"target","scope"],"sample")
                    File_Type_in=get_File_Type_data(self.params["copy_File_Types"],[transfer,"source","File_Type"],None)
                    File_Type_out=get_File_Type_data(self.params["copy_File_Types"],[transfer,"target","File_Type"],None)
                    
                    if (File_Type_in and File_Type_out)!=None:
                        if "sample" in [scope_in , scope_out]:
                            for sample in self.sample_data["samples"]:
                                if scope_in =="sample":
                                    if File_Type_in in self.sample_data[sample].keys():
                                        if scope_out=="sample":
                                            self.sample_data[sample][File_Type_out]=self.sample_data[sample][File_Type_in]
                                        else:
                                            self.sample_data[File_Type_out]=self.sample_data[sample][File_Type_in]
                                    else:
                                        raise AssertionExcept("The File_Type %s is not found in the SAMPLE level \n\t File_Types available are : %%s" % File_Type_in % self.sample_data[sample].keys() )
                                else:
                                    if File_Type_in in self.sample_data.keys():
                                        self.sample_data[sample][File_Type_out]=self.sample_data[File_Type_in]
                                    else:
                                        raise AssertionExcept("The File_Type %s is not found in the PROJECT level \n\t File_Types available are : %%s" % File_Type_in % self.sample_data.keys() )
                        else:
                            if File_Type_in in self.sample_data.keys():
                                self.sample_data[File_Type_out]=self.sample_data[File_Type_in]
                            else:
                                raise AssertionExcept("The File_Type %s is not found in the PROJECT level \n\t File_Types available are : %%s" % File_Type_in % self.sample_data.keys() )
                    else:
                        if File_Type_in==None:
                            raise AssertionExcept("The following argument/s are missing or empty in the copy_File_Types section: %s" % "source File_Type")
                        if File_Type_out==None:
                            raise AssertionExcept("The following argument/s are missing or empty in the copy_File_Types section: %s" % "target File_Type")
                else:
                    raise AssertionExcept("The following argument/s are missing in the copy_File_Types section: %s" % list(dif))
        
        if "scope" in self.params.keys():
            if "project" in self.params["scope"]:
                self.step_sample_initiation_byproject()
            else:
                self.step_sample_initiation_bysample()
        else:
            self.step_sample_initiation_bysample()
        pass
        
    def step_sample_initiation_bysample(self):
        """ A place to do initiation stages following setting of sample_data
            This set of tests is performed for sample-level
        """
        if len(get_File_Type_data(self.params,["inputs"]))>0:
            # Test if the input File_Types exists 
            for inputs in self.params["inputs"].keys():
                if get_File_Type_data(self.params["inputs"],[inputs,"File_Type"],None)==None: #Test if the user specify a File_Type for the input argument
                    raise AssertionExcept("You mast specify a File_Type argument in the input parameter: %s " % inputs)
                else:  #Test if the File_Type for the input argument exists
                    if get_File_Type_data(self.params["inputs"],[inputs,"scope"])=="project":
                        if get_File_Type_data(self.params["inputs"],[inputs,"File_Type"]) not in self.sample_data.keys():
                            raise AssertionExcept("The File_Type %s is not found in the PROJECT level \n\t File_Types available are : %%s" % get_File_Type_data(self.params["inputs"],[inputs,"File_Type"]) % self.sample_data.keys() )
                    else:
                        for sample in self.sample_data["samples"]:
                            if get_File_Type_data(self.params["inputs"],[inputs,"File_Type"]) not in self.sample_data[sample].keys():                            
                                raise AssertionExcept("The File_Type %s is not found in the SAMPLE level [in sample name %%s] \n\t File_Types available are : %%%%s" % get_File_Type_data(self.params["inputs"],[inputs,"File_Type"]) % sample % self.sample_data[sample].keys())
                if "del" in self.params["inputs"][inputs].keys():
                    self.write_warning("!!! The file/directory in the input File_Type %s will be DELETED at the end of this step!!! " % get_File_Type_data(self.params["inputs"],[inputs,"File_Type"]))
        
        if len(get_File_Type_data(self.params,["outputs"]))>0:
            # Test if the output File_Types
            for outputs in self.params["outputs"].keys():
                if get_File_Type_data(self.params["outputs"],[outputs,"File_Type"],None)!=None: #Test if the user specify a File_Type for the output argument
                    for sample in self.sample_data["samples"]:
                        if get_File_Type_data(self.params["outputs"],[outputs,"File_Type"]) in self.sample_data[sample].keys(): #Test if the File_Type for the output argument exists
                            if self.sample_data[sample][get_File_Type_data(self.params["outputs"],[outputs,"File_Type"])]==None: #Test if the File_Type was already defined in the output arguments
                                raise AssertionExcept("The output File_Type %s in the SAMPLE level was defined more the once !!! " % get_File_Type_data(self.params["outputs"],[outputs,"File_Type"]) )
                            else:
                                self.write_warning("The output File_Type %s already exists in the SAMPLE level, its content will be override !!! " % get_File_Type_data(self.params["outputs"],[outputs,"File_Type"]) )
                        else: # If the File_Type dose not exists, will generate empty File_Type 
                            self.sample_data[sample][get_File_Type_data(self.params["outputs"],[outputs,"File_Type"])]=None
                if "del" in self.params["outputs"][outputs].keys():
                    raise AssertionExcept("Output File_Types cannot be deleted")
        pass
        
    def step_sample_initiation_byproject(self):
        """ A place to do initiation stages following setting of sample_data
            This set of tests is performed for project-level 
        """
        if len(get_File_Type_data(self.params,["inputs"]))>0:
            # Test if the input File_Types exists 
            for inputs in self.params["inputs"].keys():            
                if get_File_Type_data(self.params["inputs"],[inputs,"File_Type"],None)==None: #Test if the user specify a File_Type for the input argument
                    raise AssertionExcept("You mast specify a File_Type argument in the input parameter: %s " % inputs)
                else:  #Test if the File_Type for the input argument exists
                    if get_File_Type_data(self.params["inputs"],[inputs,"scope"])=="project":
                        if get_File_Type_data(self.params["inputs"],[inputs,"File_Type"]) not in self.sample_data.keys():
                            raise AssertionExcept("The File_Type %s is not found in the PROJECT level \n\t File_Types available are : %%s" % get_File_Type_data(self.params["inputs"],[inputs,"File_Type"]) % self.sample_data.keys())
                    else:
                        for sample in self.sample_data["samples"]:
                            if get_File_Type_data(self.params["inputs"],[inputs,"File_Type"]) not in self.sample_data[sample].keys():
                                raise AssertionExcept("The File_Type %s is not found in the SAMPLE level [in sample name %%s] \n\t File_Types available are : %%%%s" % get_File_Type_data(self.params["inputs"],[inputs,"File_Type"]) % sample % self.sample_data[sample].keys() )
                if "del" in self.params["inputs"][inputs].keys():
                    self.write_warning("!!! The file/directory in the input File_Type %s will be DELETED at the end of this step!!! " % get_File_Type_data(self.params["inputs"],[inputs,"File_Type"]))
        
        if len(get_File_Type_data(self.params,["outputs"]))>0:
            # Test if the output File_Types
            for outputs in self.params["outputs"].keys():
                if get_File_Type_data(self.params["outputs"],[outputs,"File_Type"],None)!=None: #Test if the user specify a File_Type for the output argument
                    if get_File_Type_data(self.params["outputs"],[outputs,"File_Type"]) in self.sample_data.keys(): #Test if the File_Type for the output argument exists
                        if self.sample_data[get_File_Type_data(self.params["outputs"],[outputs,"File_Type"])]==None: #Test if the File_Type was already defined in the output arguments
                            raise AssertionExcept("The output File_Type %s in the PROJECT level was defined more the once !!! " % get_File_Type_data(self.params["outputs"],[outputs,"File_Type"]) )
                        else:
                            self.write_warning("The output File_Type %s already exists in the PROJECT level, it's content will be override !!! " % get_File_Type_data(self.params["outputs"],[outputs,"File_Type"]) )
                    else: # If the File_Type dose not exists, will generate empty File_Type 
                        self.sample_data[get_File_Type_data(self.params["outputs"],[outputs,"File_Type"])]=None
                if "del" in self.params["outputs"][outputs].keys():
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
        if "scope" in self.params.keys():
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
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
            
            # Name of specific script:
            self.spec_script_name = "_".join([self.step,self.name,sample])
            self.script = ""

            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)
            
            # Add the script constant args 
            self.script += self.get_script_const()
            if len(get_File_Type_data(self.params,["inputs"]))>0:
                # Adds inputs files
                for inputs in self.params["inputs"].keys(): 
                    File_Type=""
                    if get_File_Type_data(self.params["inputs"],[inputs,"scope"])=="project":
                        File_Type=self.sample_data[get_File_Type_data(self.params["inputs"],[inputs,"File_Type"])]
                    else:
                        File_Type=self.sample_data[sample][get_File_Type_data(self.params["inputs"],[inputs,"File_Type"])]
                    if inputs.startswith("empty"):
                        self.script +="%s   \\\n\t"    % File_Type
                    else:
                        self.script +="%s  %%s \\\n\t" % inputs \
                                                       % File_Type 
            if len(get_File_Type_data(self.params,["inputs"]))>0:
                # Generating delete script for input File_Types if specified 
                del_script=""
                for inputs in self.params["inputs"].keys():                                       
                    if "del" in self.params["inputs"][inputs].keys():
                        if get_File_Type_data(self.params["inputs"],[inputs,"scope"])=="project":
                            self.project_del_script.append("rm -rf %s   \n\n" % self.sample_data[get_File_Type_data(self.params["inputs"],[inputs,"File_Type"])])
                            self.project_del_script.append("echo > %s_DELETED \n\n" % self.sample_data[get_File_Type_data(self.params["inputs"],[inputs,"File_Type"])].rstrip(os.path.sep))
                        else:
                            del_script +="rm -rf %s   \n\n" % self.sample_data[sample][get_File_Type_data(self.params["inputs"],[inputs,"File_Type"])]
                            del_script +="echo > %s_DELETED  \n\n" % self.sample_data[sample][get_File_Type_data(self.params["inputs"],[inputs,"File_Type"])].rstrip(os.path.sep)

            if len(get_File_Type_data(self.params,["outputs"]))>0:
                # Add output files
                for outputs in self.params["outputs"].keys(): 
                    # Define output filename  
                    if get_File_Type_data(self.params["outputs"],[outputs,"constant_file_name"])=="":
                        prefix=get_File_Type_data(self.params["outputs"],[outputs,"prefix"])
                        suffix=get_File_Type_data(self.params["outputs"],[outputs,"suffix"])
                        if "use_base_name" not in self.params["outputs"][outputs].keys():                            
                            output_filename = "".join([use_dir ,prefix, sample , suffix])    
                        else:
                            output_filename = "".join([prefix, sample , suffix])
                    else:
                        output_filename = "".join([use_dir ,get_File_Type_data(self.params["outputs"],[outputs,"constant_file_name"])])
                    File_Type=""
                    if outputs.startswith("No_run")==False:
                        if outputs.startswith("empty"):
                            self.script +="%s   \\\n\t"    % output_filename
                        else:
                            self.script +="%s  %%s \\\n\t" % outputs \
                                                           % output_filename
                    
                    #updating the output File_Types
                    if get_File_Type_data(self.params["outputs"],[outputs,"File_Type"],None)!=None:
                        File_Type=get_File_Type_data(self.params["outputs"],[outputs,"File_Type"])
                        # Save output file location in File_Type
                        self.sample_data[sample][File_Type]=(sample_dir + os.path.basename(output_filename))
                        # Stamp the output file
                        self.stamp_file(self.sample_data[sample][File_Type])
            
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
        self.spec_script_name = "_".join([self.step,self.name,self.sample_data["Title"]])
        self.script = ""

        
        # This line should be left before every new script. It sees to local issues.
        # Use the dir it returns as the base_dir for this step.
        use_dir = self.local_start(self.base_dir)
        
         # Add the script constant args 
        self.script += self.get_script_const()
        # Adds inputs files
        if len(get_File_Type_data(self.params,["inputs"]))>0:
            for inputs in self.params["inputs"].keys(): 
                File_Type=""
                if get_File_Type_data(self.params["inputs"],[inputs,"scope"])=="project":
                    File_Type=self.sample_data[get_File_Type_data(self.params["inputs"],[inputs,"File_Type"])]
                else:
                    if len(get_File_Type_data(self.params["inputs"],[inputs,"sep"]))>0:
                        sep=get_File_Type_data(self.params["inputs"],[inputs,"sep"])
                    else:
                        if inputs.startswith("empty"):
                            sep=" "
                        else:
                            sep=" \\\n\t"+inputs+"  "
                    for sample in self.sample_data["samples"]:
                        File_Type+=self.sample_data[sample][get_File_Type_data(self.params["inputs"],[inputs,"File_Type"])]
                        File_Type+=sep
                    File_Type=File_Type.strip(sep)
                
                if inputs.startswith("empty"):
                    self.script +="%s   \\\n\t"    % File_Type
                else:
                    self.script +="%s  %%s \\\n\t" % inputs \
                                                   % File_Type


        
        if len(get_File_Type_data(self.params,["inputs"]))>0:
            # Generating delete script for input File_Types if specified 
            del_script=""
            for inputs in self.params["inputs"].keys():                                       
                if "del" in self.params["inputs"][inputs].keys():
                    if get_File_Type_data(self.params["inputs"],[inputs,"scope"])=="project":
                        self.project_del_script.append("rm -rf %s   \n\n" % self.sample_data[get_File_Type_data(self.params["inputs"],[inputs,"File_Type"])])
                        self.project_del_script.append("echo > %s_DELETED \n\n" % self.sample_data[get_File_Type_data(self.params["inputs"],[inputs,"File_Type"])].rstrip(os.path.sep))
                    else:
                        for sample in self.sample_data["samples"]:
                            del_script +="rm -rf %s   \n\n" % self.sample_data[sample][get_File_Type_data(self.params["inputs"],[inputs,"File_Type"])]
                            del_script +="echo > %s_DELETED  \n\n" % self.sample_data[sample][get_File_Type_data(self.params["inputs"],[inputs,"File_Type"])].rstrip(os.path.sep)

        if len(get_File_Type_data(self.params,["outputs"]))>0:
            # Add output files
            for outputs in self.params["outputs"].keys():
                # Define output filename  
                if get_File_Type_data(self.params["outputs"],[outputs,"constant_file_name"])=="":
                    prefix=get_File_Type_data(self.params["outputs"],[outputs,"prefix"])
                    suffix=get_File_Type_data(self.params["outputs"],[outputs,"suffix"])
                    if "use_base_name" not in self.params["outputs"][outputs].keys():                       
                        output_filename = "".join([use_dir ,prefix, self.sample_data["Title"] , suffix])
                    else:
                        output_filename = "".join([prefix, self.sample_data["Title"] , suffix])
                else:
                    output_filename = "".join([use_dir ,get_File_Type_data(self.params["outputs"],[outputs,"constant_file_name"])])            
                if outputs.startswith("No_run")==False:
                    if outputs.startswith("empty"):
                        self.script +="%s   \\\n\t"    % output_filename
                    else:
                        self.script +="%s  %%s \\\n\t" % outputs \
                                                       % output_filename

                #updating the output File_Types
                if get_File_Type_data(self.params["outputs"],[outputs,"File_Type"],None)!=None:
                    File_Type=get_File_Type_data(self.params["outputs"],[outputs,"File_Type"])
                    # Save output file location in File_Type
                    self.sample_data[File_Type] = (self.base_dir + os.path.basename(output_filename))
                    # Stamp the output file
                    self.stamp_file(self.sample_data[File_Type])
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