# -*- coding: UTF-8 -*-
""" 
``manage_types`` :sup:`*`
-----------------------------------------------------------------
:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A module for managing file type without script creation.

Supports adding, deleting, copying and moving file types.

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "operation", "add|del|mv|cp", "The operation to perform on the file type dictionary"
    "scope","project|sample",   "The scope on which to perform the operation. For 'mv' and 'cp' this is the source scope"
    "type","",   "The type on which to perform the operation. For 'mv' and 'cp' this is the source type"
    "scope_trg","project|sample",   "The destination scope for 'mv' and 'cp' operations"
    "type_trg","",   "The destination type for 'mv' and 'cp' operations"
    "path","",   "For 'add' operation, the value to insert in the file type."
    
.. Note:: If you don't pass ``scope_trg`` or ``type_trg``, they will be assumed to be the same and ``scope`` and ``type``, respectively. However, you **MUST** pass at least one of them!

.. Attention:: The operations do **NOT** operate on the actual files! They only modify internal file types index.

.. Tip:: For ``del`` operations, you can pass a list of types to remove. 

Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    manage_types1:
        module:             manage_types
        base:               STAR_bld_ind
        script_path:        
        scope:              project
        operation:          mv
        type:               trinity.contigs
        type_trg:           trinity.contigs
        scope_trg:          sample

        
"""

import os
import sys
import re
from neatseq_flow.PLC_step import Step,AssertionExcept




__author__ = "Menachem Sklarz"
__version__ = "1.1.0"



class Step_manage_types(Step):
    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        # self.file_tag = "merge"
        self.skip_scripts = True
        
        
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        
        if "scope" not in self.params:
            raise AssertionExcept("You must pass a 'scope' param!")
        if self.params["scope"] not in ["sample", "project"]:
            raise AssertionExcept("'scope' param must be 'sample' or 'project'")


        if "operation" not in self.params:
            raise AssertionExcept("You must pass an 'operation' parameter!")
        if self.params["operation"] not in ["del","mv","cp","add"]:
            raise AssertionExcept("'operation' must be one of 'del','mv','cp' and 'add'!")
##  ---------------- del ----------------------
        if self.params["operation"] == "del":
        
            if "type" not in self.params:
                raise AssertionExcept("You must pass a 'type' param!")
            type2del_list = self.params["type"]
            
            
            if not isinstance(type2del_list, list):
                type2del_list = [type2del]
            
            
            if self.params["scope"] == "sample":
                for sample in self.sample_data["samples"]:
                    for type2del in type2del_list:
                        if type2del not in self.sample_data[sample]:
                            self.write_warning("type %s does not exist for project." % type2del)
                        else:
                            del self.sample_data[sample][type2del]
            elif self.params["scope"] == "project":
                for type2del in type2del_list:
                    if type2del not in self.sample_data:
                        self.write_warning("type %s does not exist for project." % type2del)
                    else:
                        del self.sample_data[type2del]
            else:
                pass
##  ---------------- mv and cp ----------------------
        elif self.params["operation"] == "mv" or self.params["operation"] == "cp" :
            if "type" not in self.params:
                raise AssertionExcept("You must pass a 'type' param!")
            if "type_trg" not in self.params and "scope_trg" not in self.params:
                raise AssertionExcept("You must pass 'type_trg' or 'scope_trg', or both!")
            if "type_trg" not in self.params:
                self.write_warning("You did not specify 'type_trg'. Assuming same as 'type'.")
                self.params["type_trg"] = self.params["type"]
                # raise AssertionExcept("You must pass a 'type_trg' param!")
            if "scope_trg" not in self.params:
                self.write_warning("You did not specify 'scope_trg'. Assuming same as 'scope'.")
                self.params["scope_trg"] = self.params["scope"]
            if self.params["scope_trg"] not in ["sample","project"]:
                raise AssertionExcept("'scope_trg' must be either 'sample' or 'project'")
            
            type_src = self.params["type"]
            type_trg = self.params["type_trg"]
            
            if self.params["scope"] == "sample":
                if self.params["scope_trg"] == "sample":
                    for sample in self.sample_data["samples"]:
                        if type_src not in self.sample_data[sample]:
                            raise AssertionExcept("'type' does not exist in sample!",sample)
                        if type_trg in self.sample_data[sample]:
                            self.write_warning("type %s exists in sample. Overwriting!" % type_trg, sample)
                        self.sample_data[sample][type_trg] = self.sample_data[sample][type_src]
                        if self.params["operation"] == "mv":
                            del self.sample_data[sample][type_src]

                else:  # scope_trg=project
                    self.write_warning("Copying type from sample to project should be done with caution!")
                    if type_trg in self.sample_data:
                        self.write_warning("type %s exists in project. Overwriting!" % type_trg)
                    for sample in self.sample_data["samples"]:
                        self.sample_data[type_trg] = self.sample_data[sample][type_src]
                        if self.params["operation"] == "mv":
                            del self.sample_data[sample][type_src]
                        
            elif self.params["scope"] == "project":
                if type_src not in self.sample_data:
                    raise AssertionExcept("'type' does not exist in project!")
                if self.params["scope_trg"] == "sample":
                    for sample in self.sample_data["samples"]:
                        if type_trg in self.sample_data[sample]:
                            self.write_warning("type %s exists in sample. Overwriting!" % type_trg, sample)
                        self.sample_data[sample][type_trg] = self.sample_data[type_src]
                    if self.params["operation"] == "mv":
                        del self.sample_data[type_src]
                else:  # scope_trg=project
                    if type_trg in self.sample_data:
                        self.write_warning("type %s exists in project. Overwriting!" % type_trg, sample)
                    self.sample_data[type_trg] = self.sample_data[type_src]
                    if self.params["operation"] == "mv":
                        del self.sample_data[type_src]
            else:
                pass



##  ---------------- add ----------------------
        elif self.params["operation"] == "add":
            if "type" not in self.params:
                raise AssertionExcept("You must pass a 'type' param!")
            if "path" not in self.params:
                raise AssertionExcept("You must pass a 'path' param!")
            
            
            type2add = self.params["type"]
            
            if self.params["scope"] == "sample":
                for sample in self.sample_data["samples"]:
                    if type2add in self.sample_data[sample]:
                        self.write_warning("type %s exists in sample. Overwriting!" % type2add, sample)
                    self.sample_data[sample][type2add] = self.params["path"]
            elif self.params["scope"] == "project":
                if type2add in self.sample_data:
                    self.write_warning("type %s exists in project. Overwriting!" % type2add)
                self.sample_data[type2add] = self.params["path"]
            else:
                pass
            
        
        else:
            pass
        
        
        
    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        
        pass
      
    def build_scripts(self):
        
        return 
        # if 
        

                    