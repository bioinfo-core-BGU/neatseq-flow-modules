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
    

.. Attention:: The operations do **NOT** operate on the actual files! They only modify internal file types index.

.. Tip:: You can combine several operations in one module instance, by passing lists to the parameters in the table above. All lists should be of the same length, or of length 1 (`i.e.` plain strings). Plain strings will be extrapolated to all operations. `e.g.`, to delete one file type and add another, both at the project scope, pass [del,add] to the 'operation' parameter, and 'project' to the 'scope' parameter. The 'path' can also be a plain string. It will be extrapolated to 'del', as well, but will be ignored by it. See example lines below.

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

    manage_types1:
        module:             manage_types
        base:               trinity1
        script_path:        
        scope:              - project
                            - sample
                            - sample
                            - project
        operation:          - mv
                            - del
                            - cp
                            - add
        type:               - fasta.nucl
                            - fasta.nucl
                            - fastq.F
                            - bam
        type_trg:           [transcripts.nucl, None ,fastq.main, None]
        scope_trg:          sample
        path:               /path/to/mapping.bam   
        
"""

import os
import sys
import re
from neatseq_flow.PLC_step import Step,AssertionExcept

from pprint import pprint as pp



__author__ = "Menachem Sklarz"
__version__ = "1.6.0"


class Step_manage_types(Step):
    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        # self.file_tag = "merge"
        self.skip_scripts = True

    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """

        # Dealing with 'scope'
        if "scope" not in self.params:
            raise AssertionExcept("You must pass a 'scope' param!")
        
        if not isinstance(self.params["scope"], list):
            self.params["scope"] = re.split(pattern="\s*,\s*", string=self.params["scope"])

        if any([x not in ["sample", "project"] for x in self.params["scope"]]):
            raise AssertionExcept("'scope' param must be 'sample' or 'project'")

        if "operation" not in self.params:
            raise AssertionExcept("You must pass an 'operation' parameter!")
        if not isinstance(self.params["operation"], list):
            self.params["operation"] = re.split(pattern="\s*,\s*", string=self.params["operation"])
        if any([x not in ["del","mv","cp","add"] for x in self.params["operation"]]):
            raise AssertionExcept("'operation' must be one of 'del','mv','cp' and 'add'!")

        if "type" not in self.params:
            raise AssertionExcept("You must pass a 'type' param!")            
        if not isinstance(self.params["type"], list):
            self.params["type"] = re.split(pattern="\s*,\s*", string=self.params["type"])
            
        # Convert all arguments into lists, if exist
        for extra_params in (set(["scope_trg", "type_trg", "path"]) & set(self.params.keys())):
            if not isinstance(self.params[extra_params], list):
                self.params[extra_params] = re.split(pattern="\s*,\s*", string=self.params[extra_params])

        # Check all lists have len 1 or same length
        active_params = set(["scope","operation","type","scope_trg", "type_trg", "path"]) & set(self.params.keys())
        # Get those params with len>1 (i.e., lists)
        list_params = [x for x in active_params if len(self.params[x])>1]
        str_params =  [x for x in active_params if len(self.params[x])==1]
        if len(set([len(self.params[x]) for x in list_params])) > 1:
            raise AssertionExcept("More than one list with len>1 specified! (%s)" % ", ".join(list_params))
        # Extend all len-1 lists to required_len:
        if list_params:
            required_len = len(self.params[list_params[0]])
            for i in str_params:
                self.params[i] = self.params[i] * required_len

        # Now all lists are same length. Can perform operation per index of list "operation"
        for oper_ind in range(len(self.params["operation"])):
            
            operation = self.params["operation"][oper_ind]

#  ---------------- del ----------------------
            if operation == "del":
            
                type = self.params["type"][oper_ind]
                scope_src = self.params["scope"][oper_ind]

                if scope_src == "sample":
                    for sample in self.sample_data["samples"]:
                        if type not in self.sample_data[sample]:
                            self.write_warning("type %s does not exist for project." % type)
                        else:
                            del self.sample_data[sample][type]
                elif scope_src == "project":
                    if type not in self.sample_data["project_data"]:
                        self.write_warning("type %s does not exist for project." % type)
                    else:
                        del self.sample_data["project_data"][type]
                else:
                    pass

#  ---------------- mv and cp ----------------------

            elif operation in ["mv","cp"]:
            
                try:
                    type_src = self.params["type"][oper_ind]
                    scope_src = self.params["scope"][oper_ind]
                    type_trg = self.params["type_trg"][oper_ind]
                    scope_trg = self.params["scope_trg"][oper_ind]
                except KeyError:
                    raise AssertionExcept("'operation' is '%s'. Make sure you have 'type_trg' and 'scope_trg' defined." % operation)
                    
                if scope_src == "sample":
                    if scope_trg == "sample":
                        for sample in self.sample_data["samples"]:
                            if type_src not in self.sample_data[sample]:
                                raise AssertionExcept("'type' does not exist in sample!",sample)
                            if type_trg in self.sample_data[sample]:
                                self.write_warning("type %s exists in sample. Overwriting!" % type_trg, sample)
                            self.sample_data[sample][type_trg] = self.sample_data[sample][type_src]
                            if operation == "mv":
                                del self.sample_data[sample][type_src]

                    else:  # scope_trg=project
                        self.write_warning("Copying type from sample to project should be done with caution!")
                        if type_trg in self.sample_data["project_data"]:
                            self.write_warning("type %s exists in project. Overwriting!" % type_trg)
                        for sample in self.sample_data["samples"]:
                            try:
                                self.sample_data["project_data"][type_trg] = self.sample_data[sample][type_src]
                                if operation == "mv":
                                    del self.sample_data[sample][type_src]
                            except KeyError:
                                raise AssertionExcept("type '{src}' does not exist in sample!".format(src=type_src),
                                                      sample)
                            
                elif scope_src == "project":
                    if type_src not in self.sample_data["project_data"]:
                        raise AssertionExcept("'type' does not exist in project!")
                    if scope_trg == "sample":
                        for sample in self.sample_data["samples"]:
                            if type_trg in self.sample_data[sample]:
                                self.write_warning("type %s exists in sample. Overwriting!" % type_trg, sample)
                            self.sample_data[sample][type_trg] = self.sample_data["project_data"][type_src]
                        if operation == "mv":
                            del self.sample_data["project_data"][type_src]
                    else:  # scope_trg=project
                        if type_trg in self.sample_data["project_data"]:
                            self.write_warning("type %s exists in project. Overwriting!" % type_trg, sample)
                        self.sample_data["project_data"][type_trg] = self.sample_data["project_data"][type_src]
                        if operation == "mv":
                            del self.sample_data["project_data"][type_src]
                else:
                    pass

#  ---------------- add ----------------------
            elif operation == "add":

                try:
                    type = self.params["type"][oper_ind]
                    scope = self.params["scope"][oper_ind]
                    path = self.params["path"][oper_ind]
                except KeyError:
                    raise AssertionExcept("'operation' is 'add'. Make sure you have 'type', 'scope' and 'path' defined.")
                    
                if scope == "sample":
                    for sample in self.sample_data["samples"]:
                        if type in self.sample_data[sample]:
                            self.write_warning("type '%s' exists in sample. Overwriting!" % type, sample)
                        self.sample_data[sample][type] = path
                elif scope == "project":
                    if type in self.sample_data["project_data"]:
                        self.write_warning("type '%s' exists in project. Overwriting!" % type)
                    self.sample_data["project_data"][type] = path
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

        

                    