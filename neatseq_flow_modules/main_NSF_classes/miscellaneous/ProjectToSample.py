# -*- coding: UTF-8 -*-
""" 
``ProjectToSample`` :sup:`*`
-----------------------------------------------------------------
:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A utility module for moving project data to a sample, and back again.
Is useful when a module which works on sample data has to be executed on data in the project scope.

For instance, in the STAR 2 pass pipeline, the first stage involves aligning all reads to the reference in order to find splice junctions.
The reads can be merged into a project scope ``fastq.F`` and ``fastq.R`` slots, but all aligners take there reads from the sample scope!

This module overrides the sample list with a single sample containing the project slots (or a subset of the slots).
Then, the mapping modules will take the project-wide reads from the sample representing the project.

Recovering the old sample list is done by setting the ``direction`` parameter to ``smp2proj``.

See the STAR2pass workflow for the working example.

Usually, the module should be called twice, once in the ``proj2smp`` direction and the in the ``smp2proj`` direction.
Although it is possible to use the ``smp2proj`` to move data from sample ``sample_name`` to the project, it is better to do this operation with the ``manage_types`` module.

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "direction", "proj2smp|smp2proj", "Move project info to sample or vice versa"
    "type","",   "The types to operate on"
    "operation","cp|mv",   "Whether to move the slots or just copy them."
    "sample_name","",   "The name of the new sample to create or the sample to copy from. Defaults to project title"


.. Attention:: This moduel does **NOT** operate on the actual files! It only modifies internal file types index.

Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Moving from project to sample:

::

    ProjectToSample:
        module:     ProjectToSample
        base:       merge_table
        script_path:
        direction:  proj2smp
        # sample_name:    fromproj
        operation:  mv   # mv or cp
        type:       [fastq.F,fastq.R]


Copying from sample to project:

::

    SampleToProject:
        module:     ProjectToSample
        base:       STAR_map_proj
        script_path:
        direction:  smp2proj
        operation:  mv   # mv or cp
        type:       SJ.out.tab

Copying and moving from sample to project:
(Just for the example. Isn't necessarily practical)

::

    SampleToProject:
        module:     ProjectToSample
        base:       STAR_map_proj
        script_path:
        direction:  smp2proj
        operation:  [cp, mv, mv]   # mv or cp
        type:       [SJ.out.tab, fastq.F, fastq.R]


        
"""

import os
import sys
import re
from neatseq_flow.PLC_step import Step,AssertionExcept

from pprint import pprint as pp



__author__ = "Menachem Sklarz"
__version__ = "1.1.0"



class Step_ProjectToSample(Step):
    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        # self.file_tag = "merge"
        self.skip_scripts = True
        
        
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """

        # Dealing with 'scope'
        if "direction" not in self.params:
            raise AssertionExcept("You must pass a 'direction' param!")
        
        if self.params["direction"] not in ["proj2smp", "smp2proj"]:
            raise AssertionExcept("'direction' param must be 'proj2smp' or 'smp2proj'")

        if "sample_name" in self.params \
                and self.params["sample_name"] in self.sample_data["samples"]\
                and self.params["direction"] == "proj2smp":
            raise AssertionExcept("'sample_name' {sample_name} exists! To avoid overwriting, please use a different "
                                  "name".format(sample_name=self.params["sample_name"]))

        if "operation" not in self.params:
            sys.stderr.write("You did not pass an 'operation' parameter. Assuming it is 'mv'!")
        if not isinstance(self.params["operation"], list):
            self.params["operation"] = re.split(pattern="\s*,\s*", string=self.params["operation"])
        if any([x not in ["mv","cp"] for x in self.params["operation"]]):
            raise AssertionExcept("'operation' must be 'mv' or 'cp'!")

        if "type" not in self.params:
            raise AssertionExcept("You must pass a 'type' param!")            
        if not isinstance(self.params["type"], list):
            self.params["type"] = re.split(pattern="\s*,\s*", string=self.params["type"])
            
        # Convert all arguments into lists, if exist
        for extra_params in (set(["operation", "type"]) & set(self.params.keys())):
            if not isinstance(self.params[extra_params], list):
                self.params[extra_params] = re.split(pattern="\s*,\s*", string=self.params[extra_params])

        # Check all lists have len 1 or same length
        active_params = set(["operation","type"]) & set(self.params.keys())
        # Get those params with len>1 (i.e., lists)
        list_params = [x for x in active_params if len(self.params[x])>1]
        str_params =  [x for x in active_params if len(self.params[x])==1]
        if len(set([len(self.params[x]) for x in list_params])) > 1:
            raise AssertionExcept("More than one list with len>1 specified! (%s)" % ", ".join(list_params))
        if list_params:
            required_len = len(self.params[list_params[0]])
            for i in str_params:
                self.params[i] = self.params[i] * required_len

        if "sample_name" not in self.params or self.params["sample_name"] is None:
            sample_name = self.sample_data["Title"]
        else:
            sample_name = self.params["sample_name"]

        # If direction is proj2smp, create new sample list and stash old in history:
        if self.params["direction"] == "proj2smp":
            # Create dictionary for new sample:
            if sample_name not in self.sample_data["samples"]:
                self.sample_data[sample_name] = {}
            # Stash old list and create new sample:
            self.stash_sample_list(sample_name)

        # Performing operation per index of list "operation"
        for oper_ind in range(len(self.params["operation"])):
            
            operation = self.params["operation"][oper_ind]

            try:
                type_src = self.params["type"][oper_ind]
            except KeyError:
                raise AssertionExcept("Make sure you have 'type' defined.")

            if self.params["direction"] == "proj2smp":
                # Moving project data types to sample.
                # # Create dictionary for new sample:
                # if sample_name not in self.sample_data["samples"]:
                #     self.sample_data[sample_name] = {}
                #     print "in here"
                # # Stash old list and create new sample:
                # self.stash_sample_list(sample_name)
                # Copy type from project to new sample:
                self.sample_data[sample_name][type_src] = self.sample_data["project_data"][type_src]
                if operation == "mv":
                    del self.sample_data["project_data"][type_src]

                print(sample_name)
            else:   # if direction == "smp2proj":
                # Moving project data types to sample.
                # Check that sample exists:
                if sample_name not in self.sample_data["samples"]:
                    raise AssertionExcept("Sample {sample} does not exist, can't move to project".
                                          format(sample=sample_name))
                # Copy type from sample to project:
                self.sample_data["project_data"][type_src] = self.sample_data[sample_name][type_src]
                if operation == "mv":
                    del self.sample_data[sample_name][type_src]
                # Recover old sample
                try:
                    self.recover_sample_list()
                except KeyError:
                    sys.stderr.write("ProjectToSample in 'smp2proj' direction called before a 'proj2smp' instance"
                                     "not recovering sample history!")

    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        
        pass
      
    def build_scripts(self):
        
        return 
