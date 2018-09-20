# -*- coding: UTF-8 -*-
""" 
``pipe_generic`` 
----------------------------------------------------

:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

    
Requires:
~~~~~~~~~~~~~

    
Output:
~~~~~~~~~~~~~
    

Parameters that can be set        
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. .. csv-table:: 
..     :header: "Parameter", "Values", "Comments"
..     :widths: 5,10,10
..     
..     "reference", "", "A block including 'path' or 'scope', 'type' and optionally 'msh'"
..     "query", "", "A block including 'scope' (sample, project or all_samples), 'type' and optionally 'msh'"


Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


"""

import os
import sys
import re
from pprint import pprint as pp
from copy import *
from neatseq_flow.PLC_step import Step, AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.1.0"


class Step_pipe_generic(Step):
    
    def step_specific_init(self):
        """ Called on intiation
            Good place for parameter testing.
            Wrong place for sample data testing
        """
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        self.file_tag = ".mash.tbl"

        # For some reason, yaml duplicates newlines in literal multiline string.
        # This removes them:
        self.params["script_path"] = self.params["script_path"].replace("\n\n","\n")
        # pp(dict(self.params))
        # sys.exit()

        # parms = re.findall(pattern="(\{\{[io]\:.*?\}\})", string=self.params["script_path"])
        # print parms
        # inputs = [re.search("\{[io]:(.*?)\}", name).group(1) for name in parms if
        #           re.search("\{([io]):.*?\}", name).group(1) == "i"]
        # outputs = [re.search("\{[io]:(.*?)\}", name).group(1) for name in parms if
        #            re.search("\{([io]):.*?\}", name).group(1) == "o"]
        # print inputs, outputs
        # """ TODO:
        # Tests:
        # 1. For each input inp:
        #     a. inp has 'type'
        #     b. if inp has scope, scope is sample of project
        # 2.  For each output outp:
        #     a. if outp has 'form', it is either 'file' or 'dir'
        #     b. if any input is of scope 'sample', set outp scope to sample, unless specified
        #     c. if outp has 'struct':
        #         i. parse struct
        #         ii. make sure all fields are legitimate: {sample}, {project}, {module}, {instance}
        #     d. 'store'
        #         i. if does not exist,
        # """
        # print self.params["script_path"]

        # Version2:
        # Get all user defined variables in string
        variables = re.findall(pattern="(\{\{.*?\}\})", string=self.params["script_path"])
        # Set scope. If any of the variables begins with 'sample:' or is only 'sample', set scope to sample
        # Otherwise, scope is project:
        if any([re.match(string=match,pattern="\{\{sample[\:\}]") for match in variables]):
            self.params["scope"] = "sample"
        else:
            self.params["scope"] = "project"

        # sys.exit()

    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        
        # TODO: Check that files exist for the analysis
        # print [self.params["output"][outp]["name"]
        #        for outp
        #        in self.params["output"]
        #        if self.params["output"][outp]["scope"]=="project"]
        # for outp in self.params["output"]:
        #     if self.params["output"][outp]["scope"]=="project"]:
        #         self.params["output"][outp]["name"] = self.format_script_path()
        # sys.exit()
        pass

    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        pass

    def build_scripts(self):
        """ This is the actual script building function
            
        """


        if self.params["scope"] =="sample":
            for sample in self.sample_data["samples"]:

                # Name of specific script:
                self.spec_script_name = self.set_spec_script_name(sample)
                self.script = ""

                # Make a dir for the current sample:
                sample_dir = self.make_folder_for_sample(sample)

                # This line should be left before every new script. It sees to local issues.
                # Use the dir it returns as the base_dir for this step.
                use_dir = self.local_start(sample_dir)

                if "output" in self.params:
                    self.params_output = deepcopy(self.params["output"])
                else:
                    self.params_output = {}
                for outp in self.params_output:
                    self.params_output[outp]["name"] = self.format_script_path(string=self.params_output[outp]["name"],
                                                                                  use_dir=use_dir,
                                                                                  sample=sample)

                self.script = self.format_script_path(self.params["script_path"],
                                                      use_dir=use_dir,
                                                      sample=sample)
                # # Try using function to include export (setenv) etc...

                # if "output" in self.params:
                for outp in self.params_output:
                    if self.params_output[outp]["scope"] == "sample":
                        self.sample_data[sample][outp] = self.format_script_path(string=self.params_output[outp]["name"],
                                                                                 use_dir=use_dir,
                                                                                 sample=sample)
                        self.stamp_file(self.sample_data[sample][outp])
                    else:
                        self.sample_data[outp] = self.format_script_path(
                            string=self.params_output[outp]["name"],
                            use_dir=use_dir)
                        self.stamp_file(self.sample_data[outp])

                # Wrapping up function. Leave these lines at the end of every iteration:
                self.local_finish(use_dir, sample_dir)
                self.create_low_level_script()
        else:

            # Name of specific script:
            self.spec_script_name = self.set_spec_script_name()
            self.script = ""

            # Make a dir for the current sample:
            sample_dir = self.base_dir  #self.make_folder_for_sample()

            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)

            if "output" in self.params:
                self.params_output = deepcopy(self.params["output"])
            else:
                self.params_output = {}
            for outp in self.params_output:
                self.params_output[outp]["name"] = self.format_script_path(string=self.params_output[outp]["name"],
                                                                           use_dir=use_dir)

            self.script = self.format_script_path(self.params["script_path"],
                                                  use_dir=use_dir)
            # # Try using function to include export (setenv) etc...

            # if "output" in self.params:
            for outp in self.params_output:
                if self.params_output[outp]["scope"] == "sample":
                    self.write_warning("Writing sample scope output for project scope script!")
                    for sample in self.sample_data["samples"]:
                        self.sample_data[sample][outp] = self.format_script_path(
                            string=self.params_output[outp]["name"],
                            use_dir=use_dir,
                            sample=sample)
                        self.stamp_file(self.sample_data[sample][outp])
                else:
                    self.sample_data[outp] = self.format_script_path(
                        string=self.params_output[outp]["name"],
                        use_dir=use_dir)
                    self.stamp_file(self.sample_data[outp])

            # Wrapping up function. Leave these lines at the end of every iteration:
            self.local_finish(use_dir, sample_dir)
            self.create_low_level_script()

    def format_script_path(self, string, use_dir=None, sample=None):
        """

        :return:
        """

        rawstring=string
        # Try using function to include export (setenv) etc...

        variables = list(set(re.findall(pattern="(\{\{.*?\}\})", string=rawstring)))

        for variable in variables:
            # ------------------------------
            if variable == "{{project}}":
                rawstring = re.sub(pattern=re.escape(variable),
                                     repl=self.sample_data["Title"],
                                     string=rawstring)
                continue
            # ------------------------------
            if variable == "{{dir}}":
                rawstring = rawstring.replace(variable, use_dir)
                continue
            # ------------------------------
            if variable == "{{base_dir}}":
                rawstring = rawstring.replace(variable, self.base_dir)
                continue
            # ------------------------------
            re_match = re.match("\{\{project\:(.*?)\}\}", variable)
            if re_match:
                # print re_match.group(1)
                # print re.escape(variable)
                try:
                    rawstring = re.sub(pattern=re.escape(variable),
                                         repl=("{!r}".format(self.sample_data[re_match.group(1)])).strip("'"),
                                         string=rawstring)
                except KeyError:
                    raise AssertionExcept(
                        "File type '{type}' not found in project scope".format(type=re_match.group(1)))
                continue
            # ------------------------------
            if variable == "{{sample}}":
                if not sample:
                    raise AssertionExcept("Trying to parse sample in project scope script!")
                rawstring = rawstring.replace(variable, sample)
                continue
            # ------------------------------
            re_match = re.match("\{\{sample\:(.*?)\}\}", variable)
            if re_match:
                if not sample:
                    raise AssertionExcept("Trying to parse sample in project scope script!")
                # print re_match.group(1)
                # print re.escape(variable)
                try:
                    rawstring = re.sub(pattern=re.escape(variable),
                                         repl=("{!r}".format(self.sample_data[sample][re_match.group(1)])).strip("'"),
                                         string=rawstring)
                except KeyError:
                    raise AssertionExcept("File type '{type}' not found in sample".format(type=re_match.group(1)),
                                          sample)
                continue
            # ------------------------------
            re_match = re.match("\{\{o\:(.*?)\}\}", variable)
            if re_match:
                try:
                    rawstring = re.sub(pattern=re.escape(variable),
                                         repl=("{!r}".format(self.params_output[re_match.group(1)]["name"])).strip("'"),
                                         string=rawstring)
                except KeyError:
                    raise AssertionExcept("Error embedding output '{var}'".format(var=variable), sample)
                continue
            # ------------------------------
            #  variable does not match any of the expected formats:
            raise AssertionExcept('Variable {var} in script_path not identified'.format(var=variable))

        return rawstring