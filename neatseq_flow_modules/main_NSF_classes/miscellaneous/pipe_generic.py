# -*- coding: UTF-8 -*-
""" 
``pipe_generic`` 
----------------------------------------------------

:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

Description
~~~~~~~~~~~~

This module enables executing any type of shell script, including pipes and multiple steps.
Files and directories are embedded in the script by describing the file or directory in a ``{{}}`` block, as follows:

1. File names:
***************

Include 4 ``:``-separated fields: (a) scope, (b) slot, (c) separator and (d) base. For example: ``{{sample:fastq.F:,:merge1}}``
Leave the field empty if you do not want to pass a value, e.g. ``{{sample:fastq.F}}``
The following examples cover most of the options:

.. csv-table::
    :header: "File descriptor", "Result"
    :widths: 5,10

    "``{{project:fasta.nucl}}``", "The ``fasta.nucl`` slot of the project"
    "``{{sample:fastq.F}}``", "The ``fastq.F`` slot of the sample"
    "``{{sample:fastq.F:,}}``", "A comma-separated list of the ``fastq.F`` slots of all samples"
    "``{{project}}``", "The project name"
    "``{{sample}}``", "The sample name"
    "``{{sample::,}}``", "A comma-separated list of sample names"
    "``{{sample:fastq.F:,:base}}``", "A comma-separated list of the ``fastq.F`` files of all samples, taken from the sample data of step ``base``."

.. Tip:: For a colon separate list of sample names or files, use the word 'colon' in the separator slot.

.. Note:: The separator field is ignored for project-scope slots.

If a sample-scope slot is used, in the inputs or the outputs, the scripts will be sample-scope scripts. Otherwise, one project-scope script will be produced.

2. Directories
+++++++++++++++++++

Can take one of two values:

.. csv-table::
    :header: "Dir descriptor", "Result"
    :widths: 5,10

    "``{{base_dir}}``", "Returns the base directory for the step."
    "``{{dir}}``", "Returns the active directory of the script. For project-scope scripts, this is identical to ``base_dir``. For sample scope scripts, this will be a direcotry within ``base_dir`` for sample related files."

3. Outputs
+++++++++++++++

Will be replaced with the filename specified in the named output. *e.g.* ``{{o:fasta.nucl}}`` will be replced according to the specifications in the output block named ``fasta.nucl``.

Each output block must contain 2 fields: ``scope`` and ``string``. The string contains a string describing the file to be stored in the equivalent slot. In the example above, there must be a block called ``fasta.nucl`` in the ``output`` block which can be defined as shown in the example in section **Lines for parameter file** below.


Requires:
~~~~~~~~~~~~~
Customizable
    
Output:
~~~~~~~~~~~~~
Customizable


Parameters that can be set        
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table::
    :header: "Parameter", "Values", "Comments"
    :widths: 5,10,10

    "output", "", "A block including 'scope' and 'string' definining the script outputs"


Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Demonstration of embedding various files and titles in a script file::

    pipe_gen_3:
        module:             pipe_generic
        base:               pipe_gen_2
        script_path: |
            project:                    {{project}}
            project:fasta.nucl:         {{project:fasta.nucl}}
            project:fasta.nucl:merge1   {{project:fasta.nucl::merge1}}

            sample:             {{sample::,}}
            sample:fastq.F:     {{sample:fastq.F}}
            sample:fastq.F:merge1     {{sample:fastq.F::merge1}}

            output:fasta.nucl:    {{o:fasta.nucl}}
        output:
            fasta.nucl:
                scope:      project
                string:       "{{base_dir}}{{project}}_new_pipegen3.fasta"


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

        # Get all user defined variables in string
        variables = list(set(re.findall(pattern="\{\{(.*?)\}\}",
                                        string=self.params["script_path"])))
        # Find all variables in outputs:
        try:
            for outp in self.params["output"].keys():
                # Check each 'output' has a 'string' and a 'scope' defined
                try:
                    # Extract the string from the {} and append to results:
                    result = re.findall(pattern="\{\{(.*?)\}\}",
                                        string=self.params["output"][outp]["string"])
                except KeyError:
                    raise AssertionExcept("Make sure you have a 'string' and 'scope' defined "
                                          "for output {output}!".format(output=outp))
                variables.extend(result)
        except KeyError:
            self.write_warning("No 'output' section defined. Are you sure this is what you intended?")

        # print list(set(variables))
        # Default scope is project
        if "scope" not in self.params:
            self.params["scope"] = "project"

        # Check the definition of all variables
        for variable in variables:
            var_def = re.findall(pattern="([^\:]*)\:?", string=variable)
            # print variable
            # print var_def
            # If variable scope is sample and the separator field (3rd slot) is not defined, change scope to sample
            if var_def[0]=="sample" and (len(var_def)<3 or not var_def[2]):
                self.params["scope"] = "sample"

        # print "scope:\n" + self.params["scope"]

        # # Set scope. If any of the variables begins with 'sample:' or is only 'sample', set scope to sample
        # # Otherwise, scope is project:
        # if any([re.match(string=match,pattern="\{\{sample[\:\}]") for match in variables]):
        #     self.params["scope"] = "sample"
        # else:
        #     self.params["scope"] = "project"

        # sys.exit()

    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        
        # TODO: Check that files exist for the analysis
        # print [self.params["output"][outp]["string"]
        #        for outp
        #        in self.params["output"]
        #        if self.params["output"][outp]["scope"]=="project"]
        # for outp in self.params["output"]:
        #     if self.params["output"][outp]["scope"]=="project"]:
        #         self.params["output"][outp]["string"] = self.format_script_path()
        # sys.exit()
        pass

    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        pass

#         variables = list(set(re.findall(pattern="(\{\{.*?\}\})", string=rawstring)))
#
#         self.script = ""
#         for variable in variables:
#             # Splitting by ':'
#             var_def = re.findall(pattern="([^\:]*)\:?",string=variable.strip('{{').strip('}}'))
#
#             if len(var_def) < 5:
#                 var_def = var_def + [''] * (4 - len(var_def))
#
#             if var_def[4] == "del":
#                 if var_def[0]=="project":
#                     file2del = self.sample
#                 self.script += """
# rm -rf {file}
#                 """

    def build_scripts(self):
        """ This is the actual script building function
            
        """
        if self.params["scope"] == "project":
            sample_list = ["project_data"]
        elif self.params["scope"] == "sample":
            sample_list = self.sample_data["samples"]
        else:
            raise AssertionExcept("'scope' must be either 'sample' or 'project'")

        for sample in sample_list:  # Getting list of samples out of samples_hash

            # Name of specific script:
            self.spec_script_name = self.set_spec_script_name(sample)
            self.script = ""

            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)

            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)

            # Creating new copy of params output so that samples don't modify the global version
            if "output" in self.params:
                self.params_output = deepcopy(self.params["output"])
            else:
                self.params_output = {}
            for outp in self.params_output:
                self.params_output[outp]["string"] = self.format_script_path(string=self.params_output[outp]["string"],
                                                                           use_dir=use_dir,
                                                                           sample=sample)

            self.script = self.format_script_path(string=self.params["script_path"],
                                                  use_dir=use_dir,
                                                  sample=sample)
            # # Try using function to include export (setenv) etc...

            # if "output" in self.params:
            for outp in self.params_output:
                # If script and output scopes are identical:
                if self.params_output[outp]["scope"] == self.params["scope"]:
                    # Store type after formatting:
                    self.sample_data[sample][outp] = \
                        self.format_script_path(string=self.params_output[outp]["string"],
                                                use_dir=use_dir,
                                                sample=sample)
                    self.stamp_file(self.sample_data[sample][outp])
                # If script is project and output is sample:
                elif self.params_output[outp]["scope"] == "sample":
                    self.write_warning("Writing sample scope output for project scope script!")
                    for outp_sample in self.sample_data["samples"]:
                        self.sample_data[outp_sample][outp] = \
                            self.format_script_path(string=self.params_output[outp]["string"],
                                                    use_dir=use_dir,
                                                    sample=outp_sample)
                        self.stamp_file(self.sample_data[outp_sample][outp])
                # If script is sample and output is project:
                elif self.params_output[outp]["scope"] == "project":
                    self.write_warning("Writing project scope output for sample scope script!")
                    self.sample_data["project_data"][outp] = \
                        self.format_script_path(string=self.params_output[outp]["string"],
                                                use_dir=use_dir,
                                                sample=sample)
                    self.stamp_file(self.sample_data["project_data"][outp])
                else:
                    pass
            # Wrapping up function. Leave these lines at the end of every iteration:
            self.local_finish(use_dir, sample_dir)
            self.create_low_level_script()

    def format_script_path(self, string, use_dir=None, sample=None):
        """

        :return:
        """

        if sample == "project_data":
            sample = None
        # Use self.get_base_sample_data() to get historic files (4th entry in input strings)

        rawstring=string
        # Try using function to include export (setenv) etc...

        variables = list(set(re.findall(pattern="(\{\{.*?\}\})", string=rawstring)))

        for variable in variables:
            # Splitting by ':'
            var_def = re.findall(pattern="([^\:]*)\:?",string=variable.strip('{{').strip('}}'))
            # var_def = variable.split(":")

            if len(var_def) < 4:
                var_def = var_def + [''] * (4 - len(var_def))
            if var_def[2] == "colon":
                var_def[2] = ":"
            # print variable
            # print var_def
            # print "---------------------"

            # ------------------------------
            if var_def[0] == "dir":
                rawstring = rawstring.replace(variable, use_dir)
                continue
            # ------------------------------
            if var_def[0] == "base_dir":
                rawstring = rawstring.replace(variable, self.base_dir)
                continue
            # ------------------------------
            if var_def[0] == "project":
                repl_str=""
                if not var_def[1]:  # Type not defined, use title
                    repl_str = self.sample_data["Title"]
                else:                # Type defined
                    if not var_def[3]:  # Base not defined. Use current
                        try:
                            repl_str = ("{!r}".format(self.sample_data["project_data"][var_def[1]])).strip("'")
                        except KeyError:
                            raise AssertionExcept(
                                "File type '{type}' not found in project scope".format(type=var_def[1]))
                    else:               # Base defined. Use defined base
                        # print self.get_base_sample_data().keys()
                        if var_def[3] not in self.get_base_sample_data():
                            raise AssertionExcept("No base '{base}' defined!".format(base=var_def[3]))
                        try:
                            repl_str = ("{!r}".format(self.get_base_sample_data()[var_def[3]]["project_data"][var_def[1]])).strip("'")
                        except KeyError:
                            raise AssertionExcept("No file of type '{type}' in project scope for base '{base}'".
                                                  format(type=var_def[1],
                                                         base=var_def[3]))

                rawstring = re.sub(pattern=re.escape(variable),
                                   repl=repl_str,
                                   string=rawstring)

                continue
            # ------------------------------
            if var_def[0] == "sample":
                # Create local copy of sample_data. If base is defined, this will be the base sample_data
                if not var_def[3]:  # Base not defined. Use current
                    sample_data = self.sample_data
                    # print "self"
                else:  # Base defined. Use defined base
                    if var_def[3] not in self.get_base_sample_data():
                        raise AssertionExcept("No base '{base}' defined!".format(base=var_def[3]))
                    sample_data = self.get_base_sample_data()[var_def[3]]
                    # print "base"

                if var_def[2]:  # Separator is defined
                    if not var_def[1]:  # Type is not defined
                        repl_str = var_def[2].join(sample_data["samples"])
                    else:           # Type is defined
                        try:
                            repl_str = var_def[2].join([("{!r}".format(sample_data[sample][var_def[1]])).strip("'")
                                                  for sample
                                                  in sample_data["samples"]])

                        except KeyError:
                            raise AssertionExcept("File type '{type}' not found in all samples".format(type=var_def[1]))
                else:           # Separator is not defined
                    if not sample:
                        raise AssertionExcept("Trying to parse sample in project scope script!")
                    if not var_def[1]:
                        repl_str = sample
                    else:
                        try:
                            repl_str=("{!r}".format(sample_data[sample][var_def[1]])).strip("'")
                        except KeyError:
                            raise AssertionExcept("File type '{type}' not found in sample".format(type=var_def[1]),
                                                  sample)

                rawstring = re.sub(pattern=re.escape(variable),
                                   repl=repl_str,
                                   string=rawstring)

                continue
            # ------------------------------
            if var_def[0] == "o":
                try:
                    rawstring = re.sub(pattern=re.escape(variable),
                                         repl=("{!r}".format(self.params_output[var_def[1]]["string"])).strip("'"),
                                         string=rawstring)
                except KeyError:
                    raise AssertionExcept("Error embedding output '{var}'".format(var=variable), sample)
                continue
            # ------------------------------
            #  variable does not match any of the expected formats:
            raise AssertionExcept('Variable {var} in script_path not identified'.format(var=variable))
        # print "Return: ", rawstring
        return rawstring
