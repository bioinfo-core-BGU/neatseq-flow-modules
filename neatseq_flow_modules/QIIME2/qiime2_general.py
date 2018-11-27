# -*- coding: UTF-8 -*-
""" 
``qiime2_general``
-----------------------------------------------------------------
:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.



Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    


    
Output:
~~~~~~~~~~~~~



                
Parameters that can be set        
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"



    
Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::



References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"""



import os
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


class Step_qiime2_general(Step):
    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances

        # Read YAML of plugin arguments
        with open(os.path.join(os.path.dirname(os.path.realpath(__file__)),"qiime2_arguments_index.yaml"),"r") as fileh:
            filelines = fileh.readlines()

        self.qiime_args = yaml.load("".join(filelines),  Loader=yaml.Loader)

        # extract qiime path, plugin name and method/pipeline/visualization from script_path
        self.qiime_path = self.params["script_path"].split(" ")[0]
        self.plugin = self.params["script_path"].split(" ")[1]
        self.method = self.params["script_path"].split(" ")[2]

        # Check plugin and method are recognized
        if self.plugin not in self.qiime_args:
            raise AssertionExcept("Plugin '{plugin}' is not one of: {plugins}!".
                                  format(plugin=self.plugin,
                                         plugins=", ".join(self.qiime_args.keys())))
        if self.method not in self.qiime_args[self.plugin]:
            raise AssertionExcept("Plugin '{method}' is not one of: {methods}!".
                                  format(method=self.method,
                                         methods=", ".join(self.qiime_args[self.plugin].keys())))
        # Get argument index for method
        self.method_index = self.qiime_args[self.plugin][self.method]

        # Convert required into list
        if "required" in self.method_index:
            if isinstance(self.method_index["required"], str):
                self.method_index["required"] = [self.method_index["required"]]

            if not all(map(lambda x: x in self.params["redir_params"], self.method_index["required"])):
                raise AssertionExcept("The following parameters are required for method {method} of plugin {plugin}: "
                                      "{required}".format(plugin=self.plugin,
                                                          method=self.method,
                                                          required=", ".join(self.method_index["required"])))

        # print self.method_index

    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
            Here you should do testing for dependency output. These will NOT exist at initiation of this instance. They are set only following sample_data updating
        """

        self.input_dict = dict()

        # Check inputs exist in sample_data
        # Merge inputs and optional_inputs dicts:
        all_inputs = dict()
        if "inputs" in self.method_index:
            all_inputs.update(self.method_index["inputs"])
        if "optional_inputs" in self.method_index:
            all_inputs.update(self.method_index["optional_inputs"])

        # for inpflag,inptype in self.method_index["inputs"].iteritems():

        for inpflag,inptype in all_inputs.iteritems():
            # Convert single value strings into lists:
            if isinstance(inptype,str):
                inptype = [inptype]

            # Make sure at least one of the list members exists in sample_data
            # (only for required inputs!)
            if inpflag in self.method_index["inputs"]:
                if not any([itype
                            in self.sample_data["project_data"]
                            for itype
                            in inptype]):
                    raise AssertionExcept("For {meth}, you need one of '{type}' defined.".format(meth=self.method,
                                                                                             type=", ".join(inptype)))
            # inptype = list(set(inptype) or set(self.sample_data["project_data"].keys()))[0]
            # Get list of types that exist in sample_data:
            inptype = list(set(inptype) & set(self.sample_data["project_data"].keys()))
            # If type in params, try using type or fail
            if "type" in self.params:
                if not isinstance(self.params["type"],str):
                    raise AssertionExcept("'type' parameter must be a single value.")
                if self.params["type"] in inptype:
                    inptype = self.params["type"]
                else:
                    raise AssertionExcept("Type '{type}' not defined!".format(type=self.params["type"]))
            else:     # If not, and there are more than one available type, raise except
                if len(inptype)>1:
                    raise AssertionExcept("The following types are available for {flag}: {types}. Please specify which "
                                          "one to use with a 'type' parameter".format(flag=inpflag,
                                                                                      types=", ".join(inptype)))
                elif len(inptype) == 0:  # Type does not exist in project_data:
                    continue
                else:
                    # Convert single-value list inptype to string
                    inptype = inptype[0]
            self.input_dict[inpflag] = inptype

        # print "------------ %s ---------" % self.get_step_name()
        # print self.input_dict
        # sys.exit()

    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        
        pass

    def build_scripts(self):

        pass

        # Name of specific script:
        self.spec_script_name = self.set_spec_script_name()
        self.script = ""

        # Defining sample as "project_data".
        # This will enable possibly defining this on sample scope as well, in the far future
        sample = "project_data"

        # Make a dir for the current sample:
        sample_dir = self.make_folder_for_sample(sample=sample)
        # This line should be left before every new script. It sees to local issues.
        # Use the dir it returns as the base_dir for this step.
        use_dir = self.local_start(sample_dir)

        # The following enables passing types through redirects:
        # Is done like this (see also code at end of loop), to make it easier to perform on sample scope in the future
        old_redirs = dict()
        for redir in self.params["redir_params"]:
            if isinstance(self.params["redir_params"][redir],str):
                if re.match("^\{\{.*\}\}$",self.params["redir_params"][redir]):
                    old_redirs[redir] = self.params["redir_params"][redir]
                    type2use = list(set(re.findall(pattern="\{\{(.*?)\}\}", string=self.params["redir_params"][redir])))
                    if len(type2use)>1:
                        raise AssertionExcept("You are trying to pass two types through redirect '{redir}'".format(redir=redir))
                    type2use = type2use[0]

                    try:
                        self.params["redir_params"][redir] = self.sample_data[sample][type2use]
                    except KeyError:
                        raise AssertionExcept("Type '{type}' does not exist!".format(type=type2use))


        # Build inputs part:
        inputs = " \\\n\t".join([" ".join([inp, self.sample_data[sample][typ]])
                                 for inp, typ
                                 in self.input_dict.iteritems()])

        # Build inputs part:
        if "store_output" in self.params:
            output_index = {outp:typ
                            for outp, typ
                            in self.method_index["outputs"].iteritems()
                            if outp in self.params["store_output"]}
            # output_index["--output-dir"] = "results_dir"
        else:
            output_index = self.method_index["outputs"]
        # print output_index
        # sys.exit()
        outputs = " \\\n\t".join(
            [" ".join([outp,
                       "{dir}{title}.{method}.{outp}.{ext}".format(dir=use_dir,
                                                                 title=self.sample_data["Title"],
                                                                 method=self.method,
                                                                 ext="qzv" if typ=="Visualization" else "qza",
                                                                 outp=edit_qiime_params(outp))])
             for outp, typ
             # in self.method_index["outputs"].iteritems()])
             in output_index.iteritems()])

        # Create script
        self.script = """\
# Delete existing results dir
if [ -e {dir}{outdir} ]; then rm -rf {dir}{outdir}; fi

""".format(dir=use_dir,
           outdir="more_results")


        self.script += self.get_script_const()
        self.script += inputs + " \\\n\t"
        self.script += outputs
        if "store_output" in self.params:
            self.script += " \\\n\t--output-dir {dir}{outdir}".format(dir=use_dir,outdir="more_results")

        for outp, typ in output_index.iteritems():
            self.sample_data[sample][typ] = "{dir}{title}.{method}.{outp}.{ext}". \
                                                                            format(dir=sample_dir,
                                                                                   title=self.sample_data["Title"],
                                                                                   method=self.method,
                                                                                   ext="qzv" if typ == "Visualization" else "qza",
                                                                                   outp=edit_qiime_params(outp))
            self.stamp_file(self.sample_data[sample][typ])

        self.local_finish(use_dir, sample_dir)
        self.create_low_level_script()

        for redir in old_redirs:
            self.params["redir_params"][redir] = old_redirs[redir]