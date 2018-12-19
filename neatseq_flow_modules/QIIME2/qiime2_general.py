# -*- coding: UTF-8 -*-
""" 
``qiime2_general``
-----------------------------------------------------------------
:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

.. Attention:: This module is in Beta version. It is not issue-free and will be improved periodically.

A module to include any QIIME2 plugin method, pipeline or visualiation.

The required plugin and method are specified in the ``script_path`` line, as they would appear in the command line, *e.g.*::

    script_path:   qiime dada2 denoise-paired

The module will identify the required inputs for the method and extract them from the appropriate slots. If they are not found, an exception will be thrown.

If more than one type is legitimate fdr a method, and both exist in the project, NeatSeq-Flow will complain. You can
either remove the extra type with ``manage_types`` module or specify the type to use with the ``type`` parameter.

Plugins which require metadata files, passed as argument ``--m-metadata-file``, will look for a file in slot ``metadata``.
In order to specify a metadata file in the parameter file, pass the ``--m-metadata-file`` in the ``redirects`` section.

All ``redirects`` argument values are searched for in the "project_data". Thus, you can specify slots to use for
redirected arguments.

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Plugin- and method-specific.

    
Output:
~~~~~~~~~~~~~

Plugin- and method-specific.

                
Parameters that can be set        
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"

    "store_output", "list of output parameters", "These parameters will be stored as file types for use by downstream modules"
    "export", "empty or list of output parameters", "If empty, all outputs will be exported, *i.e.* unzipped with qiime tools export. If list of parameters, only those types will be exported."

    
Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DADA2 plugin, with export of stats output::

    dada2:                      # Name of this step
        module:                     qiime2_general
        base:                       import
        script_path:                qiime dada2 denoise-single #paired
        export:
            - --o-denoising-stats
        redirects:
            --p-trim-left:         10
            --p-trunc-len:         100

Classical visualization. Only ``base`` and ``script_path``::

    dada2_vis_summary:                      # Name of this step
        module:                     qiime2_general
        base:                       dada2
        script_path:                qiime feature-table summarize

Store only particular outputs in type index::

    diversity:                      # Name of this step
        module:                     qiime2_general
        base:                       phylogeny
        script_path:                qiime diversity core-metrics-phylogenetic
        export:                     --o-rarefied-table
        store_output:
            - --o-rarefied-table
            - --o-faith-pd-vector
            - --o-weighted-unifrac-distance-matrix
            - --o-weighted-unifrac-pcoa-results
            - --o-weighted-unifrac-emperor
        redirects:
            --p-sampling-depth:     50000

    taxonomy_tabulate:                      # Name of this step
        module:                     qiime2_general
        base:                       taxonomy
        script_path:                qiime metadata tabulate
        redirects:
            --m-input-file:         "{{FeatureData[Taxonomy]}}"

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
        with open(os.path.join(os.path.dirname(os.path.realpath(__file__)),"qiime2_arguments_index.yml"),"r") as fileh:
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

        # for redir,redir_val in self.params["redir_params"].iteritems():
        #     if isinstance(redir_val, str) and redir_val in self.sample_data["project_data"]:
        #         # redir_val can also be a list, in which case it can't be project_data slots.
        #         self.params["redir_params"][redir] = self.sample_data["project_data"][redir_val]

        # The following enables passing types through redirects:
        # Is done like this (see also code at end of loop), to make it easier to perform on sample scope in the future
        # old_redirs = dict()
        for redir in self.params["redir_params"]:
            if isinstance(self.params["redir_params"][redir],str):
                if re.match("^\{\{.*\}\}$",self.params["redir_params"][redir]):
                    # old_redirs[redir] = self.params["redir_params"][redir]
                    type2use = list(set(re.findall(pattern="\{\{(.*?)\}\}", string=self.params["redir_params"][redir])))
                    if len(type2use)>1:
                        raise AssertionExcept("You are trying to pass two types through redirect '{redir}'".format(redir=redir))
                    type2use = type2use[0]

                    try:
                        self.params["redir_params"][redir] = self.sample_data["project_data"][type2use]
                    except KeyError:
                        raise AssertionExcept("Type '{type}' does not exist!".format(type=type2use))

        # If a metadata file is required, is not explicit in redirects and 'metadata' slot exists, use it.
        for param_type in ["required","optional"]:
            if param_type in self.method_index and "--m-metadata-file" in self.method_index[param_type]:
                if "--m-metadata-file" not in self.params["redir_params"] \
                    and "metadata" in self.sample_data["project_data"]:
                    self.params["redir_params"]["--m-metadata-file"] = self.sample_data["project_data"]["metadata"]

        # Convert required into list and making sure all required exist
        # Done here so that -m-- files can be included automatically in redirects
        if "required" in self.method_index:
            if isinstance(self.method_index["required"], str):
                self.method_index["required"] = [self.method_index["required"]]

            if not all(map(lambda x: x in self.params["redir_params"], self.method_index["required"])):
                raise AssertionExcept("The following parameters are required for method {method} of plugin {plugin}: "
                                      "{required}".format(plugin=self.plugin,
                                                          method=self.method,
                                                          required=", ".join(self.method_index["required"])))

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
        if inputs:
            self.script += inputs + " \\\n\t"
        if outputs:
            self.script += outputs


        if "store_output" in self.params:
            self.script += " \\\n\t--output-dir {dir}{outdir}".format(dir=use_dir,outdir="more_results")

        self.script += "\n\n"

        for outp, typ in output_index.iteritems():
            self.sample_data[sample][typ] = "{dir}{title}.{method}.{outp}.{ext}". \
                                                                            format(dir=sample_dir,
                                                                                   title=self.sample_data["Title"],
                                                                                   method=self.method,
                                                                                   ext="qzv" if typ == "Visualization" else "qza",
                                                                                   outp=edit_qiime_params(outp))
            self.stamp_file(self.sample_data[sample][typ])
            if "export" in self.params: # and not outp == "--o-visualization":
                if not self.params["export"] or outp in self.params["export"]:
                    # Adding export code if requested
                    self.script += """
qiime tools export \\
    --input-path {inppath} \\
    --output-path {outpath}

""".format(inppath=self.sample_data[sample][typ],
           outpath="{dir}{title}.{method}.{outp}".format(dir=sample_dir,
                                                         title=self.sample_data["Title"],
                                                         method=self.method,
                                                         outp=edit_qiime_params(outp)))


        self.local_finish(use_dir, sample_dir)
        self.create_low_level_script()

        # for redir in old_redirs:
        #     self.params["redir_params"][redir] = old_redirs[redir]