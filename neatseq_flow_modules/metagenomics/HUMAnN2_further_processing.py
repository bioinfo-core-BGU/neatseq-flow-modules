# -*- coding: UTF-8 -*-
""" 
``HUMAnN2_further_processing``
--------------------------------

:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A module for running ``HUMAnN2`` utilities ``humann2_regroup_table`` and ``humann2_rename_table``.


Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* HUMAnN2 results:

    * ``sample_data[<sample>]["HUMAnN2.genefamilies"]``

    

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Puts the ``humann2_regroup_table`` output files in:

    * ``self.sample_data[sample]["HUMAnN2.genefamilies.uniref50_eggnog"]``
    * ``self.sample_data[sample]["HUMAnN2.genefamilies.uniref50_go"]``
    * ``self.sample_data[sample]["HUMAnN2.genefamilies.uniref50_infogo1000"]``
    * ``self.sample_data[sample]["HUMAnN2.genefamilies.uniref50_ko"]``
    * ``self.sample_data[sample]["HUMAnN2.genefamilies.uniref50_level4ec"]``
    * ``self.sample_data[sample]["HUMAnN2.genefamilies.uniref50_pfam"]``
    * ``self.sample_data[sample]["HUMAnN2.genefamilies.uniref50_rxn"]``

* Puts the ``humann2_rename_table`` output files in:

    * ``self.sample_data[sample]["HUMAnN2.genefamilies.uniref50_eggnog_names"]``
    * ``self.sample_data[sample]["HUMAnN2.genefamilies.uniref50_go_names"]``
    * ``self.sample_data[sample]["HUMAnN2.genefamilies.uniref50_infogo1000_names"]``
    * ``self.sample_data[sample]["HUMAnN2.genefamilies.uniref50_ko_names"]``
    * ``self.sample_data[sample]["HUMAnN2.genefamilies.uniref50_level4ec_names"]``
    * ``self.sample_data[sample]["HUMAnN2.genefamilies.uniref50_pfam_names"]``
    * ``self.sample_data[sample]["HUMAnN2.genefamilies.uniref50_rxn_names"]``


Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "humann2_regroup_table", "", "Block containing ``path`` to ``humann2_regroup_table``, and a ``redirects`` block if necessary."
    "humann2_rename_table", "", "Block containing ``path`` to ``humann2_rename_table``, and a ``redirects`` block if necessary."
    "scope", "sample|project", "Decides whether the scripts should operate on sample or project-wise reports"

    
Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For regrouping by all grouping options::

    HUMAnN2_uniref50_further_analysis:
        module: HUMAnN2_further_processing
        base: HUMAnN2_uniref50_hardtrimmed_reads
        script_path:
        scope: project
        qsub_params:
            -pe: shared 30
        humann2_regroup_table:
            path: humann2_regroup_table
        humann2_rename_table:
            path: humann2_rename_table

For regrouping by a specific grouping option::

    HUMAnN2_uniref50_further_analysis:
        module: HUMAnN2_further_processing
        base: HUMAnN2_uniref50_hardtrimmed_reads
        script_path:
        scope: project
        qsub_params:
            -pe: shared 30
        humann2_regroup_table:
            path: humann2_regroup_table
            mpa_dir:                     # PATH to the group mapping files directory, if not given will use the "mpa_dir" Environment Variable 
            group_mapping_files_pattern: map_{{grouping_type}}_{{protein_database}}.txt.gz # This string is used as a patten for the group mapping files
                                                                                       # The {{grouping_type}} will be replaced by the grouping type (eggnog\\go\\infogo1000\\ko\\pfam)
                                                                                       # The {{protein_database}} will be replaced by the protein database (uniref50\\uniref90)
            redirects:
                --groups:   go      # Note: NOT `uniref90_go`. The protein database name is taken from the base HUMAnN2 instance!
                                    # The possible options are eggnog\\go\\infogo1000\\ko\\pfam
        humann2_rename_table:
            path: humann2_rename_table
            group_renaming_files_pattern: map_{{grouping_type}}_name.txt.gz # This string is used as a patten for the group renaming files
                                                                                       # The {{grouping_type}} will be replaced by the grouping type (eggnog\\go\\infogo1000\\ko\\pfam)

References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`HUMAnN2 home page <http://huttenhower.sph.harvard.edu/humann2>`_


"""

import os
import sys
import re
from neatseq_flow.PLC_step import Step,AssertionExcept

import yaml
try:
    from yaml import CSafeLoader as Loader
except ImportError:
    from yaml import SafeLoader as Loader


__author__ = "Menachem Sklarz"
__version__ = "1.6.0"


class Step_HUMAnN2_further_processing(Step):

    
    def step_specific_init(self):
        """ Called on intiation
            Good place for parameter testing.
            Wrong place for sample data testing
        """
        self.shell = "bash"      # Can be set to "bash" by inheriting instances

        # Setting default scope to sample
        if "scope" not in self.params:
            self.params["scope"] = "sample"


        # Strcuture:
        # Key: grouping type (--groups in humann2_regroup_table)
        # Value[0]: Filename suffix
        # Value[1]: Equivalent value of --names in humann2_rename_table
        self.type_index = {
            "eggnog":     ["EggNog",        "eggnog"],
            "go":         ["GO",            "go"],
            "infogo1000": ["GOinfo1000",    "infogo1000"],
            "ko":         ["KO",            "kegg-orthology"],
            "level4ec":   ["level4EC",      "ec"],
            "pfam":       ["PFAM",          "pfam"],
            "rxn":        ["RXN",           "metacyc-rxn"]
        }
        
        # self.type_index = yaml.load(type_index,  Loader=yaml.Loader)
        
        self.group_mapping_files_pattern  = 'map_{grouping_type}_{protein_database}.txt.gz'
        self.group_renaming_files_pattern = 'map_{grouping_type}_name.txt.gz'
        
        grouptype = None

        if "humann2_regroup_table" not in self.params:
            raise AssertionExcept("Please include a 'humann2_regroup_table' dict in params. See docs")

        if not isinstance(self.params["humann2_regroup_table"], dict):
            raise AssertionExcept("Please pass at least path in 'humann2_regroup_table' dict. See docs.")
        if "path" not in self.params["humann2_regroup_table"]:
            raise AssertionExcept("Please pass a path to humann2_regroup_table. See docs")

        if "redirects" in self.params["humann2_regroup_table"]:

            if "-g" in self.params["humann2_regroup_table"]["redirects"]:
                grouptype = self.params["humann2_regroup_table"]["redirects"]["-g"]
            elif "--groups" in self.params["humann2_regroup_table"]["redirects"]:
                grouptype = self.params["humann2_regroup_table"]["redirects"]["--groups"]
            elif "--custom" in self.params["humann2_regroup_table"]["redirects"]:
                grouptype = self.params["humann2_regroup_table"]["redirects"]["--custom"]
            elif "-c" in self.params["humann2_regroup_table"]["redirects"]:
                grouptype = self.params["humann2_regroup_table"]["redirects"]["-c"]
            else:
                pass
                # raise AssertionExcept("Please pass -g or --groups via humann2_regroup_table redirects")
        else:
            self.params["humann2_regroup_table"]["redirects"] = dict()
            # pass
            # raise AssertionExcept("You must pass a redirects section to humann2_regroup_table with at least --groups")

        if not grouptype:
            self.grouptype = self.type_index.keys()
        elif ("--custom" or "-c") not in self.params["humann2_regroup_table"]["redirects"]:
            if grouptype not in self.type_index:
                raise AssertionExcept("Unrecognised --groups value")
            else:
                self.grouptype = [grouptype]
        else:
            self.grouptype = [grouptype]

        if "humann2_rename_table" in self.params:
            if not isinstance(self.params["humann2_rename_table"], dict):
                raise AssertionExcept("Please pass at least path in 'humann2_rename_table' dict. See docs.")
            if "redirects" in self.params["humann2_rename_table"]:
                if ({"--reverse"} & set(self.params["humann2_rename_table"]["redirects"])):
                    raise AssertionExcept("--reverse is not supported at the moment.")
                if ({"-n","--names"} & set(self.params["humann2_rename_table"]["redirects"])):
                    raise AssertionExcept("--names are determined automatically in humann2_rename_table.")
            if "path" not in self.params["humann2_rename_table"]:
                raise AssertionExcept("Please pass a path to humann2_rename_table. See docs")


    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        pass

        if "HUMAnN2.prot.db" not in self.sample_data["project_data"]:
            raise AssertionExcept("To use further processing of HUMAnN2 results, you need to have defined "
                                  "'protein-database' in the HUMAnN2 instance!")

        if self.params["scope"] == "project":
            sample_list = ["project_data"]
        elif self.params["scope"] == "sample":
            sample_list = self.sample_data["samples"]
        else:
            raise AssertionExcept("'scope' must be either 'sample' or 'project'")

        for sample in sample_list:  # Getting list of samples out of samples_hash
            if "HUMAnN2.genefamilies" not in self.sample_data[sample]:
                raise AssertionExcept("Please include a HUMAnN2 instance before this instance!", sample)


    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        pass

    def build_scripts(self):
        """ This is the actual script building function
            Most, if not all, editing should be done here 
            HOWEVER, DON'T FORGET TO CHANGE THE CLASS NAME AND THE FILENAME!
        """

        if self.params["scope"] == "project":
            sample_list = ["project_data"]
        elif self.params["scope"] == "sample":
            sample_list = self.sample_data["samples"]
        else:
            raise AssertionExcept("'scope' must be either 'sample' or 'project'")

        for sample in sample_list:  # Getting list of samples out of samples_hash

            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)



            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)

            if "redirects" in self.params["humann2_regroup_table"]:
                redirects = " \\\n\t".join(
                    [key + " " + (val if val else "")
                     for key, val
                     in self.params["humann2_regroup_table"]["redirects"].items()])
            else:
                redirects = ""

            humanntype = "genefamilies"
            
            if "mpa_dir" not in self.params.keys():
                self.params["mpa_dir"]='$mpa_dir'
            if "group_mapping_files_pattern" not in self.params["humann2_regroup_table"]:
                self.params["humann2_regroup_table"]["group_mapping_files_pattern"]=self.group_mapping_files_pattern
                
            for grouping_ty in self.grouptype:
                # Name of specific script:
                #full_grouping_ty = "_".join([self.sample_data["project_data"]["HUMAnN2.prot.db"],grouping_ty])
                if ("--custom" or "-c") not in self.params["humann2_regroup_table"]["redirects"]:
                    full_grouping_ty = self.params["humann2_regroup_table"]["group_mapping_files_pattern"].format(grouping_type=grouping_ty,
                                                                                                                protein_database=self.sample_data["project_data"]["HUMAnN2.prot.db"])
                else:
                    full_grouping_ty = grouping_ty
                    
                full_grouping_ty = os.path.join(self.params["mpa_dir"],full_grouping_ty)
                self.spec_script_name = self.set_spec_script_name(".".join([sample if not sample=="project_data"
                                                                                    else self.sample_data["Title"],
                                                                            humanntype,
                                                                            self.type_index[grouping_ty][0]]))
                self.script = ""
                inputfn = self.sample_data[sample]["HUMAnN2."+humanntype]
                outfn = ("."+self.type_index[grouping_ty][0]).join(os.path.splitext(os.path.basename(inputfn)))

                if grouping_ty!='rxn':
                    self.script += """
{setenv}{path} \\
\t--output {dir}{outfn} \\
\t--custom {group} \\
\t--input {inputfn} \\
\t{redirs}
 
""".format(setenv=self.get_setenv_part(),
           path=self.params["humann2_regroup_table"]["path"],
           dir=use_dir,
           outfn=outfn,
           inputfn=inputfn,
           group=full_grouping_ty,
           redirs=redirects)
           
                else:
                    self.script += """
{setenv}{path} \\
\t--output {dir}{outfn} \\
\t-g {group} \\
\t--input {inputfn} \\
\t{redirs}
 
""".format(setenv=self.get_setenv_part(),
           path=self.params["humann2_regroup_table"]["path"],
           dir=use_dir,
           outfn=outfn,
           inputfn=inputfn,
           group=full_grouping_ty,
           redirs=redirects)
           

                self.sample_data[sample][".".join(["HUMAnN2",humanntype,full_grouping_ty])] = (sample_dir + outfn)
                self.stamp_file(self.sample_data[sample][".".join(["HUMAnN2",humanntype,full_grouping_ty])])




#
                if "humann2_rename_table" in self.params:
                    # Adding code for normalization if required
                    # Checking redirects exists although in init this was checked already
                    if "redirects" in self.params["humann2_rename_table"]:
                        redirects = " \\\n\t".join(
                            [key + " " + (val if val else "")
                             for key, val
                             in self.params["humann2_rename_table"]["redirects"].items()])
                    else:
                        redirects = ""
                        
                    if "group_renaming_files_pattern" not in self.params["humann2_rename_table"]:
                        self.params["humann2_rename_table"]["group_renaming_files_pattern"]=self.group_renaming_files_pattern
                        
                    
                    
                                                                                                                
                    
                    self.script += """\
# Adding code for normalizing genefamilies and pathabundance tables\n\n
{path} \\
\t-i {infn} \\
\t-n {names} \\
\t--custom {group} \\
\t-o {outfn} \\
\t{redirs}

    """.format(path=self.params["humann2_rename_table"]["path"],
               infn   = use_dir+outfn,
               names=self.type_index[grouping_ty][1],
               group=full_grouping_ty,
               # names  = os.path.join(self.params["mpa_dir"],
                                     # self.params["humann2_rename_table"]["group_renaming_files_pattern"].format(grouping_type=self.type_index[grouping_ty][1])),
               outfn  = use_dir + "_names".join(os.path.splitext(outfn)),
               redirs = redirects)

                    self.sample_data[sample][".".join(["HUMAnN2",humanntype,full_grouping_ty,"names"])] = sample_dir + "_names".join(os.path.splitext(outfn))

                    self.stamp_file(self.sample_data[sample][".".join(["HUMAnN2",humanntype,full_grouping_ty,"names"])])


                # Move all files from temporary local dir to permanent base_dir
                self.local_finish(use_dir, self.base_dir)
                self.create_low_level_script()


    def make_sample_file_index(self):
        """ Make file containing samples and target file names for use by kraken analysis R script
        """
        pass