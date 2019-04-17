# -*- coding: UTF-8 -*-
""" 
``HUMAnN2_further_processing``
--------------------------------

:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

.. Note:: This module was developed as part of a study led by Dr. Jacob Moran Gilad

A module for running ``HUMAnN2``:


Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* fastq files, either forward or single:

    * ``sample_data[<sample>]["fastq.F"]``
    * ``sample_data[<sample>]["fastq.S"]``
    
    

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


self.sample_data[sample]["HUMAnN2.genefamilies.norm"]
self.sample_data[sample]["HUMAnN2.pathabundance.norm"]


* Puts the ``HUMAnN2`` output files in:  

    * ``self.sample_data[sample]["HUMAnN2.genefamilies"]``
    * ``self.sample_data[sample]["HUMAnN2.pathabundance"]``
    * ``self.sample_data[sample]["HUMAnN2.pathcoverage"]``

* If "renorm_table" is set in params:

    * ``self.sample_data[sample]["HUMAnN2.genefamilies.norm"]``
    * ``self.sample_data[sample]["HUMAnN2.pathabundance.norm"]``
    
* If "join_tables" is set in params:
    
    * ``self.sample_data["project_data"]["HUMAnN2.genefamilies"]``
    * ``self.sample_data["project_data"]["HUMAnN2.pathabundance"]``
    * ``self.sample_data["project_data"]["HUMAnN2.pathcoverage"]``

* If "join_tables" and "renorm_table" are set in params:
    
    * ``self.sample_data["project_data"]["HUMAnN2.genefamilies.norm"]``
    * ``self.sample_data["project_data"]["HUMAnN2.pathabundance.norm"]``
    * ``self.sample_data["project_data"]["HUMAnN2.pathcoverage"]``
    
Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "renorm_table",      "Empty or parameters to pass to ``humann2_renorm_table``", "Runs ``humann2_renorm_table`` on HUMAnN2 outputs with the specified parameters"
    "join_tables", "", "Runs ``humann2_join_tables`` to gather all sample tables."
    "humann2_join_tables_path", "", "Path to ``humann2_join_tables``. If not passed, will try guessing"
    "humann2_renorm_table_path", "", "Path to ``humann2_renorm_table``. If not passed, will try guessing"

    
Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::


    HUMAnN2_1:
        module: HUMAnN2
        base: trim1
        script_path: {Vars.paths.humann2}
        join_tables: 
        renorm_table: --units cpm -p
        redirects:
            --bowtie2: /path/to/bowtie2
            --gap-fill: on
            --input-format: fastq
            --metaphlan: {Vars.paths.metaphlan2}
            --minpath: on
            --nucleotide-database: {Vars.databases.chocophlan}
            --protein-database: {Vars.databases.uniref}
            --threads: 30
            
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

#         type_index = """
# uniref50_eggnog:    ["EggNog",  "eggnog"]
# uniref50_go:        ["GO",      "go"]
# uniref50_infogo1000: ["GOinfo",  "infogo1000"]
# uniref50_ko:        ["KO",      "kegg-orthology"]
# uniref50_level4ec:  ["level4EC","ec"]
# uniref50_pfam:      ["PFAM",    "pfam"]
# uniref50_rxn:       ["RXN",     "metacyc-rxn"]
# uniref90_eggnog:    ["EggNog",  "eggnog"]
# uniref90_go:        ["GO",      "go"]
# uniref90_infogo1000: ["GOinfo",  "infogo1000"]
# uniref90_ko:        ["KO",      "kegg-orthology"]
# uniref90_level4ec:  ["level4EC","ec"]
# uniref90_pfam:      ["PFAM",    "pfam"]
# uniref90_rxn:       ["RXN",     "metacyc-rxn"]
# """

        type_index = """
eggnog:    ["EggNog",  "eggnog"]
go:        ["GO",      "go"]
infogo1000: ["GOinfo",  "infogo1000"]
ko:        ["KO",      "kegg-orthology"]
level4ec:  ["level4EC","ec"]
pfam:      ["PFAM",    "pfam"]
rxn:       ["RXN",     "metacyc-rxn"]
"""


        self.type_index = yaml.load(type_index,  Loader=yaml.Loader)

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
            else:
                pass
                # raise AssertionExcept("Please pass -g or --groups via humann2_regroup_table redirects")
        else:
            self.params["humann2_regroup_table"]["redirects"] = dict()
            # pass
            # raise AssertionExcept("You must pass a redirects section to humann2_regroup_table with at least --groups")

        if not grouptype:
            self.grouptype = self.type_index.keys()
        elif grouptype not in self.type_index:
            raise AssertionExcept("Unrecognised --groups value")
        else:
            self.grouptype = [grouptype]

        if "humann2_rename_table" in self.params:
            if not isinstance(self.params["humann2_rename_table"], dict):
                raise AssertionExcept("Please pass at least path in 'humann2_rename_table' dict. See docs.")
            if "redirects" in self.params["humann2_rename_table"]:
                if ({"--custom","--reverse"} & set(self.params["humann2_rename_table"]["redirects"])):
                    raise AssertionExcept("--custom and --reverse are not supported at the moment.")
                if ({"-n","--names"} & set(self.params["humann2_rename_table"]["redirects"])):
                    raise AssertionExcept("--names are determined automatically in humann2_rename_table.")
            if "path" not in self.params["humann2_rename_table"]:
                raise AssertionExcept("Please pass a path to humann2_rename_table. See docs")


    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        pass

        if "HUMAnN2.prot.db" not in self.sample_data:
            raise AssertionExcept("To use further processing of HUMAnN2 results, you need to have defined "
                                  "'protein-database' in the HUMAnN2 instance!")


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
            for grouping_ty in self.grouptype:
                # Name of specific script:
                full_grouping_ty = "_".join([self.sample_data["HUMAnN2.prot.db"],grouping_ty])

                self.spec_script_name = self.set_spec_script_name(".".join([sample if not sample=="project_data"
                                                                                    else self.sample_data["Title"],
                                                                            humanntype,
                                                                            grouping_ty]))
                self.script = ""
                inputfn = self.sample_data[sample]["HUMAnN2."+humanntype]
                outfn = ("."+self.type_index[grouping_ty][0]).join(os.path.splitext(os.path.basename(inputfn)))


                self.script += """
{setenv}{path} \\
\t--output {dir}{outfn} \\
\t--groups {group} \\
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

                    self.script += """\
# Adding code for normalizing genefamilies and pathabundance tables\n\n
{path} \\
\t-i {infn} \\
\t-n {names} \\
\t-o {outfn} \\
\t{redirs}

    """.format(path=self.params["humann2_rename_table"]["path"],
               infn=use_dir+outfn,
               names=self.type_index[grouping_ty][1],
               outfn=use_dir + "_names".join(os.path.splitext(outfn)),
               redirs=redirects)

                    self.sample_data[sample][".".join(["HUMAnN2",humanntype,full_grouping_ty,"names"])] = sample_dir + "_names".join(os.path.splitext(outfn))

                    self.stamp_file(self.sample_data[sample][".".join(["HUMAnN2",humanntype,full_grouping_ty,"names"])])


                # Move all files from temporary local dir to permanent base_dir
                self.local_finish(use_dir, self.base_dir)
                self.create_low_level_script()


    def make_sample_file_index(self):
        """ Make file containing samples and target file names for use by kraken analysis R script
        """
        pass