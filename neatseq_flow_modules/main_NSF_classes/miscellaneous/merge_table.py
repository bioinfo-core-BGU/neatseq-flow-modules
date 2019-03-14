# -*- coding: UTF-8 -*-
"""
``merge_table``
------------------------


:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A module for merging sample tables into a single project-wide table, or into group tables by category.

The table can be with or without a header line.

Can be used for merging fasta and fastq files as well.

.. important:: When merging by category, the sample names will be set to the category level names for all subsequent steps.

.. Tip:: You can merge several types at once by passing them as a list to ``type``. If the type files have different numbers of header lines, pass a list of header line numbers with ``header``. The header list must be of length 1 or identical to the length of ``type``.

The extension of the resulting file will be the same as that of the files being merged, if they are all the same. If not, will not add an extension. To change the default behaviour, set an ``ext`` parameter with the extension to use, *e.g.* ``fna``. If several types are being merged, if ``ext`` is a string, the string will be used for all types. For a different ``ext`` for each file type, use a list of strings, in the same order as the ``type`` parameter.

.. Attention:: If you split sample-scope fasta files with ``fasta_splitter`` or ``split_fasta`` modules, the new subsamples are stored with a ``source`` category, containing the sample name from which the subsample was produced. When merging back into the sample scope, use ``scope: group`` and ``category: source``.

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* A table file in any slot:

    * ``sample_data[<sample>][<file.type>]``

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Puts output files in the following slot:

    * ``sample_data["project_data"][<file.type>]``

* Or, for merging by category, in the following slot:

    * ``sample_data[category_level][<file.type>]``


Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: Parameters that can be set:
    :header: "Parameter", "Values", "Comments"

    "type", "", "A file type that exists in all samples. Can also be a list of types, each one of which will be merged independently"
    "script_path", "", "Leave blank"
    "scope", "project|group", "Merge all samples into one project table, or merge sample tables by category."
    "category", "", "If ``scope`` is set to ``group``, you must specify the category by which to divide the samples for merging. The category must be a string containing one of the categories (columns) in the mapping file"
    "header","0","The number of header lines each table has. The header will be used for the complete table and all other headers will be removed. If there is no header line, set to 0 or leave out completely. **If set but not specified, will default to 1!**."
    "ext","","The extension to use for the merged file. If ``type`` is a list, ``ext`` will be used for all types unless ``ext`` itself is a list of the same length as ``type``."
    "add_filename", "", "If set, the source filename will be appended to each line in the resulting table."


Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Merge sample-scope tables into single project-scope table::

    merge_blast_tables:
        module:         merge_table
        base:           merge1
        script_path:
        scope:          project
        type:           blast.prot
        header:         0

Merge sample-scope tables into group-scope table, by category *country*::

    merge_blast_tables:
        module:         merge_table
        base:           merge1
        script_path:
        scope:          group
        category:       country
        type:           blast.prot
        header:         0


"""



import os
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept
from pprint import pprint as pp

__author__ = "Menachem Sklarz"
__version__ = "1.6.0"

class Step_merge_table(Step):
   
    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        # self.file_tag = "Bowtie_mapper"

        # Checking type exists, and converting to list
        if "type" not in self.params:
            raise AssertionExcept("You must supply a file type or list thereof to merge")
        else: # If passed by user, converting even single values to a list.
            if type(self.params["type"]) not in [str,list]:
                raise AssertionExcept("Type must be a string or list of strings.")
            if not isinstance(self.params["type"], list):
                self.params["type"] = [self.params["type"]]
        # Getting header parameter
        if "header" not in self.params:
            self.params["header"] = [0]
        elif not self.params["header"] and self.params["header"]!=0:
            self.params["header"] = [1]
        elif isinstance(self.params["header"], int):
            self.params["header"] = [self.params["header"]]
        elif isinstance(self.params["header"], list):
            for header in self.params["header"]:
                if not isinstance(header, int):
                    raise AssertionExcept("'header' must be an integer or a list of integers")
        else:
            raise AssertionExcept("'header' must be an integer or a list of integers")

        head_len = len(self.params["header"])
        type_len = len(self.params["type"])

        if type_len > 1 and head_len == 1:       # Many types, one header
            self.params["header"] = self.params["header"] * type_len
        elif type_len==head_len:
            pass
        else:
            raise AssertionExcept("'type' and 'header' parameters not the same length!")

        # Getting ext parameter
        if "ext" not in self.params:
            self.params["ext"] = [None] * type_len
        elif isinstance(self.params["ext"], str):
            self.params["ext"] = [self.params["ext"]] * type_len
        elif isinstance(self.params["ext"], list):
            for ext in self.params["ext"]:
                if not isinstance(ext, str):
                    raise AssertionExcept("'ext' must be a string or a list of strings")
        else:
            raise AssertionExcept("'ext' must be a string or a list of strings")

        ext_len = len(self.params["ext"])
        if type_len != ext_len:
            raise AssertionExcept("'type' and 'ext' parameters not the same length!")

    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """

        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
            for type in self.params["type"]:
                if type not in self.sample_data[sample]:
                    raise AssertionExcept("Type %s does not exist in sample" % type, sample)

        if "scope" not in self.params:
            self.write_warning("No 'scope' specified. Setting to 'project'")
            self.params["scope"] = "project"
        if self.params["scope"] == "project":
            pass
        elif self.params["scope"] == "group":
            if "category" not in self.params:
                raise AssertionExcept("For merging by group, you must supply a 'category' parameter.")
            if not isinstance(self.params["category"], str):
                raise AssertionExcept("'category' parameter must be a string.")
            for sample in self.sample_data["samples"]:  # Getting list of samples out of samples_hash
                try:
                    self.sample_data[sample]["grouping"][self.params["category"]]
                except KeyError:
                    raise AssertionExcept("category '{cat}' not defined for sample".format(cat=self.params["category"]),
                                          sample)
        else:
            raise AssertionExcept("'scope' must be either 'group' or 'project'")
        
    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        pass

    def build_scripts(self):

        if self.params["scope"] == "project":
            self.build_scripts_project()
        else:
            self.build_scripts_group()

    def build_scripts_project(self):
        """ This is the actual script building function
            Most, if not all, editing should be done here 
            HOWEVER, DON'T FORGET TO CHANGE THE CLASS NAME AND THE FILENAME!
        """
        
        # Each iteration must define the following class variables:
            # self.spec_script_name
            # self.script

        # for type in self.params["type"]:
        for i in range(len(self.params["type"])):
            type_i = self.params["type"][i]
            header_i = self.params["header"][i]
            ext_i = self.params["ext"][i]

            # If ext is not user defined, get exts of to-be-merged files. If all the same, use for merged file.
            if ext_i is None:
                base_exts = list(set([os.path.splitext(self.sample_data[sample][type_i])[1] for sample in self.sample_data["samples"]]))
                if len(base_exts)==1:
                    ext_i = base_exts[0].lstrip(".")

            # Name of specific script:
            self.spec_script_name = self.jid_name_sep.join([self.step,self.name,self.sample_data["Title"],type_i])
            self.script = ""
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(self.base_dir)

            # Define location and prefix for output files:
            output_fn = "{filename}".format(filename = ".".join(item for item
                                                                in [self.sample_data["Title"],type_i,ext_i]
                                                                if item))

            self.script += """\
SKIP={skip}
HEADER={header}
awk -v header="$HEADER" -v skip="$SKIP" \\
    '{fn_function}BEGIN{{ORS=""; headerline=0; skipline=0;}} 
    FNR==1 {{headerline=0; skipline=0}} {comment_str}
    skipline<skip {{print ""; skipline=skipline+1; next}};
    headerline<header {{
        if (NR==FNR) {{{line2print}}}
        headerline=headerline+1; next
    }}
    headerline>=header {{{line2print}}}' \\
    {infiles} \\
    > {outfile}

""".format(skip=self.params["skip"] if "skip" in self.params else 0,
           header=header_i, #.params["header"] if "header" in self.params else 0,
           line2print='printf("%s\\t%s\\n",basename(FILENAME),$0)'
                            if "add_filename" in self.params
                            else 'printf("%s\\n",$0)',
           comment_str='\n    /^{comm}/ {{print ""; next}};'.format(comm=self.params["comment"])
                                                    if "comment" in self.params
                                                    else '',
           fn_function='function basename(file) {{sub(".*/", "", file); return file}}\n    '
                                    if "add_filename" in self.params
                                    else '',
           infiles= ' \\\n    '.join([self.sample_data[sample][type_i] for sample in self.sample_data["samples"]]),
           outfile=use_dir+output_fn)

            self.sample_data["project_data"][type_i] = "%s%s" % (self.base_dir, output_fn)
            self.stamp_file(self.sample_data["project_data"][type_i])

            # Move all files from temporary local dir to permanent base_dir
            self.local_finish(use_dir,self.base_dir)
            self.create_low_level_script()

    def build_scripts_group(self):
        """ This is the actual script building function
            Most, if not all, editing should be done here
            HOWEVER, DON'T FORGET TO CHANGE THE CLASS NAME AND THE FILENAME!
        """

        # Create slots for levels of category, including types derived from samples.
        self.create_group_slots(self.params["category"])

        # Get levels of category
        cat_levels = self.get_category_levels(self.params["category"])
        for cat_lev in cat_levels:

            # for type_i in self.params["type"]:
            for i in range(len(self.params["type"])):
                type_i = self.params["type"][i]
                header_i = self.params["header"][i]
                ext_i = self.params["ext"][i]

                # If ext is not user defined, get exts of to-be-merged files. If all the same, use for merged file.
                if ext_i is None:
                    base_exts = list(set([os.path.splitext(self.sample_data[sample][type_i])[1]
                                          for sample
                                          in self.get_samples_in_category_level(self.params["category"],
                                                                                cat_lev)]))


                    if len(base_exts) == 1:
                        ext_i = base_exts[0].lstrip(".")

                self.spec_script_name = self.jid_name_sep.join([self.step, self.name, cat_lev, type_i])
                self.script = ""

                # This line should be left before every new script. It sees to local issues.
                # Use the dir it returns as the base_dir for this step.
                group_dir = self.make_folder_for_sample(cat_lev)
                use_dir = self.local_start(group_dir)

                # Define location and prefix for output files:
                output_fn = ".".join(item for item in [cat_lev, type_i, ext_i] if item)

                self.script += """\
SKIP={skip}
HEADER={header}
awk -v header="$HEADER" -v skip="$SKIP" \\
    '{fn_function}BEGIN{{ORS=""; headerline=0; skipline=0;}} 
    FNR==1 {{headerline=0; skipline=0}} {comment_str}
    skipline<skip {{print ""; skipline=skipline+1; next}};
    headerline<header {{
        if (NR==FNR) {{{line2print}}}
        headerline=headerline+1; next
    }}
    headerline>=header {{{line2print}}}' \\
    {infiles} \\
    > {outfile}

""".format(skip=self.params["skip"] if "skip" in self.params else 0,
           header=header_i, #.params["header"] if "header" in self.params else 0,
           line2print='printf("%s\\t%s\\n",basename(FILENAME),$0)' if "add_filename" in self.params else 'printf("%s\\n",$0)',
           comment_str='\n    /^{comm}/ {{print ""; next}};'.format(comm=self.params["comment"]) if "comment" in self.params else '',
           fn_function='function basename(file) {{sub(".*/", "", file); return file}}\n   ' if "add_filename" in self.params else '',
           infiles=' \\\n    '.join(
               [self.sample_data[sample][type_i]
                for sample
                in self.get_samples_in_category_level(self.params["category"], cat_lev)]),
           outfile=use_dir + output_fn)

                self.sample_data[cat_lev][type_i] = "%s%s" % (group_dir, output_fn)
                self.stamp_file(self.sample_data[cat_lev][type_i])

                # print self.script
                # Move all files from temporary local dir to permanent base_dir
                self.local_finish(use_dir, group_dir)
                self.create_low_level_script()

            # Setting new sample types:
            self.sample_data[cat_lev]["type_i"] = self.determine_sample_types(cat_lev,self.sample_data[cat_lev])
        # Setting new sample names to category levels.
        # From now on, these are the new samples.
        self.sample_data["samples"] = cat_levels
