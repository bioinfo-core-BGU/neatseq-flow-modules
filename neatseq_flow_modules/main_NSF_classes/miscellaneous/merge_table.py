# -*- coding: UTF-8 -*-
""" 
``merge_table``
------------------------


:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A module for merging sample tables into a single project wide table. The table can be with or without a header line. 

Can be used for merging fasta and fastq files as well.

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* A table file in any slot:

    * ``sample_data[<sample>][<file.type>]``

    

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Puts output files in the following slot:
        
    * ``sample_data[<file.type>]``




Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: Parameters that can be set:
    :header: "Parameter", "Values", "Comments"

    "type", "", "A file type that exists in all samples."
    "script_path", "", "Leave blank"
    "header","0","The number of header lines each table has. The header will be used for the complete table and all other headers will be removed. If there is no header line, set to 0 or leave out completely. **If not specified but not set, will default to 1!**"



Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    merge_blast_tables:
        module:         merge_table
        base:           merge1
        script_path:    
        type:           [blast,blast.prot]
        header:         1


"""



import os
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.1.0"

class Step_merge_table(Step):
   
    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        # self.file_tag = "Bowtie_mapper"
        
        if "header" not in self.params:
            self.params["header"] = 0
        elif not self.params["header"]:
            self.params["header"] = 1
        else:
            pass
        
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """

        type = set()
        
        if "type" not in self.params:
            raise AssertionExcept("You must supply a file type or list thereof to merge")
        else: # If passed by user, converting even single values to a list.
            if not isinstance(self.params["type"], list):
                self.params["type"] = [self.params["type"]]

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
                    raise AssertionExcept("'category' {cat} not defined for sample".format(cat=self.params["category"]),
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

        for type in self.params["type"]:
            # Name of specific script:
            self.spec_script_name = self.jid_name_sep.join([self.step,self.name,self.sample_data["Title"],type])
            self.script = ""
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(self.base_dir)

            # Define location and prefix for output files:
            output_fn = "{filename}".format(filename = ".".join([self.sample_data["Title"],type]))


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
           header=self.params["header"] if "header" in self.params else 0,
           line2print='printf("%s\\t%s\\n",basename(FILENAME),$0)'
                            if "add_filename" in self.params
                            else 'printf("%s\\n",$0)',
           comment_str='\n    /^{comm}/ {{print ""; next}};'.format(comm=self.params["comment"])
                                                    if "comment" in self.params
                                                    else '',
           fn_function='function basename(file) {{sub(".*/", "", file); return file}}\n    '
                                    if "add_filename" in self.params
                                    else '',
           infiles= ' \\\n    '.join([self.sample_data[sample][type] for sample in self.sample_data["samples"]]),
           outfile=use_dir+output_fn)

            self.sample_data["project_data"][type] = "%s%s" % (self.base_dir, output_fn)
            self.stamp_file(self.sample_data["project_data"][type])

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

            for type in self.params["type"]:

                self.spec_script_name = self.jid_name_sep.join([self.step, self.name, cat_lev, type])
                self.script = ""

                # This line should be left before every new script. It sees to local issues.
                # Use the dir it returns as the base_dir for this step.
                group_dir = self.make_folder_for_sample(cat_lev)
                use_dir = self.local_start(group_dir)

                # Define location and prefix for output files:
                output_fn = ".".join([cat_lev, type])

                self.script += """\
SKIP={skip}
HEADER={header}
awk -v header="$HEADER" -v skip="$SKIP" \\
    '{fn_function} BEGIN{{ORS=""; headerline=0; skipline=0;}} 
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
           header=self.params["header"] if "header" in self.params else 0,
           line2print='printf("%s\\t%s\\n",basename(FILENAME),$0)'
                            if "add_filename" in self.params
                            else 'printf("%s\\n",$0)',
           comment_str='\n    /^{comm}/ {{print ""; next}};'.format(comm=self.params["comment"])
           if "comment" in self.params
           else '',
           fn_function='function basename(file) {{sub(".*/", "", file); return file}}\n   '
           if "add_filename" in self.params
           else '$0,"\\n"',
           infiles=' \\\n    '.join(
               [self.sample_data[sample][type]
                for sample
                in self.get_samples_in_category_level(self.params["category"], cat_lev)]),
           outfile=use_dir + output_fn)

                self.sample_data[cat_lev][type] = "%s%s" % (group_dir, output_fn)
                self.stamp_file(self.sample_data[cat_lev][type])

                # print self.script
                # Move all files from temporary local dir to permanent base_dir
                self.local_finish(use_dir, group_dir)
                self.create_low_level_script()

            # Setting new sample types:
            self.sample_data[cat_lev]["type"] = self.determine_sample_types(self.sample_data[cat_lev])
        # Setting new sample names to category levels.
        # From now on, these are the new samples.
        self.sample_data["samples"] = cat_levels