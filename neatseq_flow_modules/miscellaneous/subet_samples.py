# -*- coding: UTF-8 -*-
""" 
``subet_samples``
------------------------


:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A module to subset the sample list by grouping category.

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~




Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: Parameters that can be set:
    :header: "Parameter", "Values", "Comments"

    "category", "", "The name of the category in grouping file to subset by"
    "levels", "", "A name or list of levels in category to **keep**"


Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    subset_samples_cat1:
        module:         subet_samples
        base:           merge1
        script_path:
        category:       Category1
        levels:         [A,B]


"""



import os
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.1.0"

class Step_subet_samples(Step):
   
    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        self.skip_scripts = True

        if type(self.params["category"]) not in [str,int]:
            raise AssertionExcept("Category must be a single value, not a list etc.")

    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        # Check all samples have grouping data
        bad_samples = [sample
                       for sample
                       in self.sample_data["samples"]
                       if "grouping" not in self.sample_data[sample]]
        if bad_samples:
            raise AssertionExcept("For some reason, sample '{smp}' does not have "
                                  "grouping data".format(smp=", ".join(bad_samples)))
        # Check category is in all samples
        bad_samples = [sample
                       for sample
                       in self.sample_data["samples"]
                       if self.params["category"] not in self.sample_data[sample]["grouping"]]
        if bad_samples:
            raise AssertionExcept(
                "Sample '{smp}' does not have '{cat}' category".format(smp=", ".join(bad_samples),
                                                                       cat=self.params["category"]))

        # Check all levels exist in category
        if not all(map(lambda level:  level in self.get_category_levels(self.params["category"]), self.params["levels"])):
            bad_levels = ", ".join([level for level in self.params["levels"] if level not in self.get_category_levels(self.params["category"])])
            raise AssertionExcept("Level '{lev}' is not defined for category '{cat}'".format(lev=bad_levels,
                                                                                         cat=self.params["category"]))
        # Create new sample list:
        samples=list()
        for level in self.params["levels"]:
            samples.extend(self.get_samples_in_category_level(self.params["category"],level))

        self.sample_data["samples"] = samples
        
    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        pass

    def build_scripts(self):

        pass


