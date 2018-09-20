# -*- coding: UTF-8 -*-
""" 
``del_type``
-----------------------------------------------------------------
:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A module for ...

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

..    "pipe", "", "Additional commands to be piped on the files before writing to file."

Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

        
"""

import os
import sys
import re
from neatseq_flow.PLC_step import Step,AssertionExcept

from  modules.global_defs import ZIPPED_EXTENSIONS, ARCHIVE_EXTENSIONS, KNOWN_FILE_EXTENSIONS



__author__ = "Menachem Sklarz"
__version__ = "1.1.0"



# A dict for conversion of types of sample data to positions in fasta structure:
fasta_types_dict = {"Nucleotide":"fasta.nucl","Protein":"fasta.prot"}
sam_bam_dict     = {"SAM":"sam", "BAM":"bam", "REFERENCE":"reference"}

class Step_del_type(Step):
    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        # self.file_tag = "merge"
        self.skip_scripts = True
        
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        if "type2del" not in self.params:
            raise AssertionExcept("You must pass a 'type2del' param!")
        
        type2del = self.params["type2del"]
        
        if "scope" not in self.params:
            raise AssertionExcept("You must pass a 'scope' param!")
        if self.params["scope"] == "sample":
            for sample in self.sample_data["samples"]:
                if type2del not in self.sample_data[sample]:
                    raise AssertionExcept("type %s does not exist for project." % type2del)
                else:
                    del self.sample_data[sample][type2del]
        elif self.params["scope"] == "project":
            if type2del not in self.sample_data:
                raise AssertionExcept("type %s does not exist for project." % type2del)
            else:
                del self.sample_data["project_data"][type2del]
        else:
            raise AssertionExcept("'scope' param must be 'sample' or 'project'")
        
        
    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        
        pass
      
    def build_scripts(self):
        
        return 
        # if 
        

                    