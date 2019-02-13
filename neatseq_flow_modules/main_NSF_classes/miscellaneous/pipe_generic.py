# -*- coding: UTF-8 -*-
""" 
``pipe_generic`` 
----------------------------------------------------

:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

Description
~~~~~~~~~~~~

Name changed to Fillout_Generic

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
        raise AssertionExcept("'pipe_generic has been renamed to 'Fillout_Generic'. Sorry for the inconvenience.")


    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        pass

    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        pass


    def build_scripts(self):
        """ This is the actual script building function
            
        """
        pass
