# -*- coding: UTF-8 -*-
""" 
``STAR_LoadRemoveGenome``
-----------------------------------------------------------------
:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

A module for loading a STAR genome into RAM for use by subsequent STAR mapping jobs.

.. Note:: This module saves memory and time. Set parameter ``--genomeLoad`` in the STAR mapping instance to ``LoadAndKeep``.
   This will load the genome once into memory and use it repeatedly for all instances executed on the same node.
   When all mapping jobs are completed, Scripts produced by this instance will remove the genome from RAM for all
   nodes used.

.. Tip:: Make sure you set the ``node`` parameter in ``qsub_params`` to all the nodes in use by the base ``STAR_mapper`` instance.

.. Attention:: Currently defined for project-scope or external genomes only. Not used for sample-scope genomes.

.. Note:: Loading a genome is not really required. It will be loaded by the first instance of STAR.

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* A STAR genome in:

    * ``sample_data["STAR.index"]``

Alternatively, a STAR genome index can be passed with the ``--genomeDir`` parameter.


Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

No output is created

Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "genome", "load|remove", "Load or remove genome from RAM"
    "qsub_params:node", "", "Nodes on which to load/unload genome"
    "scope", "project | sample", "The scope from which to take the genome directory. Currently not in use"

Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**For external index:**

::

    STAR_remove_genome:
        module:             STAR_LoadRemoveGenome
        base:               STAR_map
        script_path:        '{Vars.paths.STAR}STAR'
        genome:             remove
        qsub_params:
            queue:          queue.q
            node:           {Vars.nodes}
        redirects:
            --genomeDir:    /path/to/STAR/genome_directory
            
            
**For project STAR index:**

::

    STAR_remove_genome:
        module:             STAR_LoadRemoveGenome
        base:               STAR_map
        script_path:        '{Vars.paths.STAR}STAR'
        genome:             remove
        qsub_params:
            queue:          queue.q
            node:           {Vars.nodes}
    
References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Dobin, A., Davis, C.A., Schlesinger, F., Drenkow, J., Zaleski, C., Jha, S., Batut, P., Chaisson, M. and Gingeras, T.R., 2013. **STAR: ultrafast universal RNA-seq aligner**. *Bioinformatics*, 29(1), pp.15-21.

"""


import os, re
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept

__author__ = "Menachem Sklarz"
__version__ = "1.2.0"

class Step_STAR_LoadRemoveGenome(Step):

    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances

        if "qsub_params" not in self.params or \
            "node" not in self.params["qsub_params"] or \
            self.params["qsub_params"]["node"] is None:
            raise AssertionExcept("You need to define nodes on which to load the genome with 'node' in 'qsub_params'")

        if "genome" not in self.params:
            raise AssertionExcept("Please set 'genome' to 'load' or 'remove'")

        if self.params["genome"] not in ["load","remove"]:
            raise AssertionExcept("'genome' should be either 'load' or 'remove'")

    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """


        if "--genomeDir" in self.params["redir_params"]:
            if "scope" in self.params:
                raise AssertionExcept("Both 'scope' and '--genomeDir' specified!\n")

        else:
            # If --genomeDir is not set, try using project genomeDir 
            if "STAR.index" in self.sample_data["project_data"]:
                self.params["redir_params"]["--genomeDir"] = self.sample_data["project_data"]["STAR.index"]
            else:
                raise AssertionExcept("No reference exists at 'project' scope. Do you have a STAR_builder step defined?")

    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        pass

    def build_scripts(self):
        """ This is the actual script building function
            Most, if not all, editing should be done here 
            HOWEVER, DON'T FORGET TO CHANGE THE CLASS NAME AND THE FILENAME!
        """
        
        nodes_list = self.params["qsub_params"]["node"]
        for node in nodes_list:      # Getting list of samples out of samples_hash
            # Make a dir for the current sample:
            node_dir = self.make_folder_for_sample(node)

            # Name of specific script:
            self.spec_script_name = self.jid_name_sep.join([self.step,self.name,node])
            self.script = ""
            
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(node_dir)
            output_prefix = node + "_STAR_memory"

            self.script =  """
{script_path} \\
    --genomeLoad {LoadRemove}  \\
    --genomeDir {genomeDir} \\
    --outFileNamePrefix {outprefix}
""".format(script_path=self.params["script_path"],
           genomeDir=self.params["redir_params"]["--genomeDir"], 
           outprefix=use_dir + output_prefix,
           LoadRemove="Remove" if self.params["genome"]=="remove" else "LoadAndExit")

            self.params["qsub_params"]["node"] = [node]
            self.create_low_level_script()

        # Reset node list to list of nodes.
        self.params["qsub_params"]["node"] = nodes_list
