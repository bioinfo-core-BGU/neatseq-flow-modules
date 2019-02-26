# -*- coding: UTF-8 -*-
""" 
``qiime_prep``
-------------------------------

:Authors: Menachem Sklarz
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

.. Note:: This module was developed as part of a study led by Dr. Jacob Moran Gilad

A module for preparing fastq reads for analysis with QIIME (1.9):

The reads stored in each sample are optinally joined and then set it a directory in such a way the downstream, QIIME's demult can concatenate the sequences while saving the sample of origin.

The directory will contain symbolic links to the files to be used by demult in the following step.




Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* fastq files in one of the following slots:

    * ``sample_data[<sample>]["fastq.F"]``
    * ``sample_data[<sample>]["fastq.R"]``
    * ``sample_data[<sample>]["fastq.S"]``
    

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Puts directory of links to files to use with QIIME:

    * ``self.sample_data["project_data"]["qiime.prep_links_dir"]``

* If join is performed:
    
    * puts the new joined reads in:

        * ``self.sample_data[<sample>]["fastq.J"]``

    * puts the unjoined forward reads in:

        * ``self.sample_data[<sample>]["fastq.F"]``

    * puts the unjoined reverse reads in:

        * ``self.sample_data[<sample>]["fastq.R"]``


    
Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "join", "none, join (or join_cat - not implemented)", "Wheather to join paired reads."
    "unjoined", "forward, reverse, both or none", "What to do with unjoined sequences? Use only forward, only reverse, both or none. If join is none, use this parameter to indicate which reads to take for analysis."
    "join_algo", "forward, reverse, both or none", "What to do with unjoined sequences?"
    "parameters", "", "Path to QIIME parameter file to be used downstream"

    
Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    q_prep_1:
        module: qiime_prep
        base: merge1
        script_path: /path/to/join_paired_ends.py
        join: join
        unjoined: forward
        parameters: /path/to/qiime_params.txt
        redirects:
            --pe_join_method: fastq-join

References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Caporaso, J.G., Kuczynski, J., Stombaugh, J., Bittinger, K., Bushman, F.D., Costello, E.K., Fierer, N., Peña, A.G., Goodrich, J.K., Gordon, J.I. and Huttley, G.A., 2010. "QIIME allows analysis of high-throughput community sequencing data". *Nature methods*, 7(5), pp.335-336.




"""

import os
import sys
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Menachem Sklarz"
__version__ = "1.6.0"


class Step_qiime_prep(Step):
    """ A class that defines a pipeline step name (=instance).
    """
    
    def step_specific_init(self):
        self.shell = "bash"      # Can be set to "bash" by inheriting instances
        self.file_tag = "qiime_prep"
        
                
        self.links_dir = self.base_dir + "links_for_demult" + os.sep
        if not os.path.isdir(self.links_dir):
            self.write_warning("Making dir for links (%s)\n" % self.links_dir)
            os.makedirs(self.links_dir) 
        else:
            self.write_warning("Dir %s exists for links\n" % self.links_dir)
        
        try:
            self.params["join"]
        except KeyError:
            raise AssertionExcept("You must have a 'join' parameter\n")

        try:
            self.params["unjoined"]
        except KeyError:
            raise AssertionExcept("You must have an 'unjoined' parameter\n" )

        if self.params["join"].lower() == "none":  # Can be 'none', 'join' or 'join_cat', which means to just concatenate the F and R reads. Not implemented
            if self.params["unjoined"].lower() == "none":
                raise AssertionExcept("You can't set both 'join' and 'unjoined' to 'none'!")
        elif self.params["join"].lower() == "join_cat":  # Can be 'none', 'join' or 'join_cat', which means to just concatenate the F and R reads. Not implemented
            pass
        else:
            if "-m" in self.params["redir_params"]:
                self.params["join_algo"] = self.params["redir_params"]["-m"]
            elif "--pe_join_method" in self.params["redir_params"]:
                self.params["join_algo"] = self.params["redir_params"]["--pe_join_method"]
            else:
                self.write_warning("No algoritm specified. Using fastq-join")
                self.params["join_algo"] = "fastq-join"
                    
            if self.params["join_algo"] not in ["fastq-join", "SeqPrep"]:
                raise AssertionExcept("You must define the --pe_join_method parameter. Either fastq-join or SeqPrep...\n")
                        
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        
        try:
            self.sample_data["project_data"]["qiime.parameters"] = self.params["parameters"]
        except KeyError:
            self.write_warning("You did not supply a parameter file.\n" )
        try:
            self.sample_data["project_data"]["qiime.mapping"] = self.params["mapping"]
        except KeyError:
            self.write_warning("You did not supply a QIIME mapping file.\n" )
        # else:
            # if not os.path.exists(self.params["parameters"]):
                # self.write_warning("Note! The parameter file does not exist\n")
        

        
    def create_spec_wrapping_up_script(self):
        """ Add stuff to check and agglomerate the output data
        """
        
        pass
      
      
    def build_scripts(self):
        

        
        # Each iteration MUST DEFINE the following class variables:
            # spec_script_name 
            # script
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
            # General comment: If there is a parallel routine for each direction (forward, reverse), add this loop	
            
            # Name of specific script:
            self.spec_script_name = self.set_spec_script_name(sample)
            
            # Init script itself
            self.script = ""
            
            # Make a dir for the current sample:
            sample_dir = self.make_folder_for_sample(sample)
            # This line should be left before every new script. It sees to local issues.
            # Use the dir it returns as the base_dir for this step.
            use_dir = self.local_start(sample_dir)


            # Make dir in links_for_demult folder:
            if not os.path.isdir(self.links_dir + sample):
                os.makedirs(self.links_dir + sample) 

            # Define a set of directions to use for the sample
            directions_to_use = set()
            if self.sample_data[sample]["type"] == "SE":
                directions_to_use |= {"fastq.S"}
            else:   # PE or mixed
                if self.params["join"].lower() == "join":
                    directions_to_use |= {"fastq.J"}
                if self.params["unjoined"].lower() == "forward":  # Options: forward, reverse, both or none
                    directions_to_use |= {"fastq.F","fastq.S"}        # Assumption: in case of mixed PE and SE, if you want the forward then you want the single sequences too.
                elif self.params["unjoined"].lower() == "reverse":  # Options: forward, reverse, both or none
                    directions_to_use |= {"fastq.R","fastq.S"}        # Assumption: in case of mixed PE and SE, if you want the reverse then you want the single sequences too.    
                elif self.params["unjoined"].lower() == "both":  # Options: forward, reverse, both or none
                    directions_to_use |= {"fastq.F","fastq.R","fastq.S"}
                else:  # =="none"
                    if not directions_to_use:  # directions_to_use is empty!
                        raise RuntimeError("You can't pass 'none' to both 'join' and 'unjoined' parameters in step %s" % self.name)
            
            # If no join is required then only make links:
            if self.sample_data[sample]["type"] == "SE" or self.params["join"].lower() == "none":
                
                # directions contains existing files that appear in the sets to the right of the "&"
                directions = set(self.sample_data[sample].keys()) & directions_to_use  

                for direction in directions:
                    # print STDERR "$direction, $sample-->".$samples_hash->{$name}->{$sample}->{$direction}."\n";
                    link_name = "".join([self.links_dir + sample + os.sep, ".".join([sample, direction[6], "fastq"])])

                    if os.path.exists(link_name):
                        self.write_warning("Link $link_name exists. Will overwrite when script is executed!!!!\n" )
                    cmd = "ln -sf %s %s" % (self.sample_data[sample][direction], link_name)
                    
                    self.script += cmd + "\n\n";
                    
            
            else:
            # If join required, then there are 2steps: run join and link resuults in links4demult
            ############# 1
            # Add to $script joining code
                    
                self.script += self.get_script_const()        # Gets the "env", "script_path" and "redir_params" part of the script which is always the same...
                # if "env" in self.params.keys():         # Add optional environmental variables.
                    # self.script += "env %s \\\n\t" % self.params["env"]
                # self.script += "%s \\\n\t" % self.params["script_path"]
                # for key in self.params["redir_params"].keys():
                    # self.script += "%s %s \\\n\t" % (key,self.params["redir_params"][key])
                self.script += "-f " + self.sample_data[sample]["fastq.F"] + " \\\n\t";
                self.script += "-r " + self.sample_data[sample]["fastq.R"] + " \\\n\t";
                self.script += "-o " + use_dir + "\n\n";
        
                # self.sample_data[sample]["sample_dir"] = sample_dir
            
                
                ############# 2
                # Pointing to resulting files in sample structure:
                if self.params["join_algo"] == "fastq-join":
                    self.sample_data[sample]["fastq.J"] = sample_dir + "fastqjoin.join.fastq";
                    self.sample_data[sample]["fastq.F"] = sample_dir + "fastqjoin.un1.fastq";
                    self.sample_data[sample]["fastq.R"] = sample_dir + "fastqjoin.un2.fastq";
            
                elif self.params["join_algo"] == "SeqPrep":
                # Define the file names for SeqPrep. See qiime documentation
                    raise AssertionExcept("SeqPrep is not yet defined...\n")
            
                else:
                    raise AssertionExcept("You must define a join_algo. Either fastq-join or SeqPrep...\n")
            
                    ############# 3
                    # Leave space to define concatenation of R and F files:
                if self.params["join"] == "join_cat":
                    # Define qiime files as readsS (+ joined + catted)
                    raise AssertionExcept("join_cat not defined yet!!!\n")
                
                ############# 4
                # Putting links to final files 
                if self.params["join"] == "join":
                    #### 1aii. Making soft links (checking each to make sure it does not already exist:)
                    directions = set(self.sample_data[sample].keys()) & directions_to_use  
                    for direction in directions:
                        link_name = "".join([self.links_dir + sample + os.sep, ".".join([sample, direction[6], "fastq"])])
                        if os.path.exists(link_name):
                            self.write_warning("Link " + link_name + "exists. Will overwrite when script is executed!!!!\n" )

                        cmd = "ln -sf %s %s" % (self.sample_data[sample][direction], link_name)
                        self.script += cmd + "\n\n";
                    
                elif self.params["join"] == "join_cat": 
                    # Define what to do if we want join + cat.
                    raise AssertionExcept("join_cat not yet defined!!\n")

            self.sample_data["project_data"]["qiime.prep_links_dir"] = self.links_dir
            
            # Move all files from temporary local dir to permanent base_dir
            self.local_finish(use_dir,sample_dir)       # Sees to copying local files to final destination (and other stuff)

            
            # self.stamp_dir_files(sample_dir)
            
            self.create_low_level_script()
                    
            