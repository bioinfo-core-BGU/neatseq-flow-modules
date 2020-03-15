# -*- coding: UTF-8 -*-
""" 
``Tree_plot``
-----------------------

:Authors: Liron Levin
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

.. Note:: This module was developed as part of a study led by Dr. Jacob Moran Gilad

Short Description
~~~~~~~~~~~~~~~~~~~~~~~~~
    A module for plotting tree file in newick format together with MetaData information and possible additional matrix information.
    
Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    * A tree file in newick format in:
        ``self.sample_data["project_data"]["newick"]``
    * Tab delimited file with samples names in one of the columns from:
        ``self.sample_data["project_data"]["MetaData"]``
        ``self.sample_data["project_data"]["results"]``
        or from external file.
    

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    * Generate pdf file of the tree with the MetaData information:


Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "",  "", ""
    

Comments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    *  The following R packages are required:
        ``optparse``
        ``ape``
        ``ggtree``
        ``openxlsx``

Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    Step_Name:                            # Name of this step
        module: Tree_plot                 # Name of the used module
        base:                             # Name of the step [or list of names] to run after and generate a Tree plot [must be after a tree making step]
                                          # If more then one base is specified: the first overwrite the other bases overlapped slots  
        script_path:                      # Command for running the Tree plot script
                                          # If this line is empty or missing it will try using the module's associated script
        iterate_on_bases:                 # If set will iterate over the step's bases and generate a plot for each base. 
        tree_by_heatmap:                  # Generate additional tree using Hierarchical Clustering of the heatmap
        redirects:
            --layout:                     # Tree layout [fan or rectangular (default)]
            --Meta_Data:                  # Path to tab-delimited Meta Data file with header line. 
                                          # If this line is empty or missing it will try searching for results data.
            --M_Excel:                    # If the Meta_Data input is an Excel file indicate the sheet name to use
            --ID_field:                   # Column name in the Meta Data file for IDs found in the tips of the tree
            --cols_to_use:                # Columns in the Meta Data file to use and the order from the center up  
            --open.angle:                 # Tree open angle.
            --branch.length:              # Don't use branch length [cladogram]
            --conect.tip:                 # Connect the tip to its label
            --pre_spacer:                 # Space before the label text [default=0.05]
            --post_spacer:                # Space after the label text [default=0.01]
            --OTU:                        # Column name in the Meta Data file to use as OTU annotation
            --labels:                     # Use branch length labels
            --Tip_labels:                 # Show tip labels
            --heatmap:                    # Path to Data file to generate a heatmap 
                                          # If this line is empty it will try searching for results data.
            --H_Excel:                    # If the heatmap input is an Excel file indicate the sheet name to use
            --heatmap_cell_border:        # Color of heatmap cell border [default='white']
            --heatmap_lowest_value:       # Color of heatmap lowest value [default='white']
            --heatmap_highest_value:      # Color of heatmap highest value [default='red']
            --cols_to_use_heatmap:        # Columns in the heatmap Data file to use and the order from the center up
            --ID_heatmap_field:           # Column name for IDs found in the tips of the tree in the heatmap Data file
            --heatmap_variable:           # Use only variable columns in the heatmap
            --heatmap_count_by_sep:       # Count the sep in each cell to generate the values for the heatmap
            --heatmap_HC_dist:            # The heatmap Hierarchical Clustering dist method
            --heatmap_HC_agg:             # The heatmap Hierarchical Clustering agglomeration method

"""
import os
import sys
import re
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Liron Levin"
__version__= "1.2.0"
# plot tree file in newick format together with MetaData information 

class Step_Tree_plot(Step):

    def step_specific_init(self):
        """ Called on intiation
            Good place for parameter testing.
            Wrong place for sample data testing
        """
        self.shell = "bash"
        self.file_tag = ".pdf"
        import inspect
        self.module_location=os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))

    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        
        pass
        
    def create_spec_preliminary_script(self):
        """ Add script to run BEFORE all other steps
        """
        
        pass
        
    def build_scripts(self):
        """ This is the actual script building function
            Most, if not all, editing should be done here 
            HOWEVER, DON'T FORGET TO CHANGE THE CLASS NAME AND THE FILENAME!
        """
        
        # Name of specific script:
        self.spec_script_name = self.set_spec_script_name()
        self.script = ""
        sample_dir=self.base_dir
        # This line should be left before every new script. It sees to local issues.
        # Use the dir it returns as the base_dir for this step.
        use_dir = self.local_start(sample_dir)
        tree_by_heatmap_flag=True
        
        if self.params["script_path"]==None: 
            if "Tree_plot.R" in os.listdir(self.module_location): 
                self.params["script_path"]= "Rscript  %s "  % os.path.join(self.module_location,"Tree_plot.R")
            else:
                raise AssertionExcept("The file %s is not found in the Tree_plot module directory" % "Tree_plot.R" )
        
        if 'iterate_on_bases' not in list(self.params.keys()):
            if self.params["script_path"]!=None:  
                #Tree_plot main command
                if "env" in list(self.params.keys()):
                    self.script +="env %s  \\\n" % self.params["env"]
                self.script += "%s \\\n\t"  % self.params["script_path"]
                if ('--tree' or '-t') not in self.params["redir_params"]:
                    if 'newick' in self.sample_data["project_data"]:
                        self.script += "--tree %s \\\n\t" % self.sample_data["project_data"]['newick']
                    elif 'tree_by_heatmap' not in list(self.params.keys()):
                        raise AssertionExcept("\tThere is no project level tree information found [No newick slot] \n"  )
                    else:
                        self.script += "--tree_by_heatmap \\\n\t" 
                        tree_by_heatmap_flag=False
                        
                        
                if '--Meta_Data' not in self.params["redir_params"]:
                    if 'MetaData' in list(self.sample_data["project_data"].keys()):
                        self.script += "--Meta_Data %s \\\n\t" % self.sample_data["project_data"]['MetaData']  
                    elif 'results' in list(self.sample_data["project_data"].keys()):
                        self.script += "--Meta_Data %s \\\n\t" % self.sample_data["project_data"]['results']
                    else:
                        raise AssertionExcept("\tThere is no project level MetaData information found [No MetaData/results slots] \n"  )
                
                if '--heatmap' in self.params["redir_params"]:
                    if self.params["redir_params"]['--heatmap']==None:
                        if 'results' in list(self.sample_data["project_data"].keys()):
                            self.params["redir_params"]['--heatmap']= self.sample_data["project_data"]['results']
                        else:
                            raise AssertionExcept("\t You must specify external heatmap file \n"  )
                
                for par in self.params["redir_params"]:
                    if par not in ["-o","--output"]:
                        if self.params["redir_params"][par]!=None:
                            self.script += "%s %%s \\\n\t"  % par \
                                                            % self.params["redir_params"][par]
                        else:
                            self.script += "%s  \\\n\t"  % par 
                
                if tree_by_heatmap_flag:
                    self.script += "--output %s \n\n" % os.sep.join([use_dir.rstrip(os.sep),self.name+self.file_tag])
                else:
                    self.script += "--output %s \n\n" % os.sep.join([use_dir.rstrip(os.sep),'Heatmap_HC_tree'+self.file_tag])
        else:
            sample_data=self.sample_data
            if self.params["script_path"]!=None:
                for base in self.get_base_step_list():
                    if len(self.get_base_step_list())>1:
                        sample_data=self.get_base_sample_data()[base.name]
                    #Tree_plot main command
                    if "env" in list(self.params.keys()):
                        self.script +="env %s  \\\n" % self.params["env"]
                    self.script += "%s \\\n\t"  % self.params["script_path"]
                    if ('--tree' or '-t') not in self.params["redir_params"]:
                        if 'newick' in list(sample_data["project_data"].keys()):
                            self.script += "--tree %s \\\n\t" % sample_data["project_data"]['newick']
                        else:
                            raise AssertionExcept("\tThere is no project level tree information found in base step : %s [No newick slot] \n"  % base.name  )
                    else:
                        raise AssertionExcept("\t you can't give external tree file when iterate_on_bases is set !!!\n"  )
                    if '--Meta_Data' not in self.params["redir_params"]:
                        if 'MetaData' in list(self.sample_data["project_data"].keys()):
                            self.script += "--Meta_Data %s \\\n\t" % self.sample_data["project_data"]['MetaData']  
                        elif 'results' in list(self.sample_data["project_data"].keys()):
                            self.script += "--Meta_Data %s \\\n\t" % self.sample_data["project_data"]['results']
                        else:
                            raise AssertionExcept("\tThere is no project level MetaData information found [No MetaData/results slots] \n"  )
                    
                    if '--heatmap' in self.params["redir_params"]:
                        if self.params["redir_params"]['--heatmap']==None:
                            if 'results' in list(self.sample_data["project_data"].keys()):
                                self.params["redir_params"]['--heatmap']= self.sample_data["project_data"]['results']
                            else:
                                raise AssertionExcept("\t You must specify external heatmap file \n"  )
                    
                    for par in self.params["redir_params"]:
                        if par not in ["-o","--output"]:
                            if self.params["redir_params"][par]!=None:
                                self.script += "%s %%s \\\n\t"  % par \
                                                                % self.params["redir_params"][par]
                            else:
                                self.script += "%s  \\\n\t"  % par 
                    
                   
                    self.script += "--output %s \n\n" % os.sep.join([use_dir.rstrip(os.sep),base.name+self.file_tag])
                
        if 'tree_by_heatmap' in list(self.params.keys()):
            if tree_by_heatmap_flag:
                if self.params["script_path"]!=None:  
                    #Tree_plot main command
                    if "env" in list(self.params.keys()):
                        self.script +="env %s  \\\n" % self.params["env"]
                    self.script += "%s \\\n\t"  % self.params["script_path"]
                    
                    self.script += "--tree_by_heatmap \\\n\t" 
                    
                    if '--Meta_Data' not in self.params["redir_params"]:
                        if 'MetaData' in list(self.sample_data["project_data"].keys()):
                            self.script += "--Meta_Data %s \\\n\t" % self.sample_data["project_data"]['MetaData']  
                        elif 'results' in list(self.sample_data["project_data"].keys()):
                            self.script += "--Meta_Data %s \\\n\t" % self.sample_data["project_data"]['results']
                        else:
                            raise AssertionExcept("\tThere is no project level MetaData information found [No MetaData/results slots] \n"  )
                    
                    if '--heatmap' in self.params["redir_params"]:
                        if self.params["redir_params"]['--heatmap']==None:
                            if 'results' in list(self.sample_data["project_data"].keys()):
                                self.params["redir_params"]['--heatmap']= self.sample_data["project_data"]['results']
                            else:
                                raise AssertionExcept("\t You must specify external heatmap file \n"  )
                    
                    for par in self.params["redir_params"]:
                        if par not in ["-o","--output"]:
                            if self.params["redir_params"][par]!=None:
                                self.script += "%s %%s \\\n\t"  % par \
                                                                % self.params["redir_params"][par]
                            else:
                                self.script += "%s  \\\n\t"  % par 
                    
                   
                    self.script += "--output %s \n\n" % os.sep.join([use_dir.rstrip(os.sep),'Heatmap_HC_tree'+self.file_tag])
                    
        # Wrapping up function. Leave these lines at the end of every iteration:
        self.local_finish(use_dir,sample_dir)       # Sees to copying local files to final destination (and other stuff)
                  
        
        self.create_low_level_script()
