#!/fastspace/bioinfo_apps/python-2.7_SL6/bin/python
# -*- coding: UTF-8 -*-
""" 
``Roary``
-----------------------

:Authors: Liron Levin
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.

.. Note:: This module was developed as part of a study led by Dr. Jacob Moran Gilad. The Bi_clustering R script was created by Eliad Levi 

Short Description
~~~~~~~~~~~~~~~~~~~

A module for running Roary on GFF files
    
Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* For each Sample, GFF file location in:
    * ``sample_data[<sample>]["GFF"]``
* If there is a GFF directory in the following slot, no new GFF directory will be created and ONLY the GFF files in this directory will be analysed.
    * ``sample_data["GFF_dir"]``
* If the search_GFF flag is on GFF files will be searched in the last base name directory

Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
* puts output GFF directory location in the following slots:
    * ``sample_data["GFF"]``
* puts output pan_genome results directory location in the following slots:
    * ``sample_data["pan_genome_results_dir"]``
* puts output pan_genome presence_absence_matrix file location in the following slots:
    * ``sample_data["presence_absence_matrix"]``
* puts output pan_genome clustered_proteins file location in the following slots:
    * ``sample_data["clustered_proteins"]``
* puts output GWAS directory location in the following slot:
    * ``sample_data["GWAS_results_dir"]``
* puts output Biclustering directory location in the following slot:
    * ``sample_data["Bicluster_results_dir"]``
* puts output Biclustering cluster file location in the following slot:
    * ``sample_data["Bicluster_clusters"]``
* puts output Gecko directory location in the following slot:
    * ``sample_data["Gecko_results_dir"]``
* puts Accessory genes or virulence/resistance hierarchical-clustering tree file in the following slot:
    * ``self.sample_data["newick"]``

Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "",  "", ""
    

Comments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*  This module was tested on:
    * ``Roary v3.10.2``
    * ``Roary v1.006924``
    * ``Scoary v1.6.11``
    * ``Scoary v1.6.9``
    * ``Gecko3``
* For the Bi_cluster analysis the following R packages are required:
    * ``optparse``
    * ``eisa``
    * ``ExpressionView``
    * ``openxlsx``
    * ``clusterProfiler``
    * ``org.Hs.eg.db``
* To plot the pan-genome matrix the following python packages are required:
    * ``pandas``
    * ``patsy``
    * ``seaborn``
    * ``matplotlib``
    * ``numpy``
    * ``scipy``
* For the scoary analysis the following python packages are required:
    * ``pandas``
* For the Gecko analysis the following python packages are required:
    * ``pandas``
    
.. Note:: If using conda environment with R installed, the R packages will be automatically installed inside the environment.


Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::


    Step_Name:                                   # Name of this step
        module: Roary                            # Name of the module used
        base:                                    # Name of the step [or list of names] to run after [must be after a GFF file generator step like Prokka]
        script_path:                             # Command for running the Roary script 
        env:                                     # env parameters that needs to be in the PATH for running this module
        qsub_params:                             
            -pe:                                 # Number of CPUs to reserve for this analysis
        virulence_resistance_tag:                # Use the name of the db used in prokka or use "VFDB" if you used the VFDB built-in Prokka module DB 
        search_GFF:                              # Search for GFF files?
        Bi_cluster:                              # Do Bi_cluster analysis using the Roary results, if empty or this line dose not exist will not do Bi_cluster analysis 
            --Annotation:                        # location of virulence annotation file to use to annotate the clusters or use "VFDB" if you used the VFDB built-in Prokka module DB
            --ID_field:                          # The column name in the MetaData file of the samples IDs
            --cols_to_use:                       # list of the MetaData columns to use to annotate the clusters  example: '"ST","CC","source","host","geographic.location","Date"'
            --metadata:                          # location of MetaData file to use to annotate the clusters
        plot:                                    # plot gene presence/absence matrix
            format:                              # The gene presence/absence matrix plot output format. example: pdf
            Clustering_method                    # The gene presence/absence matrix plot hierarchical-clustering method. example: ward
            Tree:                                # Save s tree in newick format of the 'Accessory' genes or the 'virulence_resistance_tag' genes hierarchical-clustering
                                                 # example: Tree: Accessory 
        scoary:
            script_path:                         # Command for running the scoary script, if empty or this line dose not exist will not run scoary 
            BH_cutoff:                           # Scoary BH correction for multiple testing cut-off
            Bonferroni_cutoff:                   # Scoary Bonferroni correction for multiple testing cut-off
            metadata_file:                       # location of MetaData file to use to create the scoary traits file
            metadata_samples_ID_field:           # The column name in the MetaData file of the sample's IDs
            traits_file:                         # Path to a traits file
            traits_to_pars:               # If a traits file is not provided use a list of conditions to create the scoary traits file from MetaData file. example:"source/=='blood'"  "source/=='wound'"
                                                 # Pairs of field and operator + value to convert to boolean traits: field_name1/op_value1 .. field_nameN/op_valueN Example: "field_1/>=val_1<val_2"    "feild_2/=='str_val'"
                                                 # A Filter can be used by FILTER_field_name1/FILTER_op_value1&field_name1/op_value1
                                                 # Note that Gecko can't run if the Bi_clustering was not run
        Gecko:
            script_path:                         # Command for running the Gecko script, if empty or this line dose not exist will not run Gecko
            -d:                                  # Parameters for running Gecko
            -s:                                  # Parameters for running Gecko
            -q:                                  # Parameters for running Gecko
        redirects:
            -k:                                  # Parameters for running Roary
            -p:                                  # Parameters for running Roary
            -qc:                                 # Parameters for running Roary
            -s:                                  # Parameters for running Roary
            -v:                                  # Parameters for running Roary
            -y:                                  # Parameters for running Roary

References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    * **Roary program**: Page, Andrew J., et al. "Roary: rapid large-scale prokaryote pan genome analysis." Bioinformatics 31.22 (2015): 3691-3693.‏
    * **Scoary program**: Brynildsrud, Ola, et al. "Rapid scoring of genes in microbial pan-genome-wide association studies with Scoary." Genome biology 17.1 (2016): 238.‏
    * **Gecko program**: Winter, Sascha, et al. "Finding approximate gene clusters with Gecko 3." Nucleic acids research 44.20 (2016): 9600-9610.‏

"""
import os
import sys
import re
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Liron Levin"
__version__= "1.2.0"
# Roary on GFF files

class Step_Roary(Step):
    
    
    
    def step_specific_init(self):
        """ Called on intiation
            Good place for parameter testing.
            Wrong place for sample data testing
        """
        self.shell = "bash"
        self.file_tag = ""
        import inspect
        self.module_location=os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))
        
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        
        if "search_GFF" not in self.params.keys():                
            for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
                # Testing for existance of GFF data
                assert ("GFF" in self.sample_data[sample].keys())or ("GFF_dir" in self.sample_data.keys()), \
                    "In %s:\tThere are no GFF Annotation files to use in this step.\n" % self.get_step_name()
        pass
        
    def create_spec_preliminary_script(self):
        """ Add script to run BEFORE all other steps
        """
        
        
        GFF_dir =''
        self.script = ""
        #Find the dir for the GFF files:
        if "search_GFF" in self.params.keys():
            #Make a dir for the GFF files:
            GFF_dir = self.make_folder_for_sample("GFF")
            look_in=os.sep.join([self.pipe_data["data_dir"].rstrip(os.sep),\
                            self.get_base_step_list()[-1].step , \
                            self.get_base_step_list()[-1].name])
            self.script +="find %s -name '*.gff' -exec cp -s {} %%s \; \n\n" % look_in % GFF_dir
            
        else:
            if "GFF_dir" in self.sample_data.keys():
                GFF_dir = self.sample_data["GFF_dir"]
            if len(GFF_dir)==0:
                #Make a dir for the GFF files:
                GFF_dir = self.make_folder_for_sample("GFF")
                for sample in self.sample_data["samples"]:
                    if self.sample_data[sample]["GFF"].endswith(".gff"):
                        self.script +="cp -s %s %%s \n\n" % self.sample_data[sample]["GFF"] % GFF_dir
                    else:
                        self.script +="cp -s %s*.gff %%s \n\n" % self.sample_data[sample]["GFF"] % GFF_dir
        
        self.sample_data["GFF_dir"]=GFF_dir
        pass
        
    def build_scripts(self):
        """ This is the actual script building function
            Most, if not all, editing should be done here 
            HOWEVER, DON'T FORGET TO CHANGE THE CLASS NAME AND THE FILENAME!
        """
        GFF_dir=self.sample_data["GFF_dir"]
        sample='Pan_Genome'
        # Name of specific script:
        self.spec_script_name = "_".join([self.step,self.name,sample])
        self.script = ""

        # Make a dir for the RESULTS:
        sample_dir = self.make_folder_for_sample(sample)
        
        # This line should be left before every new script. It sees to local issues.
        # Use the dir it returns as the base_dir for this step.
        #use_dir = self.local_start(sample_dir)
            
            
        # Define output filename 
        output_filename = "".join([sample_dir , sample ])
        #Roary main command
        self.script += self.get_script_const()
        self.script += " -f %s \\\n\t" % output_filename
        self.script += " %s*.gff \n\n" % GFF_dir
        
        # Adding the results data
        self.sample_data["pan_genome_results_dir"]=output_filename
        self.sample_data["presence_absence_matrix"]=os.path.join(output_filename , "gene_presence_absence.csv") 
        self.sample_data["pan_genome_clustered_proteins"]=os.path.join(output_filename , "clustered_proteins") 
        
        # Creating the plots
        if "plot" in self.params.keys():
            if "Roary_matrix_plot.py" in os.listdir(self.module_location):
                if "env" in self.params.keys():
                    if self.shell=="bash":
                        self.script +="export env %s  \\\n\t" % self.params["env"]
                    else:
                        self.script +="env %s  \\\n\t" % self.params["env"]
                self.script += "python  %s \\\n\t"  % os.path.join(self.module_location,"Roary_matrix_plot.py")
                if is_it_dict(self.params["plot"]):
                    if "format" in self.params["plot"].keys():
                        self.script += " --format %s \\\n\t" % self.params["plot"]["format"]
                    if "virulence_resistance_tag" in self.params.keys():
                        if self.params["virulence_resistance_tag"]=="VFDB":
                            self.script += " --tag %s \\\n\t" % "Virulence_Resistance.fasta:Virulence"
                        else:
                            self.script += " --tag %s \\\n\t" % self.params["virulence_resistance_tag"]
                        if 'Tree' in self.params["plot"].keys():
                            if self.params["plot"]['Tree']=='virulence_resistance_tag':
                                self.sample_data["newick"]=os.path.join(self.sample_data["pan_genome_results_dir"],'virulence_resistance.newick')
                    if 'Tree' in self.params["plot"].keys():
                        if self.params["plot"]['Tree']=='Accessory':  
                            self.sample_data["newick"]=os.path.join(self.sample_data["pan_genome_results_dir"],'Accessory.newick')
                            
                    if "Clustering_method" in self.params["plot"].keys():
                        self.script += " -C %s \\\n\t" % self.params["plot"]["Clustering_method"]
                self.script += " -O %s \\\n\t" % self.sample_data["pan_genome_results_dir"]
                self.script += " -P %s \n\n"   % self.sample_data["presence_absence_matrix"]
            else:
                raise AssertionExcept("The file %s is not found in the Roary module directory" % "Roary_matrix_plot.py" )
        
        
        # Pan-genome wide association studies using scoary
        scoary_traits_file=''     
        gene_presence_absence_file_loc=self.sample_data["presence_absence_matrix"]
        if "scoary" in self.params.keys():
            if is_it_dict(self.params["scoary"]):
                if self.params["scoary"]["script_path"]!=None:
                    # if the traits file is provided use it
                    if "traits_file" in self.params["scoary"].keys():
                        scoary_traits_file=self.params["scoary"]["traits_file"]                
                        # Creating the result dir
                        GWAS_dir = self.make_folder_for_sample("GWAS")
                    # if the a metadata file is provided use it to create traits file
                    elif "metadata_file" in self.params["scoary"].keys():
                        if "Traits_Parser.py" in os.listdir(self.module_location):
                            if "traits_to_pars" in self.params["scoary"].keys():
                                # Creating the result dir
                                GWAS_dir = self.make_folder_for_sample("GWAS")
                                if "env" in self.params.keys():
                                    if self.shell=="bash":
                                        self.script +="export env %s  \\\n\t" % self.params["env"]
                                    else:
                                        self.script +="env %s  \\\n\t" % self.params["env"]
                                self.script +="python %s \\\n\t"  % os.path.join(self.module_location,"Traits_Parser.py")
                                self.script +=" -M %s \\\n\t"  % self.params["scoary"]["metadata_file"]
                                self.script +=" -O %s \\\n\t"  % GWAS_dir
                                # This option will create new gene presence absence file with correct samples names and the are shared with the traits file
                                self.script +=" -P %s \\\n\t"  % gene_presence_absence_file_loc                        
                                if "metadata_samples_ID_field" in self.params["scoary"].keys():
                                    self.script +=" --S_MetaData %s \\\n\t"  % self.params["scoary"]["metadata_samples_ID_field"]
                                self.script +=" --Fields_val %s \n\n"  % self.params["scoary"]["traits_to_pars"]
                                scoary_traits_file=os.path.join(GWAS_dir,'Traits_file.csv')
                                # The new gene presence absence file is in the GWAS dir and it is the input for scoary
                                gene_presence_absence_file_loc=os.path.join(GWAS_dir , "gene_presence_absence.csv")
                        else:
                            raise AssertionExcept("The file %s is not found in the Roary module directory" % "Traits_Parser.py" )
                    if len(scoary_traits_file)>0:
                        if "env" in self.params.keys():
                            if self.shell=="bash":
                                self.script +="export env %s  \\\n\t" % self.params["env"]
                            else:
                                self.script +="env %s  \\\n\t" % self.params["env"]
                        self.script += "%s \\\n\t"  % self.params["scoary"]["script_path"]
                        self.script += " -o %s \\\n\t"  % GWAS_dir
                        self.script += " -g %s \\\n\t"  % gene_presence_absence_file_loc
                        self.script += " -t %s \\\n\t"  % scoary_traits_file
                        if ("use_cluster_tree" in self.params["scoary"].keys()) & ("plot" in self.params.keys()):
                            self.script += " -n %s \\\n\t"  % os.path.join(self.sample_data["pan_genome_results_dir"],"pangenome_matrix.newick")
                        else:
                            self.script += " -u  \\\n\t"   
                        if "Bonferroni_cutoff" in self.params["scoary"].keys():
                            self.script += " -c B -p %s \\\n\t"  % self.params["scoary"]["Bonferroni_cutoff"]
                        elif "BH_cutoff" in self.params["scoary"].keys():
                            self.script += " -c BH -p %s \\\n\t"  % self.params["scoary"]["BH_cutoff"]
                        if "permutations" in self.params["scoary"].keys():
                            self.script += " -e %s  \n\n"  % self.params["scoary"]["permutations"]
                        # Adding the results data
                        self.sample_data["GWAS_results_dir"]=GWAS_dir
        self.script +=" \n\n"

        if "Bi_cluster" in self.params.keys():
            if "Biclustering.R" in os.listdir(self.module_location):
                gene_presence_absence_file_loc=self.sample_data["presence_absence_matrix"]
                # Make a dir for the results file:
                bicluster_results_dir = self.make_folder_for_sample("Bicluster")         
                #Running the bicluster script
                if "env" in self.params.keys():
                    if self.shell=="bash":
                        self.script +="export env %s  \\\n\t" % self.params["env"]
                    else:
                        self.script +="env %s  \\\n\t" % self.params["env"]
                self.script += "Rscript  %s \\\n\t"  % os.path.join(self.module_location,"Biclustering.R")
                temp_self_script=""
                if is_it_dict(self.params["Bi_cluster"]):
                    for par in self.params["Bi_cluster"].keys():                       
                        if par=="--Roary_Results":
                            self.write_warning("The '--Roary_Results' parameter in the Roary Bi_clustering analysis is ignored")
                        elif par=="-o":
                            self.write_warning("The '-o' parameter in the Roary Bi_clustering analysis is ignored")
                        elif par=="--Annotation":
                            if self.params["Bi_cluster"][par]=="VFDB":
                                if "VFDB_unified_VF_category_clustered.tsv" in os.listdir(self.module_location): 
                                    temp_self_script +="%s  %%s \\\n\t" % par \
                                                                        % os.path.join(self.module_location,"VFDB_unified_VF_category_clustered.tsv")
                                else:
                                    raise AssertionExcept("The file %s is not found in the Roary module directory" % "VFDB_unified_VF_category_clustered.tsv" )
                            else:
                                temp_self_script +="%s  %%s \\\n\t" % par \
                                                               % self.params["Bi_cluster"][par]
                        elif len(par)>0:
                            if self.params["Bi_cluster"][par]!=None:
                                temp_self_script +="%s  %%s \\\n\t" % par \
                                                               % self.params["Bi_cluster"][par]
                            else:
                                temp_self_script +="%s  \\\n\t" % par 

                self.script +="--Roary_Results %s  \\\n\t" % gene_presence_absence_file_loc
                self.script +=temp_self_script
                self.script +="-o %s  \\\n\t" % bicluster_results_dir
                self.script +=" \n\n"
                self.sample_data["Bicluster_results_dir"]=bicluster_results_dir
                self.sample_data["Bicluster_clusters"]=os.path.join(bicluster_results_dir,"Bicluster_clusters")

                # Run Gecko gene clusters analysis based on the Bi_clustering analysis
                if "Gecko" in self.params.keys():
                    if is_it_dict(self.params["Gecko"]):
                        if "script_path" in self.params["Gecko"].keys():  
                            if (self.params["Gecko"]["script_path"]!=None)&(self.params["Gecko"]["script_path"]!=''):
                                if "GFF2Gecko3.py" in os.listdir(self.module_location):
                                    Bicluster_clusters=self.sample_data["Bicluster_clusters"] 
                                    # Make a dir for the results file:
                                    Gecko_results_dir = self.make_folder_for_sample("Gecko")         
                                    #Running the GFF2Gecko3 script
                                    if "env" in self.params.keys():
                                        if self.shell=="bash":
                                            self.script +="export env %s  \\\n\t" % self.params["env"]
                                        else:
                                            self.script +="env %s  \\\n\t" % self.params["env"]
                                    self.script += "python  %s \\\n\t"  % os.path.join(self.module_location,"GFF2Gecko3.py")
                                    if "-p" in self.params["redir_params"]:
                                        self.script += "-P  %s \\\n\t"  % self.params["redir_params"]["-p"]
                                    self.script += "-D  %s \\\n\t"  % GFF_dir
                                    self.script += "-C  %s \\\n\t"  % self.sample_data["pan_genome_clustered_proteins"]
                                    self.script += "-B  %s \\\n\t"  % Bicluster_clusters
                                    self.script += "-o  %s \n\n"    % os.path.join(Gecko_results_dir,"Gecko.cog")
                                    if "env" in self.params.keys():
                                        if self.shell=="bash":
                                            self.script +="export env %s  \\\n\t" % self.params["env"]
                                        else:
                                            self.script +="env %s  \\\n\t" % self.params["env"]
                                    temp_self_script=""
                                    Gecko_pars=list()
                                    for par in self.params["Gecko"].keys():
                                            Gecko_pars.append(par)
                                            if par=="-in":
                                                self.write_warning("The '-in' parameter in the Roary Gecko analysis is ignored")
                                            elif par=="-out":
                                                self.write_warning("The '-out' parameter in the Roary Gecko analysis is ignored")
                                            elif len(par)>0:
                                                if par!="script_path":
                                                    if self.params["Gecko"][par]!=None:
                                                        temp_self_script +="%s  %%s \\\n\t" % par \
                                                                                       % self.params["Gecko"][par]
                                                    else:
                                                        temp_self_script +="%s  \\\n\t" % par 
                                            
                                      
                                    self.script +="%s  \\\n\t" % self.params["Gecko"]["script_path"]
                                    
                                    if "-r" not in Gecko_pars:
                                        temp_self_script +="-r Reference_clusters \\\n\t"
                                    if "-s" not in Gecko_pars:
                                        temp_self_script +="-s 2 \\\n\t"
                                    if "-d" not in Gecko_pars:
                                        temp_self_script +="-d 7 \\\n\t"
                                    if "-q" not in Gecko_pars:
                                        temp_self_script +="-q 2 \\\n\t"
                                    if "-rO" not in Gecko_pars:
                                        temp_self_script +="-rO zippedPdfs showFiltered %s \\\n\t" % os.path.join(Gecko_results_dir ,"Clusters.zip" )
                                    else:
                                        temp_self_script +="-rO %s %%s \\\n\t" %  self.params["Gecko"]["-rO"]\
                                                                               %  os.path.join(Gecko_results_dir ,"Clusters" )
                                    self.script +="-in %s  \\\n\t" % os.path.join(Gecko_results_dir,"Gecko.cog")
                                    self.script +="-out %s  \\\n\t" % os.path.join(Gecko_results_dir,"Gecko.gck")
                                    self.script +=temp_self_script
                                    self.script +=" \n\n"
                                    self.sample_data["Gecko_results_dir"]=Gecko_results_dir                            
                                else:
                                    raise AssertionExcept("The file %s is not found in the Roary module directory" % "GFF2Gecko3.py" )
                            else:
                                self.write_warning("No %s running command found! will not run! " % "Gecko" )
            else:
                raise AssertionExcept("The file %s is not found in the Roary module directory" % "Biclustering.R" )

        
        for sample in self.sample_data["samples"]:      # Getting list of samples out of samples_hash
            # Store Roary result location:
            self.sample_data[sample]["pan_genome_results_dir"]=self.sample_data["pan_genome_results_dir"]

        # Wrapping up function. Leave these lines at the end of every iteration:
        #self.local_finish(use_dir,sample_dir)       # Sees to copying local files to final destination (and other stuff)
                  
        
        self.create_low_level_script()


def is_it_dict(dic):
    try:
        dic.keys()
        return True
    except:
        return False