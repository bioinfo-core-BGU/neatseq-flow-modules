# -*- coding: UTF-8 -*-
""" 
``DeSeq2``
-----------------------

:Authors: Liron Levin
:Affiliation: Bioinformatics core facility
:Organization: National Institute of Biotechnology in the Negev, Ben Gurion University.


Short Description
~~~~~~~~~~~~~~~~~~~~~~~~~
    A module to preform:
    * Gene level differential expression using DeSeq2.
    * Gene annotation.
    * PCA plot.
    * Clustering of significant genes.
    * Heatmaps of significant genes by clusters.
    * Expression patterns plot by clusters
    * Enrichment analysis KEGG/GO.
    

Requires
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    * Search for count data in :
        `self.sample_data[<sample>]["RSEM"]`
        `self.sample_data[<sample>]["genes.counts"]`
        `self.sample_data[<sample>]["HTSeq.counts"]`
        `self.sample_data["project_data"]["results"]`


Parameters that can be set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: 
    :header: "Parameter", "Values", "Comments"
    :widths: 15, 10, 10

    "use_click",  "", "Will use the CLICK clustering program (Shamir et al. 2000)"
    
.. Note:: 
    If your using the **use_click** option, cite:
    Expander:
    Ulitsky I, Maron-Katz A, Shavit S, Sagir D, Linhart C, Elkon R, Tanay A, Sharan R, Shiloh Y, Shamir R. Expander: from expression microarrays to networks and functions. Nature Protocols Vol 5, pp 303 - 322, 2010
    Click:
    Shamir , R. and Sharan, R. CLICK: A Clustering Algorithm with Applications to Gene Expression Analysis. Proceedings ISMB 2000, pp.307-316 (2000)

    
Comments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    *  The following R packages are required:
        ``DESeq2``
        ``ggplot2``
        ``pheatmap``
        ``mclust``
        ``factoextra``
        ``cowplot``
        ``gridExtra``
        ``biomaRt``
        ``clusterProfiler``
        ``KEGGREST``
        ``scater``
        ``sva``
        ``rmarkdown``
        ``plotly``
        ``dt``
        ``xml2``
        ``dplyr``
        ``rcolorbrewer``
        ``colorspace``
        ``stringr``

.. Note:: 
    It is Possible to use CONDA to install all dependencies:
    
    .. code-block:: bash
        
        wget https://raw.githubusercontent.com/bioinfo-core-BGU/neatseq-flow-modules/master/neatseq_flow_modules/Liron/DeSeq2_module/DeSeq2_env_install.yaml
        conda env create -f DeSeq2_env_install.yaml
    
    Flow this `Tutorial <https://github.com/bioinfo-core-BGU/NeatSeq-Flow_Workflows/blob/master/DeSeq_Workflow/Tutorial.md>`_ for More Information.

Lines for parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    Step_Name:                              # Name of this step
        module: Proteomics                  # Name of the used module
        base:                               # Name of the step [or list of names] to run after with count results.
        script_path:                        # Command for running the a DeSeq2 script
                                            # If this line is empty or missing it will try using the module's associated script
        use_click:                          # Will use the CLICK clustering program (Shamir et al. 2000). 
        redirects:
            --SAMPLE_DATA_FILE:             # Path to Samples Information File
            --GENE_ID_TYPE:                 # The Gene ID Type i.e 'ENSEMBL'[for Bioconductor] OR 'ensembl_gene_id'/'ensembl_transcript_id' [for ENSEMBL]
            --Annotation_db:                # Bioconductor Annotation Data Base Name from https://bioconductor.org/packages/release/BiocViews.html#___OrgDb  
            --Species:                      # Species Name to Retrieve Annotation Data from ENSEMBL
            --KEGG_Species:                 # Species Name to Retrieve Annotation Data from KEGG
            --KEGG_KAAS:                    # Gene to KO file from KEGG KAAS [first column gene id, second column KO number]
            --Trinotate:                    # Path to a Trinotate annotation file in which the first column is the genes names
            --FILTER_SAMPLES:               # Filter Samples with Low Number of expressed genes OR with Small Library size using 'scater' package 
            --FILTER_GENES:                 # Filter Low-Abundance Genes using 'scater' package
            --NORMALIZATION_TYPE:           # The DeSeq2 Normalization Type To Use [VSD , RLOG] The Default is VSD
            --BLIND_NORM:                   # Perform Blind Normalization
            --DESIGN:                       # The Main DeSeq2 Design [ ~ Group ]
            --removeBatchEffect             # Will Remove Batch Effect from the Normalized counts data up to 2 
                                            # [using the limma package and only one using the sva package]
                                            # Batch Effect fields [from the Sample Data ] separated by , 
            --removeBatchEffect_method      # The method to Remove Batch Effect from the Normalized counts data using the limma or sva packages [sva is the default]
            --LRT:                          # The LRT DeSeq2 Design
            --ALPHA:                        # Significant Level Cutoff, The Default is 0.05
            --Post_statistical_ALPHA        # Post Statistical P-value Filtering
            --FoldChange:                   # Fold change Cutoff [testing for fold changes greater in absolute value], The Default is 1
            --Post_statistical_FoldChange   # Post Statistical Fold change Filtering
            --CONTRAST:                     # The DeSeq Contrast Design ["Group,Treatment,Control"] [Not For LTR] .
                                            # It is possible to define more then one contrast Design ["Group,Treatment1,Control1|Group,Treatment2,Control2|..."]
            --SPLIT_BY_CONTRAST             # Only use Samples found in the relevant contrast for Clustering and Enrichment Analysis
            --modelMatrixType:              # How the DeSeq model matrix of the GLM formula is formed [standard or expanded] ,The Default is standard
            --GENES_PLOT:                   # Genes Id To Plot count Data [separated by ','] 
            --X_AXIS:                       # The Filed In the Sample Data To Use as X Axis
            --GROUP:                        # The Filed In the Sample Data To Group By [can be two fields separated by ',']
            --SPLIT_BY:                     # The Filed In the Sample Data To Split the Analysis By.
            --FUNcluster:                   # A clustering function including [kmeans,pam,clara,fanny,hclust,agnes,diana,click]. The default is hclust
                                            # If the 'use_click' option is used the '--FUNcluster' option is set to 'click' 
            --hc_metric:                    # Hierarchical clustering metric to be used for calculating dissimilarities between observations. The default is pearson
            --hc_method:                    # Hierarchical clustering agglomeration method to be used. The default is ward.D2
            --k.max:                        # The maximum number of clusters to consider, must be at least two. The default is 20
            --nboot:                        # Number of Monte Carlo (bootstrap) samples for determining the number of clusters [Not For Mclust]. The default is 10 
            --stand:                        # The Data will be Standardized Before Clustering.
            --Mclust:                       # Use Mclust for determining the number of clusters.
            --CLICK_HOMOGENEITY:            # The HOMOGENEITY [0-1] of clusters using CLICK program (Shamir et al. 2000). The default is 0.5 
            --PCA_COLOR:                    # The Filed In the Sample Data To Determine Color In The PCA Plot
            --PCA_SHAPE:                    # The Filed In the Sample Data To Determine Shape In The PCA Plot
            --PCA_SIZE:                     # The Filed In the Sample Data To Determine Size In The PCA Plot. The default is Library Size
            --Enriched_terms_overlap:       # Test for genes overlap in enriched terms
            --USE_INPUT_GENES_AS_BACKGROUND # Use The input Genes as the Background for Enrichment Analysis
            --only_clustering               # Don't Perform Differential Analysis!!!
            --significant_genes             # Use these genes as the set of significant genes [a comma separated list]
            --collapseReplicates            # Will collapse technical replicates using a Sample Data field indicating which samples are technical replicates
"""
import os
import sys
import re
from neatseq_flow.PLC_step import Step,AssertionExcept


__author__ = "Liron Levin"
__version__= "1.2.0"
# Collect and mereg/append results from all base directories

class Step_Proteomics(Step):

    def step_specific_init(self):
        """ Called on intiation
            Good place for parameter testing.
            Wrong place for sample data testing
        """
        self.shell = "bash"
        
        import inspect
        self.module_location=os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))
        if "--COUNT_DATA_FILE" in self.params["redir_params"] :
            raise AssertionExcept("Please do not define  --COUNT_DATA_FILE")
        if "--SAMPLES" in self.params["redir_params"] :
            raise AssertionExcept("Please do not define  --SAMPLES")
        if "--outDir" in self.params["redir_params"] :
            raise AssertionExcept("Please do not define  --outDir")
    def step_sample_initiation(self):
        """ A place to do initiation stages following setting of sample_data
        """
        self.RSEM_FILES=[]
        self.RSEM_SAMPLES=[]
        
        self.HTSeq_FILES=[]
        self.HTSeq_SAMPLES=[]
        
        for sample in self.sample_data["samples"]:
            if 'RSEM' in list(self.sample_data[sample].keys()):
                self.RSEM_FILES.append(self.sample_data[sample]['RSEM']+'.genes.results' )
                self.RSEM_SAMPLES.append(sample)
            elif 'genes.counts' in list(self.sample_data[sample].keys()):
                self.RSEM_FILES.append(self.sample_data[sample]['genes.counts'])
                self.RSEM_SAMPLES.append(sample)
            elif 'genes.results' in list(self.sample_data[sample].keys()):
                self.RSEM_FILES.append(self.sample_data[sample]['genes.results'])
                self.RSEM_SAMPLES.append(sample)
            if 'HTSeq.counts' in list(self.sample_data[sample].keys()):
                self.HTSeq_FILES.append(self.sample_data[sample]['HTSeq.counts'])
                self.HTSeq_SAMPLES.append(sample)
        
        
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
        
        if self.params["script_path"]==None: 
            if "Proteomics_module.R" in os.listdir(self.module_location): 
                self.params["script_path"]= "Rscript  %s "  % os.path.join(self.module_location,"Proteomics_module.R")
            else:
                raise AssertionExcept("The file %s is not found in the DeSeq2 module directory" % "Proteomics_module.R" )

        if ("use_click" in list(self.params.keys())) and ("--CLICK_PATH" not in self.params["redir_params"]):
            if "click.exe" in os.listdir(self.module_location):
                self.params["redir_params"]["--CLICK_PATH"] = os.path.join(self.module_location,"click.exe")
                self.params["redir_params"]["--FUNcluster"] = 'click'
            else:
                raise AssertionExcept("The file %s is not found in the DeSeq2 module directory \n please use the '--CLICK_PATH' option to indicate the full file location or delete the 'use_click' option" % "click.exe" )
        
        if "--Rmarkdown" not in self.params["redir_params"]:
            if "DeSeq2_module.Rmd" in os.listdir(self.module_location):
                self.params["redir_params"]["--Rmarkdown"] = os.path.join(self.module_location,"Proteomics_module.Rmd")
        
        if self.params["script_path"]!=None:    
            # Get constant part of script:
            self.script += self.get_script_const()
            if "--COUNT_SOURCE" not in self.params["redir_params"] :
                if len(self.RSEM_FILES)>0:
                    self.script += "--COUNT_DATA_FILE %s \\\n\t" %  str(self.RSEM_FILES).replace('[','').replace(']','').replace('"','').replace(', ',',').replace("'",'') 
                    self.script += "--SAMPLES %s \\\n\t" %  str(self.RSEM_SAMPLES).replace('[','').replace(']','').replace('"','').replace(', ',',').replace("'",'') 
                    self.script += "--COUNT_SOURCE %s \\\n\t" %  'RSEM' 
                elif len(self.HTSeq_FILES)>0: 
                    self.script += "--COUNT_DATA_FILE %s \\\n\t" %  str(self.HTSeq_FILES).replace('[','').replace(']','').replace('"','').replace(', ',',').replace("'",'') 
                    self.script += "--SAMPLES %s \\\n\t" %  str(self.HTSeq_SAMPLES).replace('[','').replace(']','').replace('"','').replace(', ',',').replace("'",'') 
                    self.script += "--COUNT_SOURCE %s \\\n\t" %  'HTSEQ' 
                elif "results" in list(self.sample_data["project_data"].keys()):
                    self.script += "--COUNT_DATA_FILE %s \\\n\t" %  self.sample_data["project_data"]["results"]
                    self.script += "--COUNT_SOURCE %s \\\n\t" %  'Matrix' 
                else:
                    raise AssertionExcept("Could not fined count data [From RSEM, HTSeq or from collect results]")
            else:
                if (len(self.RSEM_FILES)>0) and (self.params["redir_params"]["--COUNT_SOURCE"] == 'RSEM'):
                    self.script += "--COUNT_DATA_FILE %s \\\n\t" %  str(self.RSEM_FILES).replace('[','').replace(']','').replace('"','').replace(', ',',').replace("'",'') 
                    self.script += "--SAMPLES %s \\\n\t" %  str(self.RSEM_SAMPLES).replace('[','').replace(']','').replace('"','').replace(', ',',').replace("'",'') 
                    self.script += "--COUNT_SOURCE %s \\\n\t" %  'RSEM' 
                elif (len(self.HTSeq_FILES)>0)  and (self.params["redir_params"]["--COUNT_SOURCE"] == 'HTSEQ'): 
                    self.script += "--COUNT_DATA_FILE %s \\\n\t" %  str(self.HTSeq_FILES).replace('[','').replace(']','').replace('"','').replace(', ',',').replace("'",'') 
                    self.script += "--SAMPLES %s \\\n\t" %  str(self.HTSeq_SAMPLES).replace('[','').replace(']','').replace('"','').replace(', ',',').replace("'",'') 
                    self.script += "--COUNT_SOURCE %s \\\n\t" %  'HTSEQ' 
                elif ("results" in list(self.sample_data["project_data"].keys())) and (self.params["redir_params"]["--COUNT_SOURCE"] == 'Matrix'):
                    self.script += "--COUNT_DATA_FILE %s \\\n\t" %  self.sample_data["project_data"]["results"]
                    self.script += "--COUNT_SOURCE %s \\\n\t" %  'Matrix' 
                else:
                    raise AssertionExcept("Could not fined %s count data " % self.params["redir_params"]["--COUNT_SOURCE"])
            
            
            if ('trino.rep' in list(self.sample_data["project_data"].keys())) and ("--Trinotate" not in list(self.params["redir_params"].keys())) :
                self.script += "--Trinotate %s \\\n\t" % self.sample_data["project_data"]['trino.rep']
            
            
            self.script += "--outDir %s \n\n" % use_dir
            
            
            
        # Wrapping up function. Leave these lines at the end of every iteration:
        self.local_finish(use_dir,sample_dir)       # Sees to copying local files to final destination (and other stuff)
                  
        
        self.create_low_level_script()
