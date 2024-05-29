import scrublet as scr
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from matplotlib import pyplot as plt
import pandas as pd
if __name__ == '__main__':
    #getting arguments from the user 
    import argparse
    parser = argparse.ArgumentParser(description='Scrublet')
    parser.add_argument('--RDS',dest='RDS',metavar="CHAR",type=str,default=None,
                        help='''Input RDS file
                                ''')
    parser.add_argument('--min_counts',dest='min_counts',metavar="CHAR",type=int,default=3,
                        help='''Used for gene filtering prior to PCA. Genes expressed at fewer than 
                               `min_counts` in fewer than `min_cells` (see below) are excluded. default=3
                                ''')
    parser.add_argument('--min_cells',dest='min_cells',metavar="CHAR",type=int,default=3,
                        help='''Used for gene filtering prior to PCA. Genes expressed at fewer than 
                               `min_counts` (see above) in fewer than `min_cells` are excluded. default=3
                                ''')
    parser.add_argument('--expected_doublet_rate',dest='expected_doublet_rate',metavar="CHAR",type= float,default=0.1,
                        help='''The estimated doublet rate for the experiment. default = 0.1 (10%)
                                ''')
    parser.add_argument('--threshold',dest='threshold',metavar="CHAR",type=float,default=0.0,
                        help='''Manual threshold for calling doublets. default = Automatically
                                ''')
    parser.add_argument('--n_neighbors',dest='n_neighbors',metavar="CHAR",type=int,default=30,
                        help='''Number of neighbors to use. 
                                ''')
    parser.add_argument('--min_gene_variability_pctl',dest='min_gene_variability_pctl',metavar="CHAR",type=int,default=85,
                        help='''Used for gene filtering prior to PCA. Keep the most highly variable genes
                                (in the top min_gene_variability_pctl percentile), as measured by 
                                 the v-statistic [Klein et al., Cell 2015]’. default=85
                                ''')
    parser.add_argument('--n_prin_comps',dest='n_prin_comps',metavar="CHAR",type=int,default=30,
                        help='''Number of principal components used to embed the transcriptomes prior
                                to k-nearest-neighbor graph construction. default=30
                                ''')
    parser.add_argument('--output_dir',dest='output_dir',metavar="CHAR",type=str,default=None,
                        help='''Output Directory .
                                ''')
                                
    args          = parser.parse_args()


Seurat    = importr('Seurat')
Seurat    = importr('SeuratObject')
readRDS   = robjects.r['readRDS']
as_matrix = robjects.r['as.matrix']
t_default = robjects.r['t.default']

df  = readRDS(args.RDS)
RNA = Seurat.GetAssayData(df, assay="RNA",slot="counts")
counts_matrix = t_default(as_matrix(RNA))

scrub = scr.Scrublet(counts_matrix,expected_doublet_rate=args.expected_doublet_rate)


doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts = args.min_counts, 
                                                          min_cells  = args.min_cells, 
                                                          min_gene_variability_pctl = args.min_gene_variability_pctl,
                                                          n_prin_comps = args.n_prin_comps)

if args.threshold>0:
    predicted_doublets = scrub.call_doublets(threshold=args.threshold)

Results = pd.DataFrame({'Predicted_Doublets':predicted_doublets , 'Doublet_Scores':doublet_scores},index=counts_matrix.rownames)

Results.loc[Results.Predicted_Doublets==False,"DF.classification_homotypic"] = "Singlet"
Results.loc[Results.Predicted_Doublets==True ,"DF.classification_homotypic"] = "Doublet"



Results.to_csv(args.output_dir+"/Classification.csv")

try:
    scrub.plot_histogram()
    plt.savefig(args.output_dir+"/Histogram.pdf", format="pdf")
except:
    scrub.threshold_ = 1
    scrub.plot_histogram()
    plt.savefig(args.output_dir+"/Histogram.pdf", format="pdf")
