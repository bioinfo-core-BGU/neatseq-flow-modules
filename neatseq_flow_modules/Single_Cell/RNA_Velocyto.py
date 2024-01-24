import scvelo as scv
from matplotlib import pyplot as plt
if __name__ == '__main__':
    #getting arguments from the user 
    import argparse
    parser = argparse.ArgumentParser(description='scvelo_py')
    parser.add_argument('--h5ad',dest='h5ad',metavar="CHAR",type=str,default=None,
                        help='''Input h5ad file
                                ''')
    parser.add_argument('--min_shared_counts',dest='min_shared_counts',metavar="CHAR",type=int,default=0,
                        help='''Minimum number of cells expressed required to pass filtering (spliced).
                                ''')
    parser.add_argument('--n_top_genes',dest='n_top_genes',metavar="CHAR",type=int,default=5000,
                        help='''Number of genes to keep.
                                ''')
    parser.add_argument('--n_pcs',dest='n_pcs',metavar="CHAR",type=int,default=30,
                        help='''Number of principal components to use.
                        If not specified, the full space is used of a pre-computed PCA, or 30 components are used when PCA is computed internally.
                                ''')
    parser.add_argument('--n_neighbors',dest='n_neighbors',metavar="CHAR",type=int,default=30,
                        help='''Number of neighbors to use.
                                ''')
    parser.add_argument('--embedding',dest='embedding',metavar="CHAR",type=str,default="umap",
                        help='''Key for embedding. If not specified, use ‘umap’.
                                ''')
    parser.add_argument('--color',dest='color',metavar="CHAR",type=str,default="seurat_clusters",
                        help='''Key for annotations of observations/cells or variables/genes. If not specified, use ‘seurat_clusters’.
                                ''')
    parser.add_argument('--output_dir',dest='output_dir',metavar="CHAR",type=str,default=None,
                        help='''Key for annotations of observations/cells or variables/genes. If not specified, use ‘seurat_clusters’.
                                ''')
                                
    args          = parser.parse_args()


adata = scv.read(args.h5ad)
scv.pp.filter_and_normalize(adata, 
                            min_shared_counts = args.min_shared_counts,
                            n_top_genes       = args.n_top_genes)
                            
scv.pp.moments(adata,
               n_pcs       = args.n_pcs,
               n_neighbors = args.n_neighbors)
               
scv.tl.velocity(adata)
scv.tl.velocity_graph(adata)

x=scv.pl.velocity_embedding_stream(adata, 
                                   basis = args.embedding,
                                   color = args.color,
                                   save  = args.output_dir + "Velocity_Embedding_Stream.png",
                                   dpi   = 600,
                                   show=False)

# scv.tl.rank_velocity_genes(adata, groupby='seurat_clusters', min_corr=.3)
# df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
# df.head()
# scv.pl.scatter(adata, df['0'][:5], ylabel='0',color="seurat_clusters")
# scv.pl.scatter(adata, df['0'][:5], ylabel='0',color="velocity")
# scv.tl.velocity_pseudotime(adata)
# scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot')

# adata.uns['neighbors']['distances'] = adata.obsp['distances']
# adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']

# scv.tl.paga(adata, groups='seurat_clusters')
# df = scv.get_df(adata, 'paga/transitions_confidence', precision=2).T


# scv.pl.paga(adata,
            # basis='umap',
            # size=50, alpha=.1, dpi=600,show=False,
            # save="/gpfs0/biores/users/gadlab/Daniel/sNuc-seq/all_lines_only_old103/01.RUN/data/Generic/CellRanger/CCHS102/CCHS102/velocyto/plot3.png",
            # min_edge_width=2, 
            # node_size_scale=1.5)