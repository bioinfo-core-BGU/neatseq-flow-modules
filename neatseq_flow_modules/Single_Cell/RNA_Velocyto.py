import scvelo as scv
from matplotlib import pyplot as plt
adata = scv.read("/gpfs0/biores/users/gadlab/Daniel/sNuc-seq/all_lines_only_old103/01.RUN/data/Generic/CellRanger/CCHS102/CCHS102/velocyto/Test.h5ad")
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
scv.tl.velocity(adata)
scv.tl.velocity_graph(adata)

 x=scv.pl.velocity_embedding_stream(adata, basis="umap",
                                    color = "seurat_clusters",
                                    save  = ,
                                    dpi   = 600,
                                    show=False)

scv.tl.rank_velocity_genes(adata, groupby='seurat_clusters', min_corr=.3)
df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
df.head()
scv.pl.scatter(adata, df['0'][:5], ylabel='0',color="seurat_clusters")
scv.pl.scatter(adata, df['0'][:5], ylabel='0',color="velocity")
scv.tl.velocity_pseudotime(adata)
scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot')

adata.uns['neighbors']['distances'] = adata.obsp['distances']
adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']

scv.tl.paga(adata, groups='seurat_clusters')
df = scv.get_df(adata, 'paga/transitions_confidence', precision=2).T


scv.pl.paga(adata,
            basis='umap',
            size=50, alpha=.1, dpi=600,show=False,
            save="/gpfs0/biores/users/gadlab/Daniel/sNuc-seq/all_lines_only_old103/01.RUN/data/Generic/CellRanger/CCHS102/CCHS102/velocyto/plot3.png",
            min_edge_width=2, 
            node_size_scale=1.5)