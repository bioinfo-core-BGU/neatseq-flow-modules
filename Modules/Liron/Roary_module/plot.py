#!/usr/bin/env python

__author__ = "Marco Galardini_modified"
__version__ = '0.1.0'

def get_options():
    import argparse

    #  create the top-level parser
    description = "Create plots from Roary outputs"
    parser = argparse.ArgumentParser(description = description,
                                     prog = 'plot.py')

    parser.add_argument('tree', action='store',
                        help='Newick Tree file', default='accessory_binary_genes.fa.newick')
    parser.add_argument('spreadsheet', action='store',
                        help='Roary gene presence/absence spreadsheet', default='gene_presence_absence.csv')
    parser.add_argument('-O', action='store',
                        help='Output location', default='gene_presence_absence.csv')
    parser.add_argument('-C', action='store',
                        help='Clustering method', default='ward')
    parser.add_argument('-T', '--tag' , action='store',
                        help='Mark genes in the matrix containing this key')
    parser.add_argument('--labels', action='store_true',
                        default=False,
                        help='Add node labels to the tree (up to 10 chars)')
    parser.add_argument('--format',
                        choices=('png',
                                 'tiff',
                                 'pdf',
                                 'svg'),
                        default='png',
                        help='Output format [Default: png]')
    parser.add_argument('-L', '--low' , action='store',default=0.01,
                        help='Low gene frequency cutoff  Default: 0.01 ')
    parser.add_argument('--version', action='version',
                         version='%(prog)s '+__version__)

    return parser.parse_args()

if __name__ == "__main__":
    options = get_options()
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import seaborn as sns

    sns.set_style('white')

    import os
    import pandas as pd
    import numpy as np
    from Bio import Phylo

    t = Phylo.read(options.tree, 'newick')

    print ' Max distance to create better plots'
    mdist = max([t.distance(t.root, x) for x in t.get_terminals()])

    print ' Load roary'
    roary = pd.read_table(options.spreadsheet,
                         sep=',',
                         low_memory=False)
    print ' Set index (group name)'
    roary.set_index('Gene', inplace=True)
    print '   Find tagged genes'
    if options.tag is not None:
        if "Inference" in roary.columns:
            roary_info=roary[["Inference"]].copy()
            roary.drop(["Inference"], axis=1, inplace=True)
            VF=roary_info[roary_info.applymap(lambda x: options.tag.upper() in x.upper()).values].index
        else:
            roary_info=roary[["Annotation"]].copy()
            VF=roary_info[roary_info.applymap(lambda x: options.tag.upper() in x.upper()).values].index
    else:
        VF=[]
        
    
    roary.loc[:,roary.columns[13:]]=roary.loc[:,roary.columns[13:]].applymap(lambda x: 0 if pd.isnull(x) else 1)
    roary.to_csv(os.path.join(options.O,'gene_presence_absence.%s' % 'tab'), sep='\t',float_format="%g",index=True)
    print ' Drop the other info columns'
    roary.drop(list(roary.columns[:13]), axis=1, inplace=True)

    # print ' Transform it in a presence/absence matrix (1/0)'
    # roary.replace('.{2,100}', 1, regex=True, inplace=True)
    # roary.replace(np.nan, 0, regex=True, inplace=True)

    print ' Sort the matrix by the sum of strains presence'
    idx = roary.sum(axis=1).sort_values(ascending=False).index
    roary_sorted = roary.ix[idx]
    

    print ' Pangenome frequency plot'
    plt.figure(figsize=(7, 5))

    plt.hist(roary.sum(axis=1), roary.shape[1],
             histtype="stepfilled", alpha=.7)

    plt.xlabel('No. of genomes')
    plt.ylabel('No. of genes')

    sns.despine(left=True,
                bottom=True)
    plt.savefig(os.path.join(options.O ,'pangenome_frequency.%s' %options.format), dpi=600)
    plt.clf()

    print ' Sort the matrix according to tip labels in the tree'
    roary_sorted = roary_sorted[[x.name for x in t.get_terminals()]]

    print ' Plot presence/absence matrix against the tree'
    from scipy.cluster.hierarchy import dendrogram, linkage, to_tree
    import seaborn as sns; sns.set()
    roary_sorted=roary_sorted.loc[roary_sorted.mean(axis=1)>options.low,]

    def change_data_val(data,col,old_v,new_v):
        for i in col:
            if i in data.columns:
                data.loc[data[i]==old_v,[i]]=new_v
        return data

    if roary.T.shape[0]<150:
        fontsize=5
    else:
        fontsize=(40./float(roary_sorted.T.shape[0]))*15

    with sns.axes_style('whitegrid'):
        fig = plt.figure(figsize=(17, 20))

        ax1=plt.subplot2grid((50,50), (13, 15), colspan=30,rowspan=35)
        ax2=plt.subplot2grid((50,50),(13, 0), colspan=10,rowspan=35, axisbg='white')
        #ax3=plt.subplot2grid((50,50),(0, 15), colspan=30,rowspan=10, axisbg='white')
        
        fig.subplots_adjust(wspace=0, hspace=0)

        ax1.set_title('Roary matrix\n(%d gene clusters)'%roary_sorted.shape[0])
        plt.rcParams["lines.linewidth"]=1
        
        Z = linkage(roary_sorted.loc[:,roary_sorted.columns].T.values,options.C)

        g=dendrogram(
            Z,
            leaf_rotation=0.,  # rotates the x axis labels
            leaf_font_size=fontsize,  # font size for the x axis labels
            labels=map(lambda x:x[0:15],roary_sorted.columns) ,
            orientation="right",
            ax=ax2,
         )
        m=g["leaves"]
        m.reverse()
        
        
        # ZZ = linkage(roary_sorted.loc[:,roary_sorted.columns].values,options.C)

        # gg=dendrogram(
            # ZZ,
            # leaf_rotation=90,  # rotates the x axis labels
            # leaf_font_size=fontsize, # font size for the x axis labels
            # labels=map(lambda x:x[0:15],roary_sorted.index) ,
            # orientation="top",
            # ax=ax3,
         # )
        # mm=gg["leaves"]
        # mm.reverse()
        # ax3.invert_xaxis()
        
        
        #roary_sorted_new=roary_sorted.loc[roary_sorted.index[mm],roary_sorted.columns[m]].copy()
        roary_sorted_new=roary_sorted.loc[:,roary_sorted.columns[m]].copy()
        roary_sorted_new.loc[VF,].to_csv(os.path.join(options.O,'pangenome_matrix.%s' % 'tab'), sep='\t',float_format="%g",index=True)
        roary_sorted_new=change_data_val(roary_sorted_new.T,VF,1,-1)
        a=ax1.imshow(roary_sorted_new, cmap=plt.cm.bwr_r,
               vmin=-1, vmax=1,
               aspect='auto',
               interpolation='none',
                )
        
        ax1.set_yticks(range(len(roary_sorted_new.index)),minor=False)
        ax2.set_xticks([])
        #ax3.set_yticks([])
        
        ax1.set_yticklabels(map(lambda x:x[0:15],roary_sorted_new.index),fontsize=fontsize)
        ax1.grid('off')
        ax2.grid('off')
        #ax3.grid('off')
        plt.savefig(os.path.join(options.O,'pangenome_matrix.%s' % options.format), dpi=600)
        plt.clf()
        
        
        

    print ' write new tree file'
    def getNewick(node, newick, parentdist, leaf_names):
        if node.is_leaf():
            return "%s:%.2f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
        else:
            if len(newick) > 0:
                newick = "):%.2f%s" % (parentdist - node.dist, newick)
            else:
                newick = ");"
            newick = getNewick(node.get_left(), newick, node.dist, leaf_names)
            newick = getNewick(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
            newick = "(%s" % (newick)
            return newick

    tree = to_tree(Z,False)
    h=open(os.path.join(options.O,'pangenome_matrix.%s' % "newick"),'w')
    h.write( getNewick(tree, "", tree.dist, map(lambda x:x[0:15],roary_sorted.columns) ))
    h.close()
    
        
    print ' Plot the pangenome pie chart'
    plt.figure(figsize=(10, 10))

    core     = roary[(roary.sum(axis=1) >= roary.shape[1]*0.99) & (roary.sum(axis=1) <= roary.shape[1]     )].shape[0]
    softcore = roary[(roary.sum(axis=1) >= roary.shape[1]*0.95) & (roary.sum(axis=1) <  roary.shape[1]*0.99)].shape[0]
    shell    = roary[(roary.sum(axis=1) >= roary.shape[1]*0.15) & (roary.sum(axis=1) <  roary.shape[1]*0.95)].shape[0]
    cloud    = roary[roary.sum(axis=1)  < roary.shape[1]*0.15].shape[0]

    total = roary.shape[0]
    
    def my_autopct(pct):
        val=int(round(pct*total/100.0))
        return '{v:d}'.format(v=val)

    a=plt.pie([core, softcore, shell, cloud],
          labels=['core\n(%d <= strains <= %d)'%(roary.shape[1]*.99,roary.shape[1]),
                  'soft-core\n(%d <= strains < %d)'%(roary.shape[1]*.95,roary.shape[1]*.99),
                  'shell\n(%d <= strains < %d)'%(roary.shape[1]*.15,roary.shape[1]*.95),
                  'cloud\n(strains < %d)'%(roary.shape[1]*.15)],
          explode=[0.1, 0.05, 0.02, 0], radius=0.9,
          colors=[(0, 0, 1, float(x)/total) for x in (core, softcore, shell, cloud)],
          autopct=my_autopct)
    plt.savefig(os.path.join(options.O,'pangenome_pie.%s' % options.format), dpi=600)
    plt.clf()
