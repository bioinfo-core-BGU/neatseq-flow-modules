#!/usr/bin/env python

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description = "Roary_matrix_plot",
                                     prog = 'Roary_matrix_plot.py')
    parser.add_argument('-P','--presence_absence', action='store',
                        help='Roary gene presence/absence csv file', default='gene_presence_absence.csv')
    parser.add_argument('-O', action='store',
                        help='Output location', default='gene_presence_absence.csv')
    parser.add_argument('-C', action='store',
                        help='Clustering method', default='ward')
    parser.add_argument('-T', '--tag' , action='store',
                        help='Mark genes in the matrix containing this key')
    parser.add_argument('--format', choices=('png','tiff','pdf','svg'),default='pdf',
                        help='Output format [Default: pdf]')
    parser.add_argument('-L', '--low' , action='store',default=0.01,
                        help='Low gene frequency cutoff  Default: 0.01 ')
    options=parser.parse_args()
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import seaborn as sns
    import os
    import pandas as pd
    import numpy as np
    sns.set_style('white')

    print(' Load roary')
    roary = pd.read_table(options.presence_absence,
                         sep=',',
                         low_memory=False)
    print(' Set index (group name)')
    roary.set_index('Gene', inplace=True)
    print('   Find tagged genes')
    if options.tag is not None:
        if "Inference" in roary.columns:
            roary_info=roary[["Inference"]].copy()
            roary.drop(["Inference"], axis=1, inplace=True)
            print(roary_info)
            VF=roary_info[roary_info.applymap(lambda x: options.tag.upper() in str(x).upper()).values].index
        else:
            roary_info=roary[["Annotation"]].copy()
            VF=roary_info[roary_info.applymap(lambda x: options.tag.upper() in str(x).upper()).values].index
    else:
        VF=[]
        
    
    roary.loc[:,roary.columns[13:]]=roary.loc[:,roary.columns[13:]].applymap(lambda x: 0 if pd.isnull(x) else 1)
    roary.to_csv(os.path.join(options.O,'gene_presence_absence.%s' % 'tab'), sep='\t',float_format="%g",index=True)
    print(' Drop the other info columns')
    roary.drop(list(roary.columns[:13]), axis=1, inplace=True)
    print(' Sort the matrix by the sum of strains presence')
    idx = roary.sum(axis=1).sort_values(ascending=False).index
    roary_sorted = roary.ix[idx]
    

    from scipy.cluster.hierarchy import dendrogram, linkage, to_tree
    import seaborn as sns; sns.set()
    roary_sorted=roary_sorted.loc[roary_sorted.mean(axis=1)>options.low,]

    def change_data_val(data,col,old_v,new_v):
        for i in col:
            if i in data.columns:
                data.loc[data[i]==old_v,[i]]=new_v
        return data
        
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

######## Generate the pangenome matrix plot and the Accessory genes tree
    temp_roary_sorted=roary_sorted.copy()        
    if roary_sorted.T.shape[0]<150:
        fontsize=2
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
            labels=[x[0:15] for x in roary_sorted.columns] ,
            orientation="left",
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
        roary_sorted_new=change_data_val(roary_sorted_new.T,VF,1,-1)
        a=ax1.imshow(roary_sorted_new, cmap=plt.cm.bwr_r,
               vmin=-1, vmax=1,
               aspect='auto',
               interpolation='none',
                )
        
        ax1.set_yticks(list(range(len(roary_sorted_new.index))),minor=False)
        ax2.set_xticks([])
        #ax3.set_yticks([])
        
        ax1.set_yticklabels([x[0:15] for x in roary_sorted_new.index],fontsize=fontsize)
        ax1.grid('off')
        ax2.grid('off')
        #ax3.grid('off')
        plt.savefig(os.path.join(options.O,'pangenome_matrix.%s' % options.format), dpi=600)
        plt.clf()
        
        print(' Write Accessory genes tree file')


        tree = to_tree(Z,False)
        h=open(os.path.join(options.O,'Accessory.%s' % "newick"),'w')
        h.write( getNewick(tree, "", tree.dist, [x for x in roary_sorted.columns] ))
        h.close()
    



######## Generate the virulence/resistance matrix plot and the virulence/resistance tree
    if len(VF)>0:
        
        roary_sorted=temp_roary_sorted.loc[VF,].copy()
        
        if roary_sorted.T.shape[0]<150:
            fontsize=3
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
                labels=[x[0:15] for x in roary_sorted.columns] ,
                orientation="left",
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
            roary_sorted_new=change_data_val(roary_sorted_new.T,VF,1,-1)
            a=ax1.imshow(roary_sorted_new, cmap=plt.cm.bwr_r,
                   vmin=-1, vmax=1,
                   aspect='auto',
                   interpolation='none',
                    )
            
            ax1.set_yticks(list(range(len(roary_sorted_new.index))),minor=False)
            ax2.set_xticks([])
            #ax3.set_yticks([])
            
            ax1.set_yticklabels([x[0:15] for x in roary_sorted_new.index],fontsize=fontsize)
            ax1.grid('off')
            ax2.grid('off')
            #ax3.grid('off')
            plt.savefig(os.path.join(options.O,'virulence_resistance.%s' % options.format), dpi=600)
            plt.clf()
            
            tree = to_tree(Z,False)
            h=open(os.path.join(options.O,'virulence_resistance.%s' % "newick"),'w')
            h.write( getNewick(tree, "", tree.dist, [x for x in roary_sorted.columns] ))
            h.close()

            roary_sorted=roary_sorted.T.copy()
            roary_sorted.index.names=['Samples']
            roary_sorted.to_csv(os.path.join(options.O,'virulence_resistance.%s' % 'tab'), sep='\t',float_format="%g",index=True)
            print(' Write virulence/resistance tree file')
