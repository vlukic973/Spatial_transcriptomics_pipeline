print('Importing packages')

import os
import numpy as np
import pandas as pd
import itertools
from itertools import permutations
import scanpy as sc
#import squidpy as sq
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import anndata as ad
from pathlib import Path
import scipy
from matplotlib.pyplot import rc_context
import os
import csv
sc.logging.print_header()
import random
import argparse
from scipy.stats import pearsonr, spearmanr
from scipy import stats
import tracemalloc
import time
#sc.logging.print_header()

# Record start time
start_time = time.time()

# Create the argument parser
parser = argparse.ArgumentParser(description='Find genes expressed in cancer cells connected to normal cells, and vice versa, when splitting up a large list of gene combinations')

# Add the --path argument
parser.add_argument('--path', type=str, help='Path to files')
#parser.add_argument('--h5ad_file', type=str, help='Name of file')
parser.add_argument('--cell_subtype', type=str, help='Normal cell subtype')
parser.add_argument('--direction', type=str, help='Cancer cell is ligand or receptor?')
parser.add_argument('--split_part', type=int, required=True, help='Split part (1, 2, 3, ...)')
parser.add_argument('--split_size', type=int, required=True, help='Size of split e.g. how many gene-pairs to consider in a chunk')
parser.add_argument('--start_index', type=int, required=True, help='Which index to start on across the gene combinations')
parser.add_argument('--end_index', type=int, required=True, help='Which index to end on across the gene combinations')
parser.add_argument('--ct_level', type=str, required=True, help='Cell type annotation')
parser.add_argument('--version', type=str, required=True, help='Run version')
parser.add_argument('--corr_coef', type=str, choices=['pearson', 'spearman'], default='spearman', help='Correlation coefficient method (pearson or spearman)')
parser.add_argument('--sample',type=str,help="Which sample")
parser.add_argument('--fake', action='store_true', help='Consider all genes')
parser.add_argument('--ligrec_known', action='store_true', help='Consider known lig-rec pairs')

args = parser.parse_args()

if args.fake:
    print("Considering non DE-genes")
    fake=True
if args.ligrec_known:
    print("Considering known ligrec_pairs")
    ligrec_known=True
else:
    print("Considering DE-genes")
    fake=False

# Access the value of the --path argument
path = args.path
#h5ad_file = args.h5ad_file 
cell_subtype = args.cell_subtype
direction = args.direction
split_part=args.split_part
split_size=args.split_size
start_index=args.start_index
end_index=args.end_index
ct_level=args.ct_level
version=args.version
sample=args.sample
fake=args.fake
ligrec_known=args.ligrec_known
corr_coef=args.corr_coef

print('Printing command line arguments')

print('path:\n')
print(path)
print('cell_subtype:\n')
print(cell_subtype)
print('direction:\n')
print(direction)
print('split part:\n')
print(split_part)
print('split size:\n')
print(split_size)
print('start_index:\n')
print(start_index)
print('end index:\n')
print(end_index)
print('ct_level:\n')
print(ct_level)
print('version:\n')
print(version)
print('fake:\n')
print(fake)
print('ligrec_known:\n')
print(ligrec_known)
print('corr_coef:\n')
print(corr_coef)
print('sample:\n')
print(sample)
#fake=args.fake

if (fake):
    fake_flag='all_genes'
if (ligrec_known):
    fake_flag='known_ligrec_genes'
else:
    fake_flag='DE_genes'

print('fake_flag')
print(fake_flag)

if corr_coef == 'pearson':
    # Use Pearson correlation coefficient
    correlation_function = pearsonr
elif corr_coef == 'spearman':
    # Use Spearman correlation coefficient
    correlation_function = spearmanr


def lig_rec_plot_normal_cancer2(adata,cell_subtype,impute=True,direction='NC',lig_rec_pair=('PDGFB', 'LRP1'),fake=False,corr_coef=spearmanr):
    if (impute):
        #print('Using imputed values')
        adata.X=adata.obsm['X_normalized_scVI']
        assert(np.sum(np.asarray(adata.X))==np.sum(np.asarray(adata.obsm['X_normalized_scVI'])))
    else:
        print('Not using imputed values')
    print('Considering '+Normal_cells[cell_subtype])
    print('Cancer_cell')
    print(Cancer_cell)
    adata_cancer_subtype=adata[(adata.obs[ct_level]==(Cancer_cell)) | (adata.obs[ct_level]==(Normal_cells[cell_subtype]))]
    assert(len(adata_cancer_subtype)==adata_cancer_subtype.obsp['connectivities'].shape[0])
    assert(len(adata_cancer_subtype)==adata_cancer_subtype.obsp['connectivities'].shape[1])
    #print("Converting connectivities to array")
    adjacency_matrix=adata_cancer_subtype.obsp['connectivities'].toarray()
    #print("Converting to data frame")
    adjacency_matrix_df=pd.DataFrame(adjacency_matrix,index=adata_cancer_subtype.obs.index, columns=adata_cancer_subtype.obs.index)
    #print('np.unique(adata_cancer_subtype.obs[ct_level])')
    #print(np.unique(adata_cancer_subtype.obs[ct_level]))
    #print('adata_cancer_subtype.obs[(adata_cancer_subtype.obs[ct_level] == Cancer_cell)]')
    #print(adata_cancer_subtype.obs[(adata_cancer_subtype.obs[ct_level] == Cancer_cell)])
    index_list1 = adata_cancer_subtype.obs[(adata_cancer_subtype.obs[ct_level] == Cancer_cell)].index
    #print('original index_list1')
    #print(index_list1)
    index_list2 = adata_cancer_subtype.obs[(adata_cancer_subtype.obs[ct_level]==Normal_cells[cell_subtype])].index
    index_list1 = index_list1.astype(str)
    index_list2 = index_list2.astype(str)
    #print('index_list1')
    #print(index_list1)
    #print('index_list2')
    #print(index_list2)
    adjacency_matrix_rows = adjacency_matrix_df.loc[index_list2, index_list1]
    row_labels = adjacency_matrix_rows.index[np.where(adjacency_matrix_rows == 1.0)[0]]
    #print('row_labels')
    #print(row_labels)
    #print('type(row_labels)')
    #print(type(row_labels))
    col_labels = adjacency_matrix_rows.columns[np.where(adjacency_matrix_rows == 1.0)[1]]
    #print('col_labels')
    #print(col_labels)
    row_labels_uc = adjacency_matrix_rows.index[np.where((adjacency_matrix_rows != 1.0).all(axis=1) & (adjacency_matrix_rows.sum(axis=1) == 0))[0]]
    col_labels_uc = adjacency_matrix_rows.columns[np.where((adjacency_matrix_rows != 1.0).all(axis=0) & (adjacency_matrix_rows.sum(axis=0) == 0))[0]]
    #rdm_sample=random.sample(range(len(col_labels_uc)), len(col_labels))
    #col_labels=col_labels.astype('str')
    #col_labels_uc=col_labels_uc.astype('str')
    #print('col_labels after conversion to string')
    #print(col_labels)
    if (direction=='CN'):
        #print('direction is CN')
        #assert(direction=='CN')
        #print('Assuming cancer is ligand and normal is receptor')
        cancer_ligand=lig_rec_pair[0]
        normal_receptor=lig_rec_pair[1]
        idx_ligand=adata.var_names.get_loc(cancer_ligand)
        idx_receptor=adata.var_names.get_loc(normal_receptor)
        #cancer_ligand_list = [adata_cancer_subtype[col_label].X[0][idx_ligand] for col_label in col_labels]
        #normal_receptor_list = [adata_cancer_subtype[row_label].X[0][idx_receptor] for row_label in row_labels]
        cancer_ligand_list = adata_cancer_subtype[col_labels, :].X[:, idx_ligand]
        normal_receptor_list = adata_cancer_subtype[row_labels, :].X[:, idx_receptor]
        #cancer_ligand_list = adata_cancer_subtype[col_labels,:].X[0][idx_ligand]
        #print('Done cancer_ligand_list successfully')
        #normal_receptor_list = adata_cancer_subtype[row_labels,:].X[0][idx_receptor]
        #print('Done normal_receptor_list successfully')
        #print('cancer_ligand_list')
        #print(cancer_ligand_list)
        #exit()
        #print('normal_receptor_list')
        #print(normal_receptor_list)
        #print('len(col_labels_uc),len(row_labels_uc)')
        #print(len(col_labels_uc),len(row_labels_uc))
        #min_length = min(len(col_labels_uc), len(row_labels_uc))
        #min_length = len(col_labels)
        #col_labels_uc_truncated = col_labels_uc[:min_length]
        #row_labels_uc_truncated = row_labels_uc[:min_length]
        #cancer_ligand_list_uc = [adata_cancer_subtype[col_label].X[0][idx_ligand] for col_label in col_labels_uc_truncated]
        #normal_receptor_list_uc = [adata_cancer_subtype[row_label].X[0][idx_receptor] for row_label in row_labels_uc_truncated]
        #cancer_ligand_list_uc = [adata_cancer_subtype[col_label].X[0][idx_ligand] for col_label in col_labels_uc]
        #normal_receptor_list_uc = [adata_cancer_subtype[row_label].X[0][idx_receptor] for row_label in row_labels_uc]
        cancer_ligand_list_uc = adata_cancer_subtype[col_labels_uc,:].X[:,idx_ligand]
        normal_receptor_list_uc = adata_cancer_subtype[row_labels_uc,:].X[:,idx_receptor]
        #cancer_ligand_list_uc_repeats=np.repeat(cancer_ligand_list_uc,len(normal_receptor_list_uc))
        #normal_receptor_list_uc_repeats=np.repeat(normal_receptor_list_uc,len(cancer_ligand_list_uc))
        cancer_ligand_list_uc = np.array(cancer_ligand_list_uc)
        normal_receptor_list_uc = np.array(normal_receptor_list_uc)
        mean_cancer_list=np.mean(cancer_ligand_list)
        mean_normal_list=np.mean(normal_receptor_list)
        mean_cancer_list_uc=np.mean(cancer_ligand_list_uc)
        mean_normal_list_uc=np.mean(normal_receptor_list_uc)
        # Number of combinations to sample
        N = 0.1*len(normal_receptor_list_uc)*len(cancer_ligand_list_uc)
        #N=100
        #print('N')
        #print(N)
        # Sample N random pairs
        sampled_indices1 = np.random.choice(len(cancer_ligand_list_uc), int(N))
        sampled_indices2 = np.random.choice(len(normal_receptor_list_uc), int(N))
        # Extract sampled elements
        #print('sampled_indices1')
        #print(sampled_indices1)
        #print('type(sampled_indices1)')
        #print(type(sampled_indices1))
        #print('cancer_ligand_list_uc')
        #print(type(cancer_ligand_list_uc))
        sampled_cancer_ligand_list_uc = cancer_ligand_list_uc[sampled_indices1]
        sampled_normal_receptor_list_uc = normal_receptor_list_uc[sampled_indices2]
        corr_unconnected=round(corr_coef(sampled_cancer_ligand_list_uc, sampled_normal_receptor_list_uc)[0], 4)
        corr_connected=round(corr_coef(cancer_ligand_list, normal_receptor_list)[0], 4)
        sig_corr_connected=round(corr_coef(cancer_ligand_list, normal_receptor_list)[1], 4)
        sig_corr_unconnected=round(corr_coef(sampled_cancer_ligand_list_uc, sampled_normal_receptor_list_uc)[1], 4)
        difference=np.round(corr_connected-corr_unconnected,5)
        #print('corr_unconnected')
        #print(corr_unconnected)
        #print('Assertion check that difference-corr_connected<0.01')
        #assert(difference-corr_connected<0.01)
        connected_df=pd.DataFrame()
        unconnected_df=pd.DataFrame()
        connected_df['normal_receptor_list']=normal_receptor_list
        connected_df['cancer_ligand_list']=cancer_ligand_list
        unconnected_df['normal_receptor_list_uc']=sampled_normal_receptor_list_uc
        unconnected_df['cancer_ligand_list_uc']=sampled_cancer_ligand_list_uc
        # Calculate the degrees of freedom
        n_con_cancer=len(cancer_ligand_list)
        n_uncon_cancer=len(sampled_cancer_ligand_list_uc)
        n_con_normal=len(normal_receptor_list)
        n_uncon_normal=len(sampled_normal_receptor_list_uc)
        ratio_cancer=np.mean(cancer_ligand_list)/np.mean(cancer_ligand_list_uc)
        ratio_normal=np.mean(normal_receptor_list)/np.mean(normal_receptor_list_uc)
        mean_cancer_ligand=np.mean(cancer_ligand_list)
        mean_normal_receptor=np.mean(normal_receptor_list)
        max_cancer_ligand=np.max(cancer_ligand_list)
        max_cancer_ligand_uc=np.max(cancer_ligand_list_uc)
        max_normal_receptor=np.max(normal_receptor_list)
        max_normal_receptor_uc=np.max(normal_receptor_list_uc)
        min_cancer_ligand=np.min(cancer_ligand_list)
        min_cancer_ligand_uc=np.min(cancer_ligand_list_uc)
        min_normal_receptor=np.min(normal_receptor_list)
        min_normal_receptor_uc=np.min(normal_receptor_list_uc)
        if ((corr_connected > 0.05) or (corr_connected < -0.05)):
            print('Plotting scatterplot for potential ligand-receptor pair')
            plt.close()
            #plt.title('C->'+Normal_cells[cell_subtype]+' '+lig_rec_pair[0]+' - '+lig_rec_pair[1])
            #plt.scatter(normal_receptor_list,cancer_ligand_list,color='blue',label='Connected pairs, correlation= '+str(corr_connected), alpha=0.5)
            plt.figure(figsize=(10, 5))
            plt.subplot(1, 2, 1)
            plt.violinplot(cancer_ligand_list,showmedians=True)
            #plt.text(0.5, -0.4, lig_rec_pair[0], horizontalalignment='center', verticalalignment='center', fontsize=12)
            plt.text(0.5, min(cancer_ligand_list) - 0.1, lig_rec_pair[0], horizontalalignment='center', verticalalignment='center', fontsize=12)
            plt.subplot(1, 2, 2)
            plt.violinplot(normal_receptor_list,showmedians=True)
            plt.text(0.5, min(normal_receptor_list) - 0.1, lig_rec_pair[1], horizontalalignment='center', verticalalignment='center', fontsize=12)
            #plt.text(0.5, -0.4, lig_rec_pair[1], horizontalalignment='center', verticalalignment='center', fontsize=12)
            #plt.xlim(0,7)
            #plt.ylim(0,7)
            #plt.legend()
            plt.savefig(path+Normal_cells[cell_subtype]+'_imgs/'+'Violin_plot_'+Normal_cells[cell_subtype]+'_'+ct_level+ '_'+ version+'_'+lig_rec_pair[0]+'-'+lig_rec_pair[1]+'.png')
            plt.close()
            plt.hist2d(normal_receptor_list, cancer_ligand_list, bins=30, cmap='viridis')
            plt.colorbar(label='Frequency')
            plt.xlabel(lig_rec_pair[1])
            plt.ylabel(lig_rec_pair[0])
            plt.title('2D Histogram of ligand-receptor expression')
            #plt.xlim(0,7)
            #plt.ylim(0,7)
            plt.savefig(path+Normal_cells[cell_subtype]+'_imgs/'+'2D_hist_'+Normal_cells[cell_subtype]+'_'+ct_level+ '_'+ version+'_'+lig_rec_pair[0]+'-'+lig_rec_pair[1]+'.png')
            #
            #plt.scatter(sampled_normal_receptor_list_uc,sampled_cancer_ligand_list_uc,color='orange',label='Unconnected pairs, correlation= '+str(corr_unconnected), alpha=0.5)
            try:
                sns.kdeplot(data=connected_df, x='normal_receptor_list',y='cancer_ligand_list', color="blue")
                #sns.kdeplot(data=unconnected_df, x='sampled_normal_receptor_list_uc',y='sampled_cancer_ligand_list_uc', color="orange")
            #    sns.kdeplot(normal_receptor_list,cancer_ligand_list,color="blue")
            #    sns.kdeplot(sampled_normal_receptor_list_uc,sampled_cancer_ligand_list_uc,color='orange')
            except ValueError as e:
                print("Error occurred while plotting KDE:", e)
                print("Gene pair responsible is:\n")
                print(lig_rec_pair[0]+' - '+lig_rec_pair[1])
            plt.close()
            plt.title('C->'+Normal_cells[cell_subtype]+' '+lig_rec_pair[0]+' - '+lig_rec_pair[1])
            plt.scatter(normal_receptor_list,cancer_ligand_list,color='blue',label='Connected pairs, correlation= '+str(corr_connected), alpha=0.5)
            plt.xlabel(Normal_cells[cell_subtype]+ ' mean = '+str(np.mean(normal_receptor_list)))
            plt.ylabel('Epithelial '+' mean = '+str(np.mean(cancer_ligand_list)))
            plt.legend()
            if (fake):
               if ((corr_connected>0.8) or (corr_connected < -0.8)):
                  plt.savefig(path+Normal_cells[cell_subtype]+'_imgs/'+'Connected_vs_unconnected_pair_'+Normal_cells[cell_subtype]+'_'+ct_level+ '_'+ version+'_'+lig_rec_pair[0]+'-'+lig_rec_pair[1]+'_cancer_ligand.png')
            else:
               plt.savefig(path+Normal_cells[cell_subtype]+'_imgs/'+'Connected_vs_unconnected_pair_'+Normal_cells[cell_subtype]+'_'+ct_level+ '_'+ version+'_'+lig_rec_pair[0]+'-'+lig_rec_pair[1]+'_cancer_ligand_DE.png')
    return difference, corr_connected, sig_corr_connected, corr_unconnected, sig_corr_unconnected, n_con_cancer, n_uncon_cancer, n_con_normal, n_uncon_normal, mean_cancer_ligand, mean_normal_receptor, ratio_cancer, ratio_normal, max_cancer_ligand, max_cancer_ligand_uc, max_normal_receptor, max_normal_receptor_uc, min_cancer_ligand, min_cancer_ligand_uc, min_normal_receptor, min_normal_receptor_uc

#file='anndata_GR_with_delaunay_40_microns.h5ad'

#print('Reading in anndata object')

h5ad_file='anndata_GR_with_delaunay_40_microns_'+ct_level+'_'+version+'.h5ad'

#print('Reading in '+path+file)

#adata=ad.read_h5ad(path+file)

adata=ad.read_h5ad(path+h5ad_file)

print('Read in the following anndata object '+ str(adata))

#Cancer_cell='Epithelial'

if (args.ct_level=='ct_level1') or (args.ct_level=='ct_level0'): #or ct_level=='annot_level0'):
    Cancer_cell='Epithelial'
    print('Doing Epithelial cell')
if args.ct_level=='annot_level0':
    Cancer_cell='Epithelial cells'
    print('Doing Epithelial cells')
if args.ct_level=='ct_clean':
    Cancer_cell=='Cancer_cell'
    print('Doing Cancer_cell, not Epithelial cell')

print('Cancer cell type is '+Cancer_cell)

#print('Now changing Cancer cell type to Epithelial')

#Cancer_cell='Epithelial'

#print('Cancer cell type is '+Cancer_cell)

Normal_cells = np.unique(adata.obs[ct_level])[[x!=Cancer_cell for x in np.unique(adata.obs[ct_level])]]

print('Normal cells:\n')
print(Normal_cells)

cell_idx=np.where(Normal_cells==cell_subtype)[0][0]

print('cell_idx')
print(cell_idx)

directory_name=path+Normal_cells[cell_idx]+'_imgs/'

print('directory_name')
print(directory_name)

# Check if the directory exists
if not os.path.exists(directory_name):
    # Create the directory if it doesn't exist
    os.makedirs(directory_name)

headers = ["Ligand", "Receptor", "difference", "corr_connected", "sig_corr_connected", "corr_unconnected","sig_corr_unconnected","n_con_cancer","n_uncon_cancer","n_con_normal","n_uncon_normal","mean_cancer_ligand","mean_normal_receptor", "ratio_cancer","ratio_normal","max_cancer_ligand","max_cancer_ligand_uc","max_normal_receptor","max_normal_receptor_uc","min_cancer_ligand","min_cancer_ligand_uc","min_normal_receptor","min_normal_receptor_uc"]

if (fake):
    print('Note: reading in all gene pairs')
    file_path=path+'all_gene_pairs_'+ct_level+'_'+version+'_'+sample+'.txt'
if (ligrec_known):
    print('Note:reading in known ligrec pairs')
    file_path=path+cell_subtype+'_gene_pairs_'+direction+'_'+ct_level+'_'+version+'_'+sample+'_known_ligrec.txt' 
else:
    print('Note: reading in true gene pairs')
    file_path=path+cell_subtype+'_gene_pairs_'+direction+'_'+ct_level+'_'+version+'_'+sample+'_cat.txt'

print('The gene combinations file is:\n')
 
print(file_path)
 
file_contents=[]
 
# Open the text file for reading
with open(file_path, 'r') as file:
    # Read each line in the file
    #lines = file.readlines()
    for line in file:
        genes = line.strip().split()
        if len(genes) >= 2:
            gene1, gene2 = genes[0],genes[1]
            file_contents.append((gene1, gene2))
        else:
            print(f"Issue in line: {line.strip()}. Not enough genes or invalid input format.")    
        # Split the line into gene1 and gene2 using a tab separator (assuming the file uses tabs)
        #gene1, gene2 = line.strip().split(' ')
        # Create a tuple with the gene pair and append it to the list
        #file_contents.append((gene1, gene2))
 
combination=list(file_contents)

#combination=combination[0:20]

print('Length of split combination:\n')
print(len(combination))

print('start_index')
print(start_index)
print('end_index')
print(end_index)

def process_split(combination,start,end):
    data=[]
    if fake:
        print('Looking at fake ligand-receptor pairs')
    else:
        print('Looking at ligand-receptor pairs from DE analysis')
    for j in range(start,end):
        print(f"doing {combination[j]}")
        difference, corr_connected, sig_corr_connected, corr_unconnected, sig_corr_unconnected, n_con_cancer, n_uncon_cancer, n_con_normal, n_uncon_normal, mean_cancer_ligand, mean_normal_receptor, ratio_cancer, ratio_normal, max_cancer_ligand, max_cancer_ligand_uc, max_normal_receptor, max_normal_receptor_uc, min_cancer_ligand, min_cancer_ligand_uc, min_normal_receptor, min_normal_receptor_uc=lig_rec_plot_normal_cancer2(adata,cell_subtype=cell_idx,impute=True,direction=direction,lig_rec_pair=combination[j],fake=fake,corr_coef=correlation_function)            #ttion_functionrue_difference_list.append(difference_true)
        data.append([combination[j][0], combination[j][1], difference, corr_connected, sig_corr_connected, corr_unconnected, sig_corr_unconnected, n_con_cancer, n_uncon_cancer, n_con_normal, n_uncon_normal,mean_cancer_ligand, mean_normal_receptor, ratio_cancer, ratio_normal, max_cancer_ligand, max_cancer_ligand_uc, max_normal_receptor, max_normal_receptor_uc, min_cancer_ligand, min_cancer_ligand_uc, min_normal_receptor, min_normal_receptor_uc])

    return data

# Process the split
data1=process_split(combination, start_index, end_index)

print(data1)

df = pd.DataFrame(data1, columns=["Ligand", "Receptor", "difference", "corr_connected", "sig_corr_connected", "corr_unconnected","sig_corr_unconnected","n_con_cancer","n_uncon_cancer","n_con_normal","n_uncon_normal","mean_cancer_ligand","mean_normal_receptor", "ratio_cancer","ratio_normal","max_cancer_ligand","max_cancer_ligand_uc","max_normal_receptor","max_normal_receptor_uc","min_cancer_ligand","min_cancer_ligand_uc","min_normal_receptor","min_normal_receptor_uc"])

print(df)
        
df.to_csv(path+Normal_cells[cell_idx]+'_imgs/'+Normal_cells[cell_idx]+'_'+ct_level+ '_'+ version+'_'+direction+'_main'+'_'+fake_flag+'_'+'part_'+str(split_part)+'_complete_cat.csv', index=False)

print('Finished run successfully')

# Record end time
end_time = time.time()

# Calculate elapsed time
elapsed_time = end_time - start_time
print(f"Elapsed time: {elapsed_time} seconds")

