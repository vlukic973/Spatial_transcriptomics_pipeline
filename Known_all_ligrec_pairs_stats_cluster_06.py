import numpy as np
import pandas as pd
import scipy
from statsmodels.stats.multitest import multipletests
from scipy import stats
import argparse
import glob
import os

parser = argparse.ArgumentParser(description='Find genes expressed in cancer cells connected to normal cells, and vice versa, when splitting up a large list of gene combinations')

# Add the --path argument
parser.add_argument('--path', type=str, help='Path to files')
parser.add_argument('--cell_subtype', type=str, help='Normal cell subtype')
parser.add_argument('--direction', type=str, help='Cancer cell is ligand or receptor?')
parser.add_argument('--ct_level', type=str, required=True, help='Cell type annotation')
parser.add_argument('--version', type=str, required=True, help='Run version')
parser.add_argument('--sample',type=str,help="Which sample")
parser.add_argument('--sig_level',type=str,help="Significance level")

args = parser.parse_args()

# Access the value of the --path argument
path = args.path
cell_subtype = args.cell_subtype
direction = args.direction
ct_level=args.ct_level
version=args.version
sample=args.sample
sig_level=args.sig_level
save_version=version

# Find all matching files for CN and NC
print('Reading in files')
files = glob.glob(os.path.join(path, f'{cell_subtype}_imgs/{cell_subtype}_{ct_level}_{version}_{direction}__main_known_ligrec_genes_part_*_subset_full_complete_connected.csv'))

print('Here is the path:\n')
print(os.path.join(path, f'{cell_subtype}_imgs/'))

print('Here are the files:\n')
print(files)

# Function to read and concatenate CSV files
def read_and_concatenate(file_list):
    df_list = [pd.read_csv(file) for file in file_list]
    print('In read_and_concatenate function')
    print('Here is df_list')
    print(df_list)
    concatenated_df = pd.concat(df_list, ignore_index=True)
    return concatenated_df

df_ligrec = read_and_concatenate(files)

print('Here is df_ligrec:\n')
print(df_ligrec)

print('Here is df_ligrec columns:\n')
print(df_ligrec.columns)

df_ligrec.to_csv(path+cell_subtype+'_'+sample+'_concatenated_'+direction+'_'+ct_level+'_'+save_version+'.csv')

##### Filter by above average expression for both ligand and receptor

if (direction=='CN'):
   column_name_ligand='mean_cancer_ligand'
   column_name_receptor='mean_normal_receptor'
   #df_ligrec = df_ligrec.dropna(subset=['mean_cancer_ligand', 'mean_normal_receptor', 'mean_normal_ligand', 'mean_cancer_ligand'])
if (direction=='NC'):
   column_name_ligand='mean_normal_ligand'
   column_name_receptor='mean_cancer_receptor'

df_ligrec = df_ligrec.dropna(subset=[column_name_ligand, column_name_receptor])

print('Successfully dropped Na from df_ligrec')

# Create ligand and receptor dataframes with unique entries
ligand_df = df_ligrec[['Ligand', column_name_ligand]].drop_duplicates(subset=['Ligand'])
receptor_df = df_ligrec[['Receptor', column_name_receptor]].drop_duplicates(subset=['Receptor'])

print('Dropped duplicates')

# Ensure valid values are present before calculating means
if not ligand_df[column_name_ligand].empty:
    mean_ligand = np.mean(ligand_df[column_name_ligand])
else:
    mean_ligand = np.nan

if not receptor_df[column_name_receptor].empty:
    mean_receptor = np.mean(receptor_df[column_name_receptor])
else:
    mean_receptor = np.nan

# Apply the filter based on the mean values if they are valid
if not np.isnan(mean_ligand) and not np.isnan(mean_receptor):
    df_ligrec = df_ligrec[
        (df_ligrec[column_name_ligand] > mean_ligand) &
        (df_ligrec[column_name_receptor] > mean_receptor)
    ]

assert(np.unique(ligand_df[ligand_df['Ligand']=='COL4A1'][column_name_ligand])[0]==ligand_df[ligand_df['Ligand']=='COL4A1'][column_name_ligand].iloc[0])

assert(np.unique(receptor_df[receptor_df['Receptor']=='ITGA1'][column_name_receptor])[0]==receptor_df[receptor_df['Receptor']=='ITGA1'][column_name_receptor].iloc[0])

##### z-test

df_ligrec['fishers_z']=0.5*np.log((1+df_ligrec.corr_connected)/(1-df_ligrec.corr_connected))

df_ligrec['z']=df_ligrec['fishers_z']/(1/np.sqrt(np.unique(df_ligrec.n_con_cancer)[0]-3))

p_values = scipy.stats.norm.sf(abs(df_ligrec['z']))*2

df_ligrec['p_value']=p_values

df_ligrec=df_ligrec[df_ligrec['p_value']<float(sig_level)]

if not df_ligrec.empty:
	df_ligrec.to_csv(path+cell_subtype+'_'+sample+'_pval_'+direction+'_'+ct_level+'_'+str(sig_level)+'_'+save_version+'_known_ligrec.csv')
	is_rejected_fdr_bh, corrected_bh_p_values, _, _ = multipletests(df_ligrec['p_value'], alpha=float(sig_level), method='fdr_bh')
	is_rejected_fdr_by, corrected_by_p_values, _, _ = multipletests(df_ligrec['p_value'], alpha=float(sig_level), method='fdr_by')
	df_ligrec['corrected_bh_p_values']=corrected_bh_p_values
	df_ligrec['corrected_by_p_values']=corrected_by_p_values
	sig_rows_bh = df_ligrec[is_rejected_fdr_bh]
	sig_rows_by = df_ligrec[is_rejected_fdr_by]
	sig_rows_bh=sig_rows_bh.drop_duplicates()
	sig_rows_bh.to_csv(path+cell_subtype+'_'+sample+'_bh_correction_pval_'+direction+'_'+ct_level+'_'+str(sig_level)+'_'+save_version+'_known_ligrec.csv')
	sig_rows_bh_filt=sig_rows_bh[(sig_rows_bh[column_name_ligand]>mean_ligand) & (sig_rows_bh[column_name_receptor]>mean_receptor)]
	sig_rows_bh_filt.to_csv(path+cell_subtype+'_'+sample+'_bh_correction_pval_mean_filtered_'+direction+'_'+ct_level+'_'+str(sig_level)+'_'+save_version+'_known_ligrec.csv')
	print('sig_rows_bh_filt:\n')
	print(sig_rows_bh_filt)
	df_ligrec.to_csv(path+cell_subtype+'_'+sample+'_corrected_pval_'+direction+'_'+ct_level+'_'+str(sig_level)+'_'+save_version+'_known_ligrec.csv')
else:
	print('No significant results after p-value filtering.')

print('Script finished successfully')
