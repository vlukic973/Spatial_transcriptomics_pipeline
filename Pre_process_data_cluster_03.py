print('Running pre-process_data_cluster script')

print('Importing packages')

import scanpy as sc
#import squidpy as sq
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import anndata as ad
from pathlib import Path
import scipy
from scipy.sparse import csr_matrix, issparse
from matplotlib.pyplot import rc_context
import argparse
#import destvi_utils
import scvi
sc.logging.print_header()

# Create the argument parser
parser = argparse.ArgumentParser(description='Pre-process GR data')

# Add the --path argument
parser.add_argument('--path', type=str, help='Path to files')
parser.add_argument('--version', type=str, help='Run version')
parser.add_argument('--max_epochs', type=int, help='Max number of epochs')

# Parse the command-line arguments
args = parser.parse_args()

# Access the value of the --path argument
path = args.path
version = args.version
max_epochs = args.max_epochs

# Use the value of the --path argument in your script
# For example:
print(f'The specified path is: {path}')

file='adata_'+version+'.h5ad'

adata=ad.read_h5ad(path+file)
#sample='sample_2'

print('adata:\n')
print(adata)

adata.X=adata.layers['counts']
print('Using adata.X=adata.layers[counts]')

print('Checking if adata.X is sparse')
is_sparse=issparse(adata.X)
if (is_sparse):
   print('adata.X is sparse')
   print(np.sum(adata.X.toarray()==0)/(adata.shape[0]*adata.shape[1]))
else:
   print('adata.X is not sparse')
   print(print(np.sum(np.asarray(adata.X)==0)/(adata.shape[0]*adata.shape[1])*100))

scvi.model.SCVI.setup_anndata(adata)
vae = scvi.model.SCVI(adata)

print('Training vae using vae.train(train_size=0.9,max_epochs=400,early_stopping_patience=10,enable_checkpointing=True,early_stopping_monitor=\'elbo_validation\')\n')
vae.train(train_size=0.9,max_epochs=max_epochs,early_stopping_patience=10,enable_checkpointing=True,early_stopping_monitor='elbo_validation')

lower_case_genes=[x.lower() for x in adata.var_names]

training_history = vae.history_['elbo_train']

print('Plotting training loss')
epochs = range(1, len(training_history) + 1)
plt.plot(epochs, training_history, label='Training Loss')
plt.xlabel('Epochs')
plt.ylabel('Loss')
plt.title('Training Loss Over Epochs')
plt.legend()
plt.savefig(path+'Training_loss_SCVI_impute_'+version+'.png')

adata.obsm["X_normalized_scVI"] = vae.get_normalized_expression(library_size='latent')

print('Here min of adata.obsm["X_normalized_scVI"]')
print(np.min(adata.obsm["X_normalized_scVI"]))
print('Here max of adata.obsm["X_normalized_scVI"]')
print(np.max(adata.obsm["X_normalized_scVI"]))

print('Writing imputed anndata object')
print(path+'SCVI_imputed_anndata_GR_'+version+'.h5ad')
adata.write_h5ad(path+
    'SCVI_imputed_anndata_GR_'+version+'.h5ad'
)
