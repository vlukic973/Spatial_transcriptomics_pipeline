#!/bin/bash
#SBATCH --job-name=ligrec_known
#SBATCH --output=%x.o%j
#SBATCH --error=%x.e%j   # Add this line to specify the error file
#SBATCH --ntasks=1
#SBATCH --time=08:00:00
#SBATCH --mem=120GB
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu_long

# To clean and load modules defined at the compile and link phases
module purge
module load anaconda3/2020.02/gcc-9.2.0

# Activate the conda environment
source activate scvi-env
#source activate scvi-env-cluster
#pip install scanpy

# Define the total number of tasks
total_tasks="$1"

# Define the number of gene pairs per task (10 gene pairs per task)
gene_pairs_per_task="$2"

path="$3"

cell_subtype="$4"

ct_level="$5"

version="$6"

corr_coef="$7"

sample="$8"

hfad_file="anndata_GR_with_delaunay_40_microns_${ct_level}_${version}.h5ad"

# Calculate the range for the task array
start_task=$SLURM_ARRAY_TASK_ID
end_task=$start_task
step=1

# Calculate the start and end gene pair indices for this task
start_index=$((($start_task - 1) * $gene_pairs_per_task))
end_index=$((($end_task * $gene_pairs_per_task)))

echo "task_id: $SLURM_ARRAY_TASK_ID"
echo "start_index: $start_index"
echo "end_index: $end_index"
echo "path: $path"
echo "cell_subtype: $cell_subtype"
echo "corr_coef: $corr_coef"
echo "sample: $sample"
echo "h5ad_file: $h5ad_file"

# Execute the Python script with task-specific arguments
python Ligrec_split_run_gene_pair_stats_ctlevel_connected_unconnected_all_dirs_expression_public_v7.py --path "$path" --cell_subtype "$cell_subtype" --ct_level "$ct_level" --version "$version" --corr_coef "$corr_coef" --sample "$sample" --ligrec_known --split_part "$SLURM_ARRAY_TASK_ID" --split_size "$gene_pairs_per_task" --start_index "$start_index" --end_index "$end_index" --hfad_file "$hfad_file"


