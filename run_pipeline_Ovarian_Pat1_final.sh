#!/bin/bash
#SBATCH --job-name=OC_pipeline_1
#SBATCH --output=%x.o%j
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --mem=120GB
#SBATCH --cpus-per-task=4
#SBATCH --partition=cpu_long

# Config
paths=('/path/to/file/HumanOvarianCancerPatient1/')
output_log="pipeline_HumanOvarianCancerPatient1_v_f1.txt"
ct_level="ct_level0"
version="v1"
prefix="OvarianCancerPatient1_"
corr_coef="spearman"
sample="ovarian"
max_epochs=100
TOTAL_CELLS=358485
SPLIT_SIZE=23899

START_CELL=$((SLURM_ARRAY_TASK_ID * SPLIT_SIZE))
END_CELL=$((START_CELL + SPLIT_SIZE))

module purge
module load anaconda3/2020.02/gcc-9.2.0

for data_path in "${paths[@]}"; do
    
    source activate tangram_v3
    echo "Running Get_celltype_annotation_tangram_01.py --path $data_path ..."
    python Get_celltype_annotation_tangram_01.py \
        --path "$data_path" \
        --prefix "$prefix" \
        --file_sc "OCA_UNB_10X_E-MTAB-8107.h5ad" \
        --start_cell "$START_CELL" \
        --end_cell "$END_CELL" \
        >> "$output_log" 2>&1

    source activate scvi-env-cluster
    echo "==== Aggregate data and remap cells ====" | tee -a "$output_log"
    python Aggregate_adata_remap_cells_02.py --path "$data_path" --prefix "$prefix" --version "$version" >> "$output_log" 2>&1

    echo "==== Pre-process ====" | tee -a "$output_log"
    python Pre_process_data_cluster_03.py --path "$data_path" --version "$version" --max_epochs "$max_epochs" >> "$output_log" 2>&1

    echo "==== Build Delaunay Graphs ====" | tee -a "$output_log"
    python Generate_delaunay_graph_value_counts_04.py --path "$data_path" --ct_level "$ct_level" --version "$version" --prefix "$prefix" >> "$output_log" 2>&1
    
    echo "==== Run ligand-receptor pair detection across all cell types ====" | tee -a "$output_log"

    echo "==== Doing Endothelial cell type ====" | tee -a "$output_log"
    python Submit_ligand_receptor_celltype_scripts_05.py \
        --path "$data_path" --cell_subtype "Endothelial" --ct_level "$ct_level" \
        --version "$version" --corr_coef "$corr_coef" --sample "$sample" --ligrec_known >> "$output_log" 2>&1

    echo "==== Doing Fibroblast cell type ====" | tee -a "$output_log"
    python Submit_ligand_receptor_celltype_scripts_05.py \
        --path "$data_path" --cell_subtype "Fibroblast" --ct_level "$ct_level" \
        --version "$version" --corr_coef "$corr_coef" --sample "$sample" --ligrec_known >> "$output_log" 2>&1

    echo "==== Doing Myeloid cell type ====" | tee -a "$output_log"
    python Submit_ligand_receptor_celltype_scripts_05.py \
        --path "$data_path" --cell_subtype "Myeloid" --ct_level "$ct_level" \
        --version "$version" --corr_coef "$corr_coef" --sample "$sample" --ligrec_known >> "$output_log" 2>&1
        
    echo "==== Doing Lymphoctye cell type ====" | tee -a "$output_log"
    python Submit_ligand_receptor_celltype_scripts_05.py \
        --path "$data_path" --cell_subtype "Lymphocyte" --ct_level "$ct_level" \
        --version "$version" --corr_coef "$corr_coef" --sample "$sample" --ligrec_known >> "$output_log" 2>&1

    echo "==== Aggregate known ligrec pairs ====" | tee -a "$output_log"
#    python Known_all_ligrec_pairs_stats_cluster_v3.py --path "$data_path" --ct_level "$ct_level" --version "$version" >> "$output_log" 2>&1

    python Known_all_ligrec_pairs_stats_cluster_06.py --path "$data_path" --cell_subtype "Lymphocyte" --ct_level "$ct_level" --version "$version" --sample "$sample" --sig_level "0.05" --direction "CN" >> "$output_log" 2>&1

    python Known_all_ligrec_pairs_stats_cluster_06.py --path "$data_path" --cell_subtype "Lymphocyte" --ct_level "$ct_level" --version "$version" --sample "$sample" --sig_level "0.05" --direction "NC" >> "$output_log" 2>&1

    python Known_all_ligrec_pairs_stats_cluster_06.py --path "$data_path" --cell_subtype "Fibroblast" --ct_level "$ct_level" --version "$version" --sample "$sample" --sig_level "0.05" --direction "CN" >> "$output_log" 2>&1

    python Known_all_ligrec_pairs_stats_cluster_06.py --path "$data_path" --cell_subtype "Fibroblast" --ct_level "$ct_level" --version "$version" --sample "$sample" --sig_level "0.05" --direction "NC" >> "$output_log" 2>&1

    python Known_all_ligrec_pairs_stats_cluster_06.py --path "$data_path" --cell_subtype "Myeloid" --ct_level "$ct_level" --version "$version" --sample "$sample" --sig_level "0.05" --direction "CN" >> "$output_log" 2>&1

    python Known_all_ligrec_pairs_stats_cluster_06.py --path "$data_path" --cell_subtype "Myeloid" --ct_level "$ct_level" --version "$version" --sample "$sample" --sig_level "0.05" --direction "NC" >> "$output_log" 2>&1

    python Known_all_ligrec_pairs_stats_cluster_06.py --path "$data_path" --cell_subtype "Endothelial" --ct_level "$ct_level" --version "$version" --sample "$sample" --sig_level "0.05" --direction "CN" >> "$output_log" 2>&1

    python Known_all_ligrec_pairs_stats_cluster_06.py --path "$data_path" --cell_subtype "Endothelial" --ct_level "$ct_level" --version "$version" --sample "$sample" --sig_level "0.05" --direction "NC" >> "$output_log" 2>&1

done
