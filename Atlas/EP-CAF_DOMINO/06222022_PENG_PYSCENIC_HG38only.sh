#!/bin/bash
#SBATCH --job-name=PENG_PYSCENIC_HG38only
#SBATCH --time=72:00:00
#SBATCH --partition=defq
#SBATCH --nodes=1
# number of tasks (processes) per node
#SBATCH --ntasks-per-node=10
#SBATCH --mem=160G
#SBATCH --mail-type=end
#SBATCH --mail-user=jtandur1@jhu.edu

# directory for Inputs and results
DOMINO_DIR="Domino"
REFERENCE_DIR="Domino/Feather_Files/HG38"
RESULT_DIR="PDAC_EpCAF/PENG_HG38only"
COUNTS_LOOM="PDAC_EpCAF/cds_PENG_EpCAF.loom"

# run step1
singularity exec "${DOMINO_DIR}/aertslab-pyscenic-0.11.0.sif" pyscenic grn \
	$COUNTS_LOOM \
	"${DOMINO_DIR}/allTFs_hg38.txt" \
	-o "${RESULT_DIR}/adj_PDAC_PENG_EpCAF_hg38.tsv" \
	--num_workers 10

echo "Step one complete"

# run step2
singularity exec "${DOMINO_DIR}/aertslab-pyscenic-0.11.0.sif" pyscenic ctx \
	"${RESULT_DIR}/adj_PDAC_PENG_EpCAF_hg38.tsv" \
	"${REFERENCE_DIR}/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather" \
	"${REFERENCE_DIR}/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather" \
	--annotations_fname "${REFERENCE_DIR}/motifs-v9-nr.hgnc-m0.001-o0.0.tbl" \
	--expression_mtx_fname "${COUNTS_LOOM}" \
	--mode "dask_multiprocessing" \
	--output "${RESULT_DIR}/regulons_PDAC_PENG_EpCAF_hg38.csv" \
	--num_workers 10

echo "Step two complete"

# run step3
singularity exec "${DOMINO_DIR}/aertslab-pyscenic-0.11.0.sif" pyscenic aucell \
	"${COUNTS_LOOM}" \
	"${RESULT_DIR}/regulons_PDAC_PENG_EpCAF_hg38.csv" \
	-o "${RESULT_DIR}/auc_PDAC_PENG_EpCAF_hg38.csv"

echo "Step three complete"

echo "Finished with job $SLURM_JOBID"
