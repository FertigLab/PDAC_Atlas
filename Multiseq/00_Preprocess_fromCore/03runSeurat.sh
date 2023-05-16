#!/bin/sh
#$ -N seurat_umap
#$ -S /bin/sh
#$ -cwd
#$ -q ngsc
#$ -l mem_free=25G,h_vmem=30G

. /etc/profile.d/modules.sh

module load sharedapps conda 

export PATH="/home/agupta/miniconda3/bin:$PATH"
source activate /home/agupta/miniconda3/envs/R_v4.0.3

samplename=`echo $in | sed -e s/^..._//`

Rscript seurat.R -i /data/ngsc4data/JZ02JHU501/JZ02JHU501_000_analysis/cellranger/count/${samplename}_GE/filtered_feature_bc_matrix
