#!/bin/sh
#$ -N makeReadtable
#$ -S /bin/sh
#$ -cwd
#$ -q ngsc,ngsclargememory
#$ -l mem_free=50G,h_vmem=60G

. /etc/profile.d/modules.sh

module load sharedapps python

samplename=`echo $in | sed -e s/^..._//`

python /data/ngsc2data/NGSC_softwares/Multiseq10x-master/scripts/makeReadtable.py -C /data/ngsc4data/JZ02JHU501/JZ02JHU501_000_analysis/cellranger/count/${samplename}_GE/filtered_feature_bc_matrix/barcodes.txt -R1 /data/ngsc4data/JZ02JHU501/Multiseq_raw_data/${in}_barcodes/${samplename}_barcodes_R1.fastq.gz -R2 /data/ngsc4data/JZ02JHU501/Multiseq_raw_data/${in}_barcodes/${samplename}_barcodes_R2.fastq.gz -O readTable.csv
