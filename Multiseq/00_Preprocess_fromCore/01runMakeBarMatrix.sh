#!/bin/sh
#$ -N makeBarMatrix
#$ -S /bin/sh
#$ -cwd
#$ -q ngsc,ngsclargememory
#$ -l mem_free=50G,h_vmem=60G

. /etc/profile.d/modules.sh

module load sharedapps python

samplename=`echo $in | sed -e s/^..._//`

python /data/ngsc2data/NGSC_softwares/Multiseq10x-master/scripts/makeBarMatrix.py -C /data/ngsc4data/JZ02JHU501/JZ02JHU501_000_analysis/cellranger/count/${samplename}_GE/filtered_feature_bc_matrix/barcodes.txt -R readTable.csv -B ssSeq_Barcode.txt -O barMatrix.csv
