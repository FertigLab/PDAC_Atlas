#!/bin/bash
#$ -N SG_Novogene_Alignment
#$ -pe local 8
#$ -l h_vmem=5G,mem_free=5G
#$ -m e
#$ -M jtandur1@jhu.edu

module load salmon/1.9.0

salmon quant -i SG_Novogene_12192022/Salmon_Output/salmon_partial_sa_index__default -l A -1 SG_Novogene_12192022/usftp21_novogene_com/RawData/JHH348CAF/JHH348CAF_1.fq.gz -2 SG_Novogene_12192022/usftp21_novogene_com/RawData/JHH348CAF/JHH348CAF_2.fq.gz -p 8 --validateMappings -o SG_Novogene_12192022/Salmon_Output/quants/JHH348CAF
salmon quant -i SG_Novogene_12192022/Salmon_Output/salmon_partial_sa_index__default -l A -1 SG_Novogene_12192022/usftp21_novogene_com/RawData/JHH348CoCAF/JHH348CoCAF_1.fq.gz -2 SG_Novogene_12192022/usftp21_novogene_com/RawData/JHH348CoCAF/JHH348CoCAF_2.fq.gz -p 8 --validateMappings -o SG_Novogene_12192022/Salmon_Output/quants/JHH348CoCAF
salmon quant -i SG_Novogene_12192022/Salmon_Output/salmon_partial_sa_index__default -l A -1 SG_Novogene_12192022/usftp21_novogene_com/RawData/JHH348CoOrg/JHH348CoOrg_1.fq.gz -2 SG_Novogene_12192022/usftp21_novogene_com/RawData/JHH348CoOrg/JHH348CoOrg_2.fq.gz -p 8 --validateMappings -o SG_Novogene_12192022/Salmon_Output/quants/JHH348CoOrg
salmon quant -i SG_Novogene_12192022/Salmon_Output/salmon_partial_sa_index__default -l A -1 SG_Novogene_12192022/usftp21_novogene_com/RawData/JHH348Org/JHH348Org_1.fq.gz -2 SG_Novogene_12192022/usftp21_novogene_com/RawData/JHH348Org/JHH348Org_2.fq.gz -p 8 --validateMappings -o SG_Novogene_12192022/Salmon_Output/quants/JHH348Org
salmon quant -i SG_Novogene_12192022/Salmon_Output/salmon_partial_sa_index__default -l A -1 SG_Novogene_12192022/usftp21_novogene_com/RawData/JHH362CAF/JHH362CAF_1.fq.gz -2 SG_Novogene_12192022/usftp21_novogene_com/RawData/JHH362CAF/JHH362CAF_2.fq.gz -p 8 --validateMappings -o SG_Novogene_12192022/Salmon_Output/quants/JHH362CAF
salmon quant -i SG_Novogene_12192022/Salmon_Output/salmon_partial_sa_index__default -l A -1 SG_Novogene_12192022/usftp21_novogene_com/RawData/JHH362CoCAF/JHH362CoCAF_1.fq.gz -2 SG_Novogene_12192022/usftp21_novogene_com/RawData/JHH362CoCAF/JHH362CoCAF_2.fq.gz -p 8 --validateMappings -o SG_Novogene_12192022/Salmon_Output/quants/JHH362CoCAF
salmon quant -i SG_Novogene_12192022/Salmon_Output/salmon_partial_sa_index__default -l A -1 SG_Novogene_12192022/usftp21_novogene_com/RawData/JHH362CoOrg/JHH362CoOrg_1.fq.gz -2 SG_Novogene_12192022/usftp21_novogene_com/RawData/JHH362CoOrg/JHH362CoOrg_2.fq.gz -p 8 --validateMappings -o SG_Novogene_12192022/Salmon_Output/quants/JHH362CoOrg
salmon quant -i SG_Novogene_12192022/Salmon_Output/salmon_partial_sa_index__default -l A -1 SG_Novogene_12192022/usftp21_novogene_com/RawData/JHH362Org/JHH362Org_1.fq.gz -2 SG_Novogene_12192022/usftp21_novogene_com/RawData/JHH362Org/JHH362Org_2.fq.gz -p 8 --validateMappings -o SG_Novogene_12192022/Salmon_Output/quants/JHH362Org
salmon quant -i SG_Novogene_12192022/Salmon_Output/salmon_partial_sa_index__default -l A -1 SG_Novogene_12192022/usftp21_novogene_com/RawData/JHH380CAF/JHH380CAF_1.fq.gz -2 SG_Novogene_12192022/usftp21_novogene_com/RawData/JHH380CAF/JHH380CAF_2.fq.gz -p 8 --validateMappings -o SG_Novogene_12192022/Salmon_Output/quants/JHH380CAF
salmon quant -i SG_Novogene_12192022/Salmon_Output/salmon_partial_sa_index__default -l A -1 SG_Novogene_12192022/usftp21_novogene_com/RawData/JHH380CoCAF/JHH380CoCAF_1.fq.gz -2 SG_Novogene_12192022/usftp21_novogene_com/RawData/JHH380CoCAF/JHH380CoCAF_2.fq.gz -p 8 --validateMappings -o SG_Novogene_12192022/Salmon_Output/quants/JHH380CoCAF
salmon quant -i SG_Novogene_12192022/Salmon_Output/salmon_partial_sa_index__default -l A -1 SG_Novogene_12192022/usftp21_novogene_com/RawData/JHH380Org/JHH380Org_1.fq.gz -2 SG_Novogene_12192022/usftp21_novogene_com/RawData/JHH380Org/JHH380Org_2.fq.gz -p 8 --validateMappings -o SG_Novogene_12192022/Salmon_Output/quants/JHH380Org
salmon quant -i SG_Novogene_12192022/Salmon_Output/salmon_partial_sa_index__default -l A -1 SG_Novogene_12192022/usftp21_novogene_com/RawData/mix_JHH380CoOrg/mix_JHH380CoOrg_1.fq.gz -2 SG_Novogene_12192022/usftp21_novogene_com/RawData/mix_JHH380CoOrg/mix_JHH380CoOrg_2.fq.gz -p 8 --validateMappings -o SG_Novogene_12192022/Salmon_Output/quants/mix_JHH380CoOrg
