library(tools)

# create shell script to run salmon on all samples from a folder
shellScriptName = 'salmonQuant.sh'
# path to file with previously created genome index
indexFile = 'SG_Novogene_12192022/Salmon_Output/salmon_partial_sa_index__default'
# path to folder with fastq files to run
dataFolder = 'SG_Novogene_12192022/usftp21_novogene_com/RawData'
# suffixes for the first end reads
read_suffix1 = c('_1')
# suffics for the second end reads
read_suffix2 = c('_2')
# file extansion (fastq or fastq.gz)
ext = '.fq.gz'
# output folder
out = 'SG_Novogene_12192022/Salmon_Output/quants'
# check if this folder exists, create if not
if(!file.exists(file.path(out))) dir.create(file.path(out))

##Extract samp names from foldernames
names <- list.dirs(dataFolder,full.names = F,recursive = F)

# open shell script
sink(file.path("SG_Novogene_12192022/Salmon_Output/",shellScriptName))
cat('#!/bin/bash\n')

cat('module load salmon/1.9.0')

# iterate sample names
for (i in names)
{
  # create paths to files for one end files
  read1_files = file.path(dataFolder,paste0(i,"/",i,read_suffix1,ext))
  # create paths to files for the second end files
  read2_files = file.path(dataFolder,paste0(i,"/",i,read_suffix2,ext))
  # create path to output files
  outFile = file.path(out,i)
  cmdln = paste('salmon quant -i', indexFile, # path to genome index file
                '-l A', # detect type of library automatically
                '-1', paste(read1_files, collapse = ' '), # the first read files
                '-2', paste(read2_files, collapse = ' '), # the second read files
                '-p 8 --validateMappings -o',outFile)  # -p is the number of threads and -o is output folder
  cat(cmdln,'\n')
}
# close shell script
sink()