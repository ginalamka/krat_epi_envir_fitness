#!/bin/bash

##  Set username
USER=glamka

## Set dir name
DIR=THIS2_ms_kratepi_02_nfcore_methylseq  #_rerun

## Set project names
PROJ1=bismark  #bwameth

## Set error location
#ERR = work/


## --------------------------------
## Activate conda environment and load modules
## (For some reason, I get java-related errors if I rely on the config file for 
## these steps)

module load python/anaconda/3.10.9

source activate nextflow
module load java/15.0.1
module load singularity/3.8.4


## --------------------------------
## Download test data

### commented out by GFL 4/9/24 cuz already got the folders and such made
#mkdir /scratch/${USER}/${DIR}/
#chmod -R 777 /scratch/${USER}/${DIR}/
#mkdir /scratch/${USER}/${DIR}/${PROJ1}/
cd /scratch/${USER}/${DIR}/${PROJ1}/

#mkdir data
#cp -r ~/kratepi/ref_fasta .
#cp ~/kratepi/kratepi_sample_names_byhand.csv ./data/
#cp ~/kratepi/easleyconfig_gfl.conf .

#BELOW are old (new locations) and some seemingly not necessary as of 4/15/24
###cp -r /scratch/avrilh/krat_wgs_files/filtered_vcf .                                         #need to grab the vcf file
###cp -r /scratch/avrilh/krat_wgs_files/ref_fasta .                                            #need to grab the same ref Avril used
###cp /home/amh0254/krat_roh_analyses/sample_lists/krat_ms_samplesheet.csv ./data/              #permission denied -- grab_sample_names.sh grabs the names instead
###cp /home/amh0254/krat_roh_analyses/scripts/files/easleyconfig_amh.conf .                     #permission denied
###cp /home/gfl0003/krat_wgs_analyses/files/easleyconfig_amh.conf .                            #grab the file from elsewhere instead


## --------------------------------
## Run nf-score methylseq pipeline

#mkdir results

###removing 4/9/24 cuz doesnt seem to work and looks unnecessary
## aligning with bismark -- requires a manual fix after this step
nextflow run nf-core/methylseq \
--input data/kratepi_sample_names_byhand.csv \
--outdir results \
--fasta ref_fasta/dspec_genbank_assem.fa.gz \
--save_reference \
-resume \
-bg \
-c easleyconfig_gfl.conf \
-profile singularity 



## for some reason, ref .fa can't be used when it's just a symbolic link? idk. this fixes
## it tho, alongside specif.
#cd /scratch/${USER}/${DIR}/${PROJ1}/${ERR}
#rm dspec_genbank_assem.fa
#cp /scratch/avrilh/kratroh_01_assembindex/dspec_genbank_assem.fa

#cd /scratch/glamka/ms_kratepi_02_nfcore_methylseq/bismark/work/d5/082d525c5435b5851311157f54b344/BismarkIndex # /scratch/glamka/ms_kratepi_02_nfcore_methylseq_rerun/bismark/work/04/f271e46329a2d092e200bc2c09d84c/BismarkIndex
#rm dspec_genbank_assem.fa.gz
#cp /home/gfl0003/kratepi/ref_fasta/dspec_genbank_assem.fa /scratch/glamka/ms_kratepi_02_nfcore_methylseq/bismark/work/d5/082d525c5435b5851311157f54b344/BismarkIndex/dspec_genbank_assem.fa  # /scratch/glamka/ms_kratepi_02_nfcore_methylseq_rerun/bismark/work/04/f271e46329a2d092e200bc2c09d84c/BismarkIndex/dspec_genbank_assem.fa.gz 
###cd /scratch/avrilh/ms_kratroh_02_nfcore_methylseq/bismark/work/fc/e4a108a5ce95a58106631b57f1dd5b/BismarkIndex/
###rm dspec_genbank_assem.fa
###cp /scratch/glamka/kratroh_02_nfcore_methylseq/bismark_TEST/ref_fasta/dspec_genbank_assem.fa .             #/scratch/avrilh/kratroh_01_assembindex/dspec_genbank_assem.fa

# cd /scratch/${USER}/${DIR}/${PROJ1}/ 

#nextflow run nf-core/methylseq \
#--input data/kratepi_sample_names_byhand.csv \
#--outdir results \
#--fasta ref_fasta/dspec_genbank_assem.fa \
#--save_reference \
#-resume \
#-bg \
#-c easleyconfig_gfl.conf \
#-profile singularity \
#--aligner bwameth 

#consider if want to turn on the "comprehensive" switch - recommended for small genomes only


## below should be covered by config file
# conda deactivate
