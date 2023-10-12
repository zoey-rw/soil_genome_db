#!/bin/bash -l


# Set SCC project
#$ -P dietzelab
#$ -l h_rt=4:00:00  # Specify hard time limit for the job.
#$ -pe omp 28 -l mem_per_core=4G
#$ -j y # Merge stderr into the stdout file, to reduce clutter
#$ -o run_kraken_test.log


cd /projectnb/talbot-lab-data/zrwerbin/soil_genome_db/Struo2
module load miniconda
module load kraken2
conda activate struo2

NTHREADS=28
DBDIR=/projectnb/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/databases/fungal_db/test/kraken2


READ1=/projectnb2/talbot-lab-data/metabolic_models/scripts/metaGEM/qfiltered/HARV_004-O-20170731-COMP-DNA1/HARV_004-O-20170731-COMP-DNA1_R1.fastq.gz
READ2=/projectnb2/talbot-lab-data/metabolic_models/scripts/metaGEM/qfiltered/HARV_004-O-20170731-COMP-DNA1/HARV_004-O-20170731-COMP-DNA1_R2.fastq.gz
REPORT_PATH=/projectnb/talbot-lab-data/zrwerbin/soil_genome_db/kraken_output/test_HARV_004-O-20170731-COMP-DNA1.kreport
OUTPUT_PATH=/projectnb/talbot-lab-data/zrwerbin/soil_genome_db/kraken_output/test_HARV_004-O-20170731-COMP-DNA1.output


kraken2 --db $DBDIR --report $REPORT_PATH --output $OUTPUT_PATH --paired $READ1 $READ2 --threads $NTHREADS


READ1=/projectnb2/talbot-lab-data/metabolic_models/scripts/metaGEM/qfiltered/HARV_001-O-20160803-COMP-DNA1/HARV_001-O-20160803-COMP-DNA1_R1.fastq.gz
READ2=/projectnb2/talbot-lab-data/metabolic_models/scripts/metaGEM/qfiltered/HARV_001-O-20160803-COMP-DNA1/HARV_001-O-20160803-COMP-DNA1_R2.fastq.gz
REPORT_PATH=/projectnb/talbot-lab-data/zrwerbin/soil_genome_db/kraken_output/test_HARV_001-O-20160803-COMP-DNA1.kreport
OUTPUT_PATH=/projectnb/talbot-lab-data/zrwerbin/soil_genome_db/kraken_output/test_HARV_001-O-20160803-COMP-DNA1.output

kraken2 --db $DBDIR --report $REPORT_PATH --output $OUTPUT_PATH --paired $READ1 $READ2 --threads $NTHREADS


READ1=/projectnb2/talbot-lab-data/metabolic_models/scripts/metaGEM/qfiltered/HARV_010-O-20160803-COMP-DNA1/HARV_010-O-20160803-COMP-DNA1_R1.fastq.gz
READ2=/projectnb2/talbot-lab-data/metabolic_models/scripts/metaGEM/qfiltered/HARV_010-O-20160803-COMP-DNA1/HARV_010-O-20160803-COMP-DNA1_R2.fastq.gz
REPORT_PATH=/projectnb/talbot-lab-data/zrwerbin/soil_genome_db/kraken_output/test_HARV_010-O-20160803-COMP-DNA1.kreport
OUTPUT_PATH=/projectnb/talbot-lab-data/zrwerbin/soil_genome_db/kraken_output/test_HARV_010-O-20160803-COMP-DNA1.output

kraken2 --db $DBDIR --report $REPORT_PATH --output $OUTPUT_PATH --paired $READ1 $READ2 --threads $NTHREADS

