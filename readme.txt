#### command line steps for creating, modifying, and using fungal genome database on BU's SCC ####

cd /projectnb/talbot-lab-data/zrwerbin/soil_genome_db/Struo2
OUTDIR=./fungal_genomes/

#### CREATE STRUO2 ENVIRONMENT ####
module load miniconda
conda activate mamba
mamba create -n struo2 -c bioconda
mamba env update -n struo2 --file /projectnb/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/conda_env.yaml
conda activate struo2

# Note: to speed things up I modified the kraken snakefile rules to use 28 threads, and a --fast-build parameter

#### DOWNLOAD LOTS OF GENOMES TO CREATE DATABASE ####
conda install -c bioconda ncbi-genome-download

#### RUN create_fungal_genome_list.R to create list of NCBI taxon IDs for download ####

# Checking with a single genome download
ncbi-genome-download -F fasta -o /projectnb/talbot-lab-data/zrwerbin/soil_genome_db/fungal_genomes -p 1 -r 3 -A GCA_001642055.1 -s genbank archaea,bacteria,fungi

# Downloading everything from corinne's list (~260 taxa)
ncbi-genome-download -F fasta -o /projectnb/talbot-lab-data/zrwerbin/soil_genome_db/fungal_genomes --taxids /projectnb/talbot-lab-data/zrwerbin/soil_genome_db/fungal_genomes/taxids.txt -s genbank archaea,bacteria,fungi

# Download everything (~2400 taxa)
ncbi-genome-download -F fasta -o /projectnb/talbot-lab-data/zrwerbin/soil_genome_db/fungal_genomes --taxids /projectnb/talbot-lab-data/zrwerbin/soil_genome_db/fungal_genomes/taxids_full.txt -s genbank archaea,bacteria,fungi --metadata-table /projectnb/talbot-lab-data/zrwerbin/soil_genome_db/fungal_genomes/ncbi_metadata_table.csv

wget --directory-prefix $OUTDIR https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz

#### RUN create_fungal_struo_input.R to reformat for Struo2 ####

# Use the taxonomizr R package to link NCBI accessions and taxonomy?
#install.packages("taxonomizr")
#library(taxonomizr)
#prepareDatabase('accessionTaxa.sql')

#### AFTER CREATING TSV AND DOWNLOADING FILES - RUN STRUO CREATE ####
# Takes a few hours - make sure you request enough computational resources

qrsh -l h_rt=18:00:00 -pe omp 28 -l mem_per_core=16G -n
cd /projectnb/talbot-lab-data/zrwerbin/soil_genome_db/Struo2
OUTDIR=./fungal_genomes/
module load miniconda
conda activate struo2

# Dry run, as a preview
snakemake --use-conda -j 1 --dry-run --configfile config.yaml --printshellcmds --rerun-incomplete
snakemake --use-conda -j 28 --configfile config.yaml --printshellcmds --rerun-incomplete --local-cores 28

# Nvm, instead run as batch job since it takes too long
qsub snakemake_sge.sh config.yaml 28


#### NOW USE DATABASE WITH METAGENOME FILES ####
qrsh -l h_rt=4:00:00 -pe omp 28 -l mem_per_core=4G -now n

module load kraken2
DBDIR=/projectnb/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/databases/fungal_db/test/kraken2
READ1=/projectnb2/talbot-lab-data/metabolic_models/scripts/metaGEM/qfiltered/HARV_020-O-20170731-COMP-DNA1/HARV_020-O-20170731-COMP-DNA1_R1.fastq.gz
READ2=/projectnb2/talbot-lab-data/metabolic_models/scripts/metaGEM/qfiltered/HARV_020-O-20170731-COMP-DNA1/HARV_020-O-20170731-COMP-DNA1_R2.fastq.gz
REPORT_PATH=/projectnb/talbot-lab-data/zrwerbin/soil_genome_db/kraken_output/test_HARV_020-O-20170731-COMP-DNA1.kreport
OUTPUT_PATH=/projectnb/talbot-lab-data/zrwerbin/soil_genome_db/kraken_output/test_HARV_020-O-20170731-COMP-DNA1.output

kraken2 --db $DBDIR --report $REPORT_PATH --output $OUTPUT_PATH --paired $READ1 $READ2 --threads $NTHREADS

#### RUN parse_kraken_output.r to combine multiple output files ####

#### RUN add_environmental_metadata.r to link with soil properties ####

