library(tidyverse)

to_download_ids_in <- readLines("/projectnb/talbot-lab-data/zrwerbin/soil_genome_db/fungal_genomes/taxids.txt")

genome_info_in <- read_tsv("/projectnb/talbot-lab-data/zrwerbin/soil_genome_db/fungal_genomes/acc_table.tsv")

genome_info_in$file_path = file.path("/projectnb/talbot-lab-data/zrwerbin/soil_genome_db/fungal_genomes/genbank/fungi", genome_info_in$ncbi_genbank_assembly_accession)


already_downloaded <- list.files("/projectnb/talbot-lab-data/zrwerbin/soil_genome_db/fungal_genomes/genbank/fungi",full.names = T, #include.dirs = T,
																 recursive = T, pattern = ".fna.gz")

accession_no = basename(dirname(already_downloaded))

in_csv = accession_no[accession_no %in% genome_info_in$ncbi_genbank_assembly_accession]
in_csv_filepath = already_downloaded[accession_no %in% genome_info_in$ncbi_genbank_assembly_accession]

downloaded_file_info = genome_info_in[match(in_csv, genome_info_in$ncbi_genbank_assembly_accession),]
identical(downloaded_file_info$ncbi_genbank_assembly_accession, in_csv)
downloaded_file_info$fasta_file_path=in_csv_filepath

not_in_csv = accession_no[!accession_no %in% genome_info_in$ncbi_genbank_assembly_accession]

to_write = downloaded_file_info %>% select(ncbi_species_taxid = `ORGANISM NCBI TAX ID`,
																					 ncbi_organism_name = `ORGANISM NAME`,
																					 species = `ORGANISM NCBI SPECIES`,
																					 accession = ncbi_genbank_assembly_accession,
																					 fasta_file_path,ncbi_taxonomy)

# Remove any taxon ID with multiple genomes - THIS IS NOT A GOOD METHOD! INSTEAD USE THE NEWEST GENOME, OR BEST GENOME
to_write <- to_write %>% distinct(ncbi_species_taxid, .keep_all = TRUE)
#ncbi_taxonomy

dim(to_write)
write_tsv(to_write,  "/projectnb/talbot-lab-data/zrwerbin/soil_genome_db/fungal_genomes/fungal_struo.tsv")


library(taxonomizr)

taxonomizr::accessionToTaxa(to_write$accession, "accessionTaxa.sql")
taxonomizr::accessionToTaxa("GCA_000002515.1", "accessionTaxa.sql")
read.accession2taxid(list.files('.','accession2taxid.gz$'),'accessionTaxa.sql',indexTaxa=TRUE,overwrite=TRUE)
getAccessions(3702,'accessionTaxa.sql',limit=10)

getTaxonomy(to_write$ncbi_species_taxid,'accessionTaxa.sql')


already_downloaded <- list.files("/projectnb/talbot-lab-data/zrwerbin/soil_genome_db/fungal_genomes/",full.names = T, #include.dirs = T,
																 recursive = T, pattern = ".csv")
