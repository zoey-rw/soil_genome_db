# For GOLD data from 20 Jun, 2022
# Downloaded from https://gold.jgi.doe.gov/downloads

library(readxl)
library(tidyverse)

gold_data_readme <- readxl::read_xlsx("/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/goldData.xlsx", sheet = 1)
gold_data_study <- readxl::read_xlsx("/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/goldData.xlsx", sheet = 2)
gold_data_biosample <- readxl::read_xlsx("/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/goldData.xlsx",  sheet = 3)
gold_data_organism <- readxl::read_xlsx("/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/goldData.xlsx",  sheet = 4)
gold_data_sequencing <- readxl::read_xlsx("/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/goldData.xlsx",  sheet = 5)
gold_data_analysis <- readxl::read_xlsx("/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/goldData.xlsx",  sheet = 6)

ecosystem_paths <- readxl::read_xlsx("/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/GOLDs5levelEcosystemClassificationPaths.xlsx")



# Subset to fungi
only_fungi_organism <- gold_data_organism %>% filter(`ORGANISM NCBI KINGDOM` == "Fungi")
only_fungi_organism_ids <- only_fungi_organism %>% select(`ORGANISM GOLD ID`) %>% unique() %>% unlist()

# Filter to genomes (analysis projects) from fungi
only_fungi_sequencing <- gold_data_sequencing %>% filter(`ORGANISM GOLD ID` %in% only_fungi_organism_ids)
only_fungi_sequencing_ids <- only_fungi_sequencing %>% select(`PROJECT GOLD ID`) %>% unlist()

# Filter to genomes (sequencing metadata) from soil organisms
# This file does not have all the completeness columns we want :(
only_fungi_ap <- gold_data_analysis %>% filter(`AP PROJECT GOLD IDS` %in% only_fungi_sequencing_ids)

genome_info <- merge(only_fungi_ap, only_fungi_sequencing, by.x = "AP PROJECT GOLD IDS", by.y = "PROJECT GOLD ID")
# Most of these columns can be removed, either before or after merging

genome_info <- merge(genome_info, only_fungi_organism, by = "ORGANISM GOLD ID")


# alternaria = genome_info_taxonomy[which(genome_info_taxonomy$`ORGANISM NCBI SPECIES`=="Alternaria alternata"),]
# genome_info <- alternaria

# Next: add accession info as a column
#accession_info <- genome_info_taxonomy$`AP GENBANK`
accession_info_test <- "[{\"genbankId\":\"ARFV00000000\",\"assemblyAccession\":\"GCA_000374445.1\"}]"
get_accession <- function(accession_info_test) {
	accession_info <- jsonlite::parse_json(accession_info_test)
	if (length(accession_info)==0) {
		return(NA)
	}
	out <- accession_info[[1]]
	print(out$assemblyAccession)
	to_return <- ifelse(is.null(out$assemblyAccession), NA, out$assemblyAccession)
	return(to_return)
}
get_accession(accession_info_test)

# Approach 1
genome_info$NCBI_accession1 <- lapply(genome_info$`AP GENBANK`, get_accession)

# Approach 2
genome_info$NCBI_accession <- NA
for (i in 1:nrow(genome_info)){
	print(i)
	genome_info$NCBI_accession[i] <- get_accession(genome_info$`AP GENBANK`[i])
}

# Check why ~1000 are NA
accession_NAs <- genome_info[is.na(genome_info$NCBI_accession),]

genome_info$ncbi_genbank_assembly_accession = genome_info$NCBI_accession

genome_info$ncbi_taxonomy = paste0("k__", genome_info$`ORGANISM NCBI KINGDOM`,
																	 ";p__",genome_info$`ORGANISM NCBI PHYLUM`,
																	 ";c__",genome_info$`ORGANISM NCBI CLASS`,
																	 ";o__",genome_info$`ORGANISM NCBI ORDER`,
																	 ";f__",genome_info$`ORGANISM NCBI FAMILY`,
																	 ";g__",genome_info$`ORGANISM NCBI GENUS`,
																	 ";s__",genome_info$`ORGANISM NCBI SPECIES`)
genome_info <- genome_info %>% select("ORGANISM GOLD ID",
																			"ORGANISM NAME",
																			"ORGANISM NCBI TAX ID" ,
																			"ORGANISM NCBI SPECIES",
																			"UTIL_STATUS" = "AP JGI DATA UTILIZATION STATUS\n (For projects sequenced at JGI only)",
																			"NCBI BIOSAMPLE ACCESSION",
																			"NCBI BIOPROJECT ACCESSION",
																			"ORGANISM NCBI SPECIES",
																			"STUDY GOLD ID",
																			"ncbi_genbank_assembly_accession","ncbi_taxonomy")

write_tsv(genome_info, "/projectnb/talbot-lab-data/zrwerbin/soil_genome_db/fungal_genomes/acc_table.tsv")

unique(df$ORGANISM.NCBI.TAX.ID)
TAXIDS

ncbi-genome-download -F fasta -o /projectnb/talbot-lab-data/zrwerbin/soil_genome_db/fungal_genomes -p 1 -r 3 -t 5599 -s genbank archaea,bacteria,fungi -d

corinne_list <- read_csv("/projectnb/talbot-lab-data/zrwerbin/soil_genome_db/corinne_files/23.7.26_WP21_genome_download_links__noSppDups_UPDATED.csv")

corinne_list_manual <- corinne_list %>% filter(Go.Download.Link=="manual")

#to_download <- genome_info_taxonomy %>% filter(`ORGANISM NCBI SPECIES` %in% corinne_list$Taxa)
to_download <- genome_info %>% filter(`ORGANISM NCBI SPECIES` %in% corinne_list_manual$Name)
to_download_ids <- unique(to_download$`ORGANISM NCBI TAX ID`) %>% unique() %>%  paste(sep="", collapse="\n")#
	#paste(sep="", collapse=",")

writeLines(to_download_ids,"/projectnb/talbot-lab-data/zrwerbin/soil_genome_db/fungal_genomes/taxids.txt")

#2385 IDs
to_download_ids_full <- unique(genome_info$`ORGANISM NCBI TAX ID`) %>% unique() %>%  paste(sep="", collapse="\n")#
writeLines(to_download_ids_full,"/projectnb/talbot-lab-data/zrwerbin/soil_genome_db/fungal_genomes/taxids_full.txt")


cmd <- paste0("ncbi-genome-download -F fasta -o /projectnb/talbot-lab-data/zrwerbin/soil_genome_db/fungal_genomes -p 1 -r 3 -t ", to_download_ids, " -s genbank fungi")

taxaId<-accessionToTaxa(c("LN847353.1","AL079352.3"),"accessionTaxa.sql")

getTaxonomy(taxaId,'accessionTaxa.sql')
