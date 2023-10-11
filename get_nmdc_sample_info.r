# Testing how to download data from NMDC API without losing sample names
library(curl)
library(tidyverse)
library(httr)

# All NEON samples (from study ID)
#https://api.microbiomedata.org/docs#/find/find_data_objects_for_study_data_objects_study__study_id__get
test <- curl_download("https://api.microbiomedata.org/data_objects/study/nmdc%3Asty-11-34xj1150", destfile = "nmdc_samples.json")

read_samples_in <- jsonify::from_json("nmdc_samples.json") %>%  flatten
sample_names <- names(read_samples_in)
headers = c(`accept` = "application/json")

# Takes a couple minutes to loop through all 4443
sample_info_list <- list()
for (i in 1:length(sample_names)){
	print(i)
	sample_name_simple = gsub("nmdc:bsm-","",sample_names[[i]])
	url_to_call <- paste0("https://api.microbiomedata.org/biosamples/nmdc%3Absm-",sample_name_simple)
  sample_json <- httr::GET(url = url_to_call, httr::add_headers(.headers=headers))
  sample_info <- content(sample_json, encoding = "UTF-8") %>% unlist %>% cbind %>% t %>% as.data.frame()
  sample_info_list[[i]] <- sample_info
}
sample_info_df = plyr::rbind.fill(sample_info_list)

# Filter by processing institution - only want <100 samples from JGI (deeply-sequenced)
sample_info_jgi <- sample_info_df %>% filter(!is.na(gold_biosample_identifiers))

# Save for later
write_csv(sample_info_jgi, "nmdc_sample_info.csv")

# Now use API to pull downstream analysis products, a la https://microbiomedata.github.io/nmdc-runtime/nb/get_data/
# Not working with original python code either
id_biosample = "dobj-11-g33nhf73"
id_biosample = "bsm-11-9dcapv51"
rs_ompro = content(httr::GET(paste0("https://api.microbiomedata.org/activities?filter=type:nmdc:ReadBasedTaxonomyAnalysisActivity,has_input:",id_biosample)))





sample_info_nmdc <- read.csv("nmdc_sample_info.csv")

nmdc_kraken_reports <- list.files("/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/kraken_output/NMDC_pipeline_output")

nmdc_biosample <- gsub(".1_kraken2_report.tsv","",nmdc_kraken_reports)

recode_list <- sample_info_nmdc$id
names(recode_list) <- sample_info_nmdc$name

recode_list2 <- sample_info_nmdc$name
names(recode_list2) <- sample_info_nmdc$id
sampleID <- recode(nmdc_biosample, !!!recode_list2)


