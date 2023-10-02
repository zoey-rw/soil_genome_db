# Create list of high-abundance taxa, from mapping NEON metagenomes to SoilGenome DB



metadata_in <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/metagenomics/all_metagenome_metadata_2021.rds")


full_report_out <- list()

soil_sample_dir <- "/projectnb2/talbot-lab-data/zrwerbin/toolchest/output/reports/standard/"
samp_names <- list.files(soil_sample_dir, recursive=F)
samp_names <- gsub("_k2_report_standard.txt","",samp_names)
for (samp in samp_names){

	print(samp)
	standard_path <-	file.path("/projectnb2/talbot-lab-data/zrwerbin/toolchest/output/reports/standard", paste0(samp, "_k2_report_standard.txt"))

	if (file.exists(standard_path)){

		report_standard <- read_report2(standard_path) %>% mutate(sampleID = samp)


		summary_standard <- summarize_report_custom(report_standard) %>% mutate(database = "standard")
		#samples_summary <- data.table::rbindlist(list(summary_standard, summary_pluspf)) %>% mutate(name = samp)
		samples_summary <- summary_standard %>% mutate(sampleID = samp)
		#		samples_summary <- rbind(summary_standard, summary_pluspf) %>% mutate(name = samp)
		colnames(samples_summary) <- pavian:::beautify_string(colnames(samples_summary))

		raw_reads_column <- 1
		classified_reads_column <- 2
		microbial_reads_column <- 6
		samples_summary[, classified_reads_column:(ncol(samples_summary)-2)] <-
			signif(100 * sweep(samples_summary[, classified_reads_column:(ncol(samples_summary)-2)], 1, samples_summary[, raw_reads_column], `/`),3)
		out[[samp]] <- samples_summary
		full_report_out[[samp]] <- report_standard
	} else {
		print(paste0("skipping ", samp))
		next()
	}
}


all_samp_reports <- do.call(rbind, full_report_out)
saveRDS(all_samp_reports, "/projectnb2/talbot-lab-data/zrwerbin/toolchest/kraken_reports_merged.rds")




