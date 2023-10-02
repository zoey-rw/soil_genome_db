# code for live-testing toolchest installation and data transfer
setwd("/projectnb2/talbot-lab-data/zrwerbin/toolchest/")
#devtools::install_github("trytoolchest/toolchest-client-r", force = T)

library(toolchest)
library(tidyverse)

toolchest::set_key("YzAxY2U.MDExZmY2MmQtYjNlYy00MTcxLWI2MGUtNjhjNTA5YzBmYjg2") # key from demo email

soil_sample_dir <- "/projectnb2/talbot-lab-data/jlopezna/second_metaGEM/qfiltered"
samp_names <- list.files(soil_sample_dir, recursive=F)
#file_names <- list.files(soil_sample_dir, pattern = "fastq", recursive=T)

out <- list()

samp_names <- c("HARV_001-M-25-38-20130709-gen","HARV_001-O-20160803-COMP-DNA1_mms","HARV_002-O-20180711-COMP-DNA1", "HARV_004-O-20131122-comp","HARV_004-O-31-8-20131122-gen","HARV_010-O-20160803-COMP-DNA1_mms")

samp_names <- c("HARV_004-O-20131122-comp","HARV_010-O-20160803-COMP-DNA1_mms")
samp = samp_names[[1]]

library(doParallel)
cl <- makeCluster(28, outfile="", type="PSOCK")
registerDoParallel(cl)

foreach(samp=samp_names, .errorhandling = 'pass') %dopar% {
	print(paste0("Starting Kraken2 loop for:", samp))

		library(pkgload)
		suppressMessages(library(toolchest))
		library(tidyverse)

		toolchest::set_key("YzAxY2U.MDExZmY2MmQtYjNlYy00MTcxLWI2MGUtNjhjNTA5YzBmYjg2") # key from demo email

	input_r1 <- file.path("/projectnb2/talbot-lab-data/zrwerbin/metagenomics/metaGEM/qfiltered", samp, paste0(samp, "_R1.fastq.gz"))
	input_r2 <- file.path("/projectnb2/talbot-lab-data/zrwerbin/metagenomics/metaGEM/qfiltered", samp, paste0(samp, "_R2.fastq.gz"))
	output_standard <- file.path("/projectnb2/talbot-lab-data/zrwerbin/toolchest/output/", paste0(samp, "_standard.tar.gz" ))
	output_pluspf <- file.path("/projectnb2/talbot-lab-data/zrwerbin/toolchest/output/", paste0(samp, "_pluspf.tar.gz" ))

	output_dir_standard <- file.path("/projectnb2/talbot-lab-data/zrwerbin/toolchest/output/kraken_standard/")
	output_dir_pluspf <- file.path("/projectnb2/talbot-lab-data/zrwerbin/toolchest/output/kraken_pluspf")

	kraken2(read_one = input_r1, read_two = input_r2,
					output_path = output_dir_standard, tool_args = c("report"),
					database_name = "standard")

	if(!file.exists(output_standard)){
		kraken2(read_one = input_r1, read_two = input_r2,
						output_path = output_standard, tool_args = c("report"),
						database_name = "standard")
	} else print(paste0("Already ran standard DB for:", samp))

	if(!file.exists(output_pluspf)){
	kraken2(read_one = input_r1, read_two = input_r2,
					output_path = output_pluspf, tool_args = c("report"),
					database_name = "PlusPF")
	} else print(paste0("Already ran plusPF DB for:", samp))

	return(paste0("Loop complete for:", samp))
}
stopCluster(cl)


# report <- untar(output_standard,files="k2_report.txt")
# report <- read.csv("/projectnb2/talbot-lab-data/zrwerbin/toolchest/output/k2_report.txt", sep="\t", header=F)
# untar(output_standard)
#
# #### read in outputs for percent classified ----
# out <- list()
# out_full <- list()
# for (samp in samp_names){
# 	output_standard <- file.path("/projectnb2/talbot-lab-data/zrwerbin/toolchest/output/", paste0(samp, "_standard.txt" ))
# 	output_pluspf <- file.path("/projectnb2/talbot-lab-data/zrwerbin/toolchest/output/", paste0(samp, "_pluspf.txt" ))
#
# 	standard <- data.table::fread(output_standard, header = F, nThread = 8)
# 	pluspf <- data.table::fread(output_pluspf, header = F, nThread = 8)
#
# 	standard_class <- table(standard$V1) %>% as.data.frame() %>%
# 		mutate(db = "Standard", freq = Freq/nrow(standard))
# 	pluspf_class <- table(pluspf$V1) %>% as.data.frame() %>%
# 		mutate(db = "PlusPF", freq = Freq/nrow(pluspf))
# 	#out[[samp]] <- do.call(rbind, c(standard_class, pluspf_class, pluspfp_class))
# 	out[[samp]] <- 	data.table::rbindlist(list(standard_class, pluspf_class))
# 	out[[samp]]$sampleID <- samp
# }
#
# df <- do.call(rbind, out)
# saveRDS(df, "/projectnb2/talbot-lab-data/zrwerbin/toolchest/kraken_pct_classified.rds")



