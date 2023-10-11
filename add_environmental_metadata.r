# Read in Kraken assignment data and combine with soil metadata
source("/projectnb/dietzelab/zrwerbin/microbialForecasts/microbialForecast/R/helperFunctions.r")
library(tidyverse)
library(broom)

#file_names <- list.files(soil_sample_dir, pattern = "fastq", recursive=T)
all_samp_reports <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/toolchest/kraken_reports_merged.rds")

samp_info <- parseNEONsampleIDs(as.character(all_samp_reports$sampleID))
all_samp_reports$plot_date = samp_info$plot_date

NEON_soil_phys <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/random_data/time_series_data/soil_phys_allsites.rds") %>%
	select("sampleID", "siteID", "plotID", "soilInCaClpH", "litterDepth", "soilTemp", "soilMoisture")
phys_info <- parseNEONsampleIDs(as.character(NEON_soil_phys$sampleID))
NEON_soil_phys$plot_date = phys_info$plot_date

NEON_soil_chem <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/random_data/time_series_data/soil_chem_allsites.rds") %>%
	select("siteID", "plotID", "organicCPercent", "nitrogenPercent", "CNratio", "sampleID")
chem_info <- parseNEONsampleIDs(as.character(NEON_soil_chem$sampleID))
NEON_soil_chem$plot_date = chem_info$plot_date

relevant_soil_data = merge(NEON_soil_phys, NEON_soil_chem) %>%
	filter(plot_date %in% samp_info$plot_date)

soil_info_df <- merge(relevant_soil_data, all_samp_reports, by="plot_date")

write.csv(soil_info_df, "/projectnb/talbot-lab-data/zrwerbin/soil_genome_db/example_soil_db_output.csv")

# Determine species that vary strongly with environment (just an example)
soil_info_df <- read.csv("/projectnb/talbot-lab-data/zrwerbin/soil_genome_db/example_soil_db_output.csv")

# Limit to species-level identifications and filter to species above .1%
species_abun <- soil_info_df[soil_info_df$taxRank=="S",] %>%
	filter(percentage > .01)

#make the dataframe wide so that correlations are easier
species_wide <- species_abun %>%
	select("plot_date", "siteID", "plotID",
				 "soilInCaClpH", "litterDepth", "soilTemp", "soilMoisture", "organicCPercent",
				 "nitrogenPercent", "CNratio", "percentage", "taxID", "name", "sampleID.y") %>%
	pivot_wider(names_from = name, values_from = percentage)

# Compute correlations
cor_fun <- function(df) cor.test(df$percentage, df$soilInCaClpH, method = "spearman", exact = FALSE) %>% tidy()
species_nest = species_abun %>% group_by(name, taxID) %>% nest()
data_nest <- mutate(species_nest, model = map(data, cor_fun))
# Filter by correlation strength
corr_pr <- select(data_nest, -data) %>% unnest(cols = c(model)) %>% filter(p.value <0.05 & abs(estimate) > .4)
species_of_interest <- corr_pr$name

# Visualize abundances that correlate with pH
ggplot(species_abun %>% dplyr::filter(name %in% species_of_interest[1:4])) +
	geom_point(aes(x = soilInCaClpH, y = percentage, color=name), alpha=.5, position=position_jitter(), size=2) +
	facet_wrap(~name, scales = "free") + theme_bw()

