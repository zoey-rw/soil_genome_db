#custom pavian functions to work with struo kraken output


read_report3 <- function(myfile,collapse=TRUE,keep_taxRanks=c("D","K","P","C","O","F","G","S"),min.depth=0,filter_taxon=NULL,
												 has_header=NULL,add_taxRank_columns=FALSE) {

	first.line <- readLines(myfile,n=1)
	isASCII <-  function(txt) all(charToRaw(txt) <= as.raw(127))
	if (!isASCII(first.line)) {
		dmessage(myfile," is no valid report - not all characters are ASCII")
		return(NULL)
	}
	if (is.null(has_header)) {
		has_header <- grepl("^[a-zA-Z]",first.line)
	}

	if (has_header) {
		report <- utils::read.table(myfile,sep="\t",header = T,
																quote = "",stringsAsFactors=FALSE, comment.char="#")
		#colnames(report) <- c("percentage","cladeReads","taxonReads","taxRank","taxID","n_unique_kmers","n_kmers","perc_uniq_kmers","name")

		## harmonize column names. TODO: Harmonize them in the scripts!
		colnames(report)[colnames(report)=="clade_perc"] <- "percentage"
		colnames(report)[colnames(report)=="perc"] <- "percentage"

		colnames(report)[colnames(report)=="n_reads_clade"] <- "cladeReads"
		colnames(report)[colnames(report)=="n.clade"] <- "cladeReads"

		colnames(report)[colnames(report)=="n_reads_taxo"] <- "taxonReads"
		colnames(report)[colnames(report)=="n.stay"] <- "taxonReads"

		colnames(report)[colnames(report)=="rank"] <- "taxRank"
		colnames(report)[colnames(report)=="tax_rank"] <- "taxRank"

		colnames(report)[colnames(report)=="taxonid"] <- "taxID"
		colnames(report)[colnames(report)=="tax"] <- "taxID"

	} else {
		report <- utils::read.table(myfile,sep="\t",header = F,
																col.names = c("percentage","cladeReads","taxonReads","taxRank","taxID","name"),
																quote = "",stringsAsFactors=FALSE, comment.char="#")
	}

	report$depth <- nchar(gsub("\\S.*","",report$name))/2
	report$name <- gsub("^ *","",report$name)
	report$name <- paste(tolower(report$taxRank),report$name,sep="_")


	## Only stop at certain taxRanks
	## filter taxon and further up the tree if 'filter_taxon' is defined
	kraken.tree <- pavian:::build_kraken_tree(report)
	#report <- pavian:::collapse.taxRanks(kraken.tree,keep_taxRanks=keep_taxRanks,filter_taxon=filter_taxon)

	## Add a metaphlan-style taxon string
	# if (add_taxRank_columns) {
	# 	report[,keep_taxRanks] <- NA
	# }
	report$taxLineage = report$name
	rows_to_consider <- rep(FALSE,nrow(report))

	report$percentage <- round(report$cladeReads/sum(report$taxonReads),6) * 100

	for (column in c("taxonReads", "cladeReads"))
		if (all(floor(report[[column]]) == report[[column]]))
			report[[column]] <- as.integer(report[[column]])

	if ('n_unique_kmers'  %in% colnames(report))
		report$kmerpercentage <- round(report$n_unique_kmers/sum(report$n_unique_kmers,na.rm=T),6) * 100
	#report$taxRankperc <- 100/taxRank(report$cladeReads)

	rownames(report) <- NULL

	report
}

summarize_report_custom <- function(my_report) {
	protist_taxids <- c("-_Diplomonadida"=5738,
											"-_Amoebozoa"=554915,
											"-_Alveolata"=33630)
	rownames(my_report)[rownames(my_report) == "r_root"] <- "-_root"
	my_report <- my_report[!duplicated(my_report$name),]
	row.names(my_report) <- my_report[["name"]]
	unidentified_reads <- my_report["u_unclassified","cladeReads"]
	identified_reads <- my_report["r_root","cladeReads"]
	artificial_reads <- pavian:::zero_if_na1(my_report["s_synthetic construct",
																										 "cladeReads"])
	human_reads <- pavian:::zero_if_na1(my_report["s_Homo sapiens", "cladeReads"])
	chordate_reads <- pavian:::zero_if_na1(my_report["p_Chordata", "cladeReads"])
	root_reads <- pavian:::zero_if_na1(my_report["r_root", "taxonReads"])
	out <- data.frame(number_of_raw_reads = unidentified_reads + identified_reads,
										classified_reads = identified_reads, chordate_reads = chordate_reads,
										artificial_reads = artificial_reads, unclassified_reads = unidentified_reads,
										microbial_reads = identified_reads - chordate_reads -
											artificial_reads - root_reads, bacterial_reads = pavian:::zero_if_na1(my_report["d_Bacteria",
																																																			"cladeReads"]) + pavian:::zero_if_na1(my_report["k_Bacteria",
																																																																											"cladeReads"]), viral_reads = pavian:::zero_if_na1(my_report["d_Viruses",
																																																																																																									 "cladeReads"]) + pavian:::zero_if_na1(my_report["k_Viruses",
																																																																																																									 																								"cladeReads"]), fungal_reads = pavian:::zero_if_na1(my_report["k_Fungi",
																																																																																																									 																																																							"cladeReads"]), protozoan_reads = sum(pavian:::zero_if_na1(my_report[names(protist_taxids),
																																																																																																									 																																																																																									 "cladeReads"])))
}
