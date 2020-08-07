# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(magrittr)
library(stringr)
library(readr)
library(httr)
library(seqinr)
library(Bethany)


load("intermediate_data/chromosome_features.Rdata")
load("intermediate_data/human_stress_granules.Rdata")

h_seq <- human_stress_granules %>%
	plyr::alply(1, function(df){
		cat("Getting fasta for gene '", df$gene_name, "' ...", sep="")
		r <- httr::GET(
			url="http://www.uniprot.org/uniprot/",
				httr::user_agent("httr mattjomeara@gmail.com"),
				query=list(
					query=paste0("organism:9606 and gene:", df$gene_name[[1]]),
					format='fasta',
					compress='yes')) %>%
			httr::content()
		if(r %>% is.null){
			cat(" MISSING\n")
			return(NULL)
		}
		cat(" GOT IT\n")
		seq_con <- r %>% rawConnection()
			sequences <- seq_con %>%
				gzcon() %>%
				seqinr::read.fasta(seqtype='AA')
		seq_con %>% close()

		cat("  n_sequences: ", length(sequences), "\n", sep="")

		is_swiss_prot <- sequences %>% seqinr::getName() %>% stringr::str_detect("^sp") %>% sum
		cat("  n_swiss_prot: ", is_swiss_prot %>% sum, "\n", sep="")
		sequences <- sequences[is_swiss_prot]
	})

# what proteins were we unable to find sequences?
z <- h_seq %>% purrr::map(~length(.)) %>% unlist
missing <- human_stress_granules %>% dplyr::slice(which(!z))

h_seq[human_stress_granules$gene_name == 'ALG13'] <- h_seq2[missing$gene_name == "ALG13"]
h_seq[human_stress_granules$gene_name == 'CNOT9'] <- h_seq2[missing$gene_name == "CNOT9"]
h_seq[human_stress_granules$gene_name == 'DNAJA1'] <- h_seq2[missing$gene_name == "DNAJA1"]
h_seq[human_stress_granules$gene_name == 'DUSP12'] <- h_seq2[missing$gene_name == "DUSP12"]
h_seq[human_stress_granules$gene_name == 'L1RE1'] <- h_seq2[missing$gene_name == "L1RE1"]
h_seq[human_stress_granules$gene_name == 'MKRN2'] <- h_seq2[missing$gene_name == "MKRN2"]

human_stress_granule_sequences <- human_stress_granules %>%
	dplyr::mutate(
		fasta = h_seq) %>%
			dplyr::filter(fasta %>% purrr::map_lgl(~!is.null(.))) %>%
	dplyr::mutate(
		seq_label = fasta %>% seqinr::getName()) %>%
	tidyr::separate(
		col="seq_label",
		into=c("database", "uniprot_accn", "uniprot_entry"),
		sep="[|]")
save(human_stress_granule_sequences, file="intermediate_data/human_stress_granule_sequences.Rdata")

seqinr::write.fasta(
	sequences=human_stress_granule_sequences$fasta,
	names=human_stress_granule_sequences$uniprot_accn,
	file.out="intermediate_data/human_stress_granule_sequences.fasta")

# read it back into to check that it looks ok
fasta <- seqinr::read.fasta(
	file="intermediate_data/human_stress_granule_sequences.fasta",
	seqtype="AA")

ca_to_hu_stress <- Bethany::blastp(
	ref="raw_data/C_albicans_SC5314_A22_current_default_protein.fasta",
	query="intermediate_data/human_stress_granule_sequences.fasta",
	run_id="ca-vs-hu_stress",
	verbose=TRUE) %>%
	dplyr::select(
		feature_name = ref_target,
		hu_uniprot_accn = query_target,
		bit_score,
		EValue) %>%
	dplyr::left_join(
		chromosome_features %>%
			dplyr::filter(feature_class=="ORF"),
		by="feature_name")


blastp_results_fname <- "/scratch/momeara/RtmpM31adH/file3ad7472f309a_blastp_ca-vs-hu_stress/blast_results.csv"
ca_to_hu_stress <- readr::read_tsv(
	file = blastp_results_fname,
	col_names = c(
		"query_target",
		"ref_target",
		"bit_score",
		"EValue"),
	col_types = readr::cols(
		query_target = readr::col_character(),
		ref_target = readr::col_character(),
		bit_score = readr::col_double(),
		EValue = readr::col_double())) %>%
	dplyr::select(
		feature_name = ref_target,
		hu_uniprot_accn = query_target,
		bit_score,
		EValue) %>%
	dplyr::inner_join(
		chromosome_features %>%
		dplyr::select(
			feature_name,
			gene_name,
			sac_ortholog,
			description,
			feature_type,
			feature_class,
			feature_status),
		by="feature_name") %>%
	dplyr::inner_join(
		human_stress_granule_sequences %>%
			dplyr::select(
				hu_uniprot_accn = uniprot_accn,
				hu_uniprot_entry = uniprot_entry,
				hu_gene_name = gene_name,
				hu_description = description,
				hu_in_p_bodies = in_p_bodies),
		by="hu_uniprot_accn")
save(ca_to_hu_stress, file="intermediate_data/ca_to_hu_stress.Rdata")

ca_to_hu_stress_best <- ca_to_hu_stress %>%
	dplyr::filter(EValue < 1e-80) %>%
	dplyr::arrange(EValue) %>%
	dplyr::group_by(feature_name) %>%
	dplyr::slice(1) %>%
	dplyr::ungroup()
readr::write_tsv(ca_to_hu_stress_best,"product/ca_to_hu_stress_best.tsv")


