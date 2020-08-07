# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(magrittr)
library(stringr)
library(readr)

library(Bethany)

load("intermediate_data/chromosome_features.Rdata")

system("cd raw_data && \\
  wget http://www.candidagenome.org/download/sequence/C_albicans_SC5314/Assembly22/current/C_albicans_SC5314_A22_current_default_protein.fasta.gz && \\
  gunzip C_albicans_SC5314_A22_current_default_protein.fasta.gz")


ca_blastp <- Bethany::blastp(
	ref="raw_data/C_albicans_SC5314_A22_current_default_protein.fasta",
	query="raw_data/C_albicans_SC5314_A22_current_default_protein.fasta",
	run_id="ca-vs-ca",
	verbose=TRUE) %>%
	dplyr::select(
		feature_name_1 = ref_target,
		feature_name_2 = query_target,
		bit_score,
		EValue) %>%
	dplyr::semi_join(
		chromosome_features %>%
			dplyr::filter(feature_class=="ORF"),
		by=c("feature_name_1"="feature_name")) %>%
	dplyr::semi_join(
		chromosome_features %>%
			dplyr::filter(feature_class=="ORF"),
		by=c("feature_name_2"="feature_name"))

save(ca_blastp, file="intermediate_data/ca_blastp.Rdata")
