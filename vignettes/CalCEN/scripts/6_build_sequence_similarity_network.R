# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(EGAD)

load("intermediate_data/ca_blastp.Rdata")
load("intermediate_data/ca_genes.Rdata")

ca_blastp_network <- ca_blastp %>%
	dplyr::select(feature_name_1, feature_name_2, bit_score) %>%
	as.data.frame() %>%
	EGAD::build_weighted_network(ca_genes)

save(ca_blastp_network, file="intermediate_data/ca_blastp_network.Rdata")

ca_blastp_network %>%
	as.data.frame() %>%
	readr::write_tsv("product/ca_blastp_network_20201113.Rdata")


ca_blastp_rank_network <- ca_blastp %>%
		dplyr::mutate(
				rank_score = rank(bit_score, na.last = "keep", ties.method = "average") / dplyr::n()) %>%
		dplyr::select(feature_name_1, feature_name_2, rank_score) %>%
		as.data.frame() %>%
		EGAD::build_weighted_network(ca_genes)

ca_blastp_rank_network %>%
	as.data.frame() %>%
	readr::write_tsv("product/ca_blastp_rank_network_20201124.tsv")

save(ca_blastp_rank_network, file = "intermediate_data/ca_blastp_rank_network.Rdata")


# save it in the long format
ca_blastp_rank_network_full <- ca_blastp_rank_network %>%
		data.frame() %>%
		tibble::rownames_to_column("feature_name_1") %>%
		tidyr::gather(key="feature_name_2", value="score", -feature_name_1)

save(ca_blastp_rank_network_full, file = "intermediate_data/ca_blastp_rank_network_full.Rdata")

blastp_rank_network_full %>%
		as.data.frame() %>%
		readr::write_tsv("product/ca_blastp_rank_network_full_20201124.tsv")
