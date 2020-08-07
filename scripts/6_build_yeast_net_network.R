# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(stringr)
library(readr)
library(reshape2)
library(EGAD)
library(glue)

load("intermediate_data/ca_genes.Rdata")
load("intermediate_data/chromosome_features.Rdata")
load("intermediate_data/sac_chromosome_features.Rdata")
load("intermediate_data/yeast_net.Rdata")
load("intermediate_data/yeast_net_CC.Rdata")
load("intermediate_data/yeast_net_CX.Rdata")
load("intermediate_data/yeast_net_DC.Rdata")
load("intermediate_data/yeast_net_GN.Rdata")
load("intermediate_data/yeast_net_GT.Rdata")
load("intermediate_data/yeast_net_HT.Rdata")
load("intermediate_data/yeast_net_LC.Rdata")
load("intermediate_data/yeast_net_PG.Rdata")
load("intermediate_data/yeast_net_TS.Rdata")


ca_to_sac <- chromosome_features %>%
	dplyr::filter(
		feature_name %in% ca_genes,
		!is.na(sac_ortholog)) %>%
	dplyr::select(
		ca_feature_name=feature_name,
		sac_ortholog) %>%
	dplyr::inner_join(
		sac_chromosome_features %>%
			dplyr::transmute(
				sac_ortholog = ifelse(is.na(standard_gene_name), feature_name, standard_gene_name),
				sac_feature_name = feature_name),
		by="sac_ortholog") %>%
	dplyr::select(
		ca_feature_name,
		sac_feature_name)

build_yeast_net_network <- function(tag, data){

	cat("Building network for '", tag, "' network\n", sep="")
 	ca_data <- data %>%
		dplyr::inner_join(
			ca_to_sac %>% dplyr::rename(gene1=sac_feature_name, feature_name_1=ca_feature_name),
			by="gene1") %>%
		dplyr::inner_join(
			ca_to_sac %>% dplyr::rename(gene2=sac_feature_name, feature_name_2=ca_feature_name),
			by="gene2") %>%
		dplyr::select(feature_name_1, feature_name_2, score) %>%
		dplyr::distinct(feature_name_1, feature_name_2, .keep_all=TRUE)

	n_nodes <- c(
		ca_data %>% magrittr::extract2("feature_name_1"),
		ca_data %>% magrittr::extract2("feature_name_2")) %>%
		unique %>%
		length
	cat("\tn nodes: ", n_nodes , "\n", sep="")
	cat("\tn edges: ", ca_data %>% nrow, "\n", sep="")

	network <- ca_data %>%
		as.data.frame() %>%
		EGAD::build_weighted_network(ca_genes)
}


yeast_net_network <- build_yeast_net_network("full", yeast_net)
yeast_net_CC_network <- build_yeast_net_network("CC", yeast_net_CC)
yeast_net_CX_network <- build_yeast_net_network("CX", yeast_net_CX)
yeast_net_DC_network <- build_yeast_net_network("DC", yeast_net_DC)
yeast_net_GN_network <- build_yeast_net_network("GN", yeast_net_GN)
yeast_net_GT_network <- build_yeast_net_network("GT", yeast_net_GT)
yeast_net_HT_network <- build_yeast_net_network("HT", yeast_net_HT)
yeast_net_LC_network <- build_yeast_net_network("LC", yeast_net_LC)
yeast_net_PG_network <- build_yeast_net_network("PG", yeast_net_PG)
yeast_net_TS_network <- build_yeast_net_network("TS", yeast_net_TS)

save(yeast_net_network, file="intermediate_data/yeast_net_network.Rdata")
save(yeast_net_CC_network, file="intermediate_data/yeast_net_CC_network.Rdata")
save(yeast_net_CX_network, file="intermediate_data/yeast_net_CX_network.Rdata")
save(yeast_net_DC_network, file="intermediate_data/yeast_net_DC_network.Rdata")
save(yeast_net_GN_network, file="intermediate_data/yeast_net_GN_network.Rdata")
save(yeast_net_GT_network, file="intermediate_data/yeast_net_GT_network.Rdata")
save(yeast_net_HT_network, file="intermediate_data/yeast_net_HT_network.Rdata")
save(yeast_net_LC_network, file="intermediate_data/yeast_net_LC_network.Rdata")
save(yeast_net_PG_network, file="intermediate_data/yeast_net_PG_network.Rdata")
save(yeast_net_TS_network, file="intermediate_data/yeast_net_TS_network.Rdata")


