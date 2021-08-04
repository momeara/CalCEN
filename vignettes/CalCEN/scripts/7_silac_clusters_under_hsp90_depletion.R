# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(stringr)
library(readr)
library(tidyr)
library(magrittr)


load("intermediate_data/ca_silac_clusters.Rdata")
load("intermediate_data/ca_silac_hsp90_fold_change.Rdata")


ca_silac_fc <- ca_silac_clusters %>%
	dplyr::inner_join(
		ca_silac_hsp90_fold_change,
		by=c("feature_name", "uniprot_accn"))


mean_cluster_fc <- ca_silac_fc %>%
	dplyr::group_by(
		complex_index, feature_name) %>%
	dplyr::summarize(
		mean_pharm=mean(pharmacological_log_fc),
		mean_gene =mean(genetic_log_fc)) %>%
	dplyr::summarize(
		cluster_size = n(),
		mean_cluster_pharm = mean(mean_pharm),
		mean_cluster_gene  = mean(mean_gene))



	dplyr::group_by(cluster_index) %>%
	dplyr::
