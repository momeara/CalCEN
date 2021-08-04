# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(reshape2)
library(tidyr)

source("scripts/embed_network.R")

load("intermediate_data/ca_go_propagated_filtered.Rdata")
load("intermediate_data/ca_genes.Rdata")

ca_go_positive <-  ca_go_propagated_filtered %>%
	dplyr::filter(is.na(qualifier) || qualifier != "NOT") %>%
	dplyr::distinct(feature_name, go_id, .keep_all=TRUE)

co_annotation <-ca_go_positive %>%
	dplyr::select(
		feature_name_1 = feature_name,
		go_id) %>%
	dplyr::inner_join(
		ca_go_positive %>%
			dplyr::select(
				feature_name_2 = feature_name,
				go_id),
		by=c("go_id"))

z <- co_annotation %>%
	dplyr::filter(feature_name_1 < feature_name_2) %>%
	dplyr::group_by(feature_name_1, feature_name_2) %>%
	dplyr::summarize(n=length(go_id)) %>%
	dplyr::ungroup()
save(co_annotation, file="intermediate_data/co_annotation.Rdata")


co_annotation_network <- co_annotation %>%
	reshape2::acast(
		feature_name_1 ~ feature_name_2,
		value.var="go_id",
		fun.aggregate=length) %>%
	embed_network(ca_genes)
save(co_annotation_network, file="intermediate_data/co_annotation_network.Rdata")

co_annotation_count <- co_annotation_network %>%
	data.frame() %>%
	tibble::rownames_to_column("feature_name_1") %>%
	tidyr::gather(key="feature_name_2", value="score", -feature_name_1) %>%
	dplyr::filter(score!=0)
save(co_annotation_count, file="intermediate_data/co_annotation_count.Rdata")


