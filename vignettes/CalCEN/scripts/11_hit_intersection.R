# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(stringr)
library(readr)
library(tidyr)
library(magrittr)

source("scripts/annotate.R")

load("intermediate_data/apms_tap.Rdata")
load("intermediate_data/apms_gfp.Rdata")
load("intermediate_data/genes_of_interest.Rdata")

apms_tap_annotation <- apms_tap %>%
	dplyr::rename(
		feature_name_1 = bait_feature_name,
		feature_name_2 = prey_feature_name) %>%
	annotate_associations()

apms_tap_annotation %>%
	readr::write_tsv("product/apms_tap_annotation_180718.tsv")

####

apms_gfp_annotation <- apms_gfp %>%
	dplyr::mutate(
		feature_name_1 = bait_feature_name,
		feature_name_2 = prey_feature_name) %>%
	annotate_associations()

apms_gfp_annotation %>%
	readr::write_tsv("product/apms_gfp_annotation_180717.tsv")


####

genes_of_interest_annotation <-
	genes_of_interest %>%
	dplyr::select(-gene_name, -sac_ortholog, -description) %>%
	annotate_genes()

genes_of_interest_annotation %>%
	dplyr::filter(set != "Predicted Stress Granule") %>%
	readr::write_tsv("product/genes_of_interest_annotation_180718.tsv")

genes_of_interest_annotation %>%
	dplyr::filter(set == "R2TP Complex") %>%
	dplyr::arrange(set, gene_name, source_organism, mode, method, score) %>%
	dplyr::transmute(
		set,
		feature_name,
		gene_name,
		sac_ortholog,
		description = description %>% stringr::str_sub(1,50),
		feature_name_2,
		gene_name_2,
		sac_ortholog_2,
		description_2 = description_2 %>% stringr::str_sub(1,50),
		source_organism,
		mode,
		method,
		score) %>%
	data.frame
