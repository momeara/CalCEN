# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(stringr)
library(readr)
library(tidyr)
library(magrittr)
library(readxl)

load("intermediate_data/apms_tap_cytoscape.Rdata")
load("intermediate_data/apms_gfp.Rdata")
load("intermediate_data/cgd_id_2_uniprot_accn.Rdata")

source("scripts/annotate.R")


apms_tap_cytoscape <- apms_tap_cytoscape %>%
	dplyr::filter(
		!is.na(bait_feature_name),
		!is.na(prey_feature_name)) %>%
	dplyr::mutate(
		feature_name_1 = bait_feature_name,
		feature_name_2 = prey_feature_name) %>%
	annotate_associations() %>%
	dplyr::filter(
		bait_feature_name != prey_feature_name) %>%
	dplyr::filter(BFDR <= .02) %>%
	dplyr::mutate(
		is_coexp = coexp_score > .99,
		has_genetic_interaction = sac_gi_score > .95,
		is_novel = !is.na(sac_ppi_interactions) | !is.na(ca_biogrid_ppi),
		is_coannotated = !is.na(n_shared_go_terms))

nodes <- rbind(
	apms_tap_cytoscape %>%
		dplyr::select(
			feature_name = bait_feature_name,
			gene_name = bait_gene,
			sac_ortholog = sac_ortholog_1,
			description = description_1) %>%
		dplyr::mutate(
			node_type="bait"),
	apms_tap_cytoscape %>%
		dplyr::select(
			feature_name = prey_feature_name,
			gene_name = prey_gene,
			sac_ortholog = sac_ortholog_2,
			description = description_2) %>%
		dplyr::mutate(
			node_type="prey")) %>%
	dplyr::distinct(feature_name, .keep_all=TRUE) %>%
	dplyr::left_join(
		cgd_id_2_uniprot_accn %>%
			dplyr::filter(!is.na(feature_name)) %>%
			dplyr::select(
				feature_name,
				uniprot_accn),
		by="feature_name")

nodes %>%
	readr::write_tsv("product/apms_tap_cytoscape_nodes_181010.tsv")

edges <- apms_tap_cytoscape %>%
	dplyr::mutate(
		type="apms_tap") %>%
	dplyr::select(
		type,
		bait_feature_name,
		prey_feature_name,
		AvgSpec,
		BFDR,
		coexp_score,
		blastp_EValue,
		ca_biogrid_ppi,
		sac_ppi_interactions,
		sac_gi_score,
		n_shared_go_terms,
		is_coexp,
		has_genetic_interaction,
		is_novel,
		is_coannotated)

edges %>%
	readr::write_tsv("product/apms_tap_cytoscape_edges_181010.tsv")



######

# N unique proteins:
nodes %>% nrow
# 195

# N significant interactions
edges %>%
	dplyr::anti_join(
		edges, by=c(
			"bait_feature_name" = "prey_feature_name",
			"prey_feature_name" = "bait_feature_name")) %>%
	dplyr::distinct(
		bait_feature_name, prey_feature_name) %>%
	nrow
# 352

# N already annotated interactions
edges %>%
	dplyr::anti_join(
		edges, by=c(
			"bait_feature_name" = "prey_feature_name",
			"prey_feature_name" = "bait_feature_name")) %>%
	dplyr::filter(
		!is.na(ca_biogrid_ppi) | !is.na(sac_ppi_interactions)) %>%
	dplyr::distinct(
		bait_feature_name, prey_feature_name) %>%
	nrow
# 75



supp_table_1 <- edges %>%
	dplyr::left_join(
		nodes %>%
			dplyr::transmute(
				bait_feature_name=feature_name,
				bait_gene_name=gene_name,
				bait_sac_ortholog=sac_ortholog,
				bait_uniprot_accn=uniprot_accn),
		by="bait_feature_name") %>%
	dplyr::left_join(
		nodes %>%
			dplyr::transmute(
				prey_feature_name=feature_name,
				prey_gene_name=gene_name,
				prey_sac_ortholog=sac_ortholog,
				prey_uniprot_accn=uniprot_accn),
		by="prey_feature_name") %>%
	dplyr::select(
		bait_feature_name,
		bait_gene_name,
		bait_sac_ortholog,
		bait_uniprot_accn,
		prey_feature_name,
		prey_gene_name,
		prey_sac_ortholog,
		prey_uniprot_accn,
		AvgSpec,
		BFDR,
		blastp_EValue,
		sac_ppi_interactions,
		is_coannotated)

supp_table_1 %>% readr::write_tsv("product/supp_table_1.tsv")


apms_by_cond_avg_spec <- apms_gfp %>%
	dplyr::select(
		Prey, prey_gene, prey_feature_name, condition, AvgSpec) %>%
	tidyr::spread(key=condition, value=AvgSpec) %>%
	dplyr::select(
		Prey, prey_gene, prey_feature_name,
		NoDrug_AvgSpec=NoDrug, Caspo_AvgSpec=Caspo, Flu_AvgSpec=Flu) %>%

apms_by_cond_bfdr <- apms_gfp %>%
	dplyr::select(
		Prey, prey_gene, prey_feature_name, condition, BFDR) %>%
	tidyr::spread(key=condition, value=BFDR) %>%
	dplyr::select(
		Prey, prey_gene, prey_feature_name,
		NoDrug_BFDR=NoDrug, Caspo_BFDR=Caspo, Flu_BFDR=Flu)

apms_by_cond <- apms_by_cond_bfdr %>%
	dplyr::full_join(
		apms_by_cond_avg_spec,
		by=c("Prey", "prey_gene", "prey_feature_name")) %>%
	dplyr::mutate(
		Caspo_fc = Caspo_AvgSpec/NoDrug_AvgSpec, Flu_fc = Flu_AvgSpec/NoDrug_AvgSpec)



apms_by_cond %>%
	readr::write_tsv("product/apms_gfp_AvgSpec_fc_181110.tsv")





