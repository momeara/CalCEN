# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(tidyr)
library(tibble)
library(EGAD)

source("scripts/evaluate_network_gba.R")
source

load("product/ca_coexp_network_180612.Rdata")
load("product/ca_coexp_network_full_180704.Rdata")
load("intermediate_data/ca_blastp_network.Rdata")
#load("intermediate_data/ca_ppi_network.Rdata")
load("intermediate_data/ca_silac_network.Rdata")


#load("intermediate_data/ca_sac_ortholog_ppi_network.Rdata")
load("intermediate_data/ca_sac_ortholog_genetic_ppi_network.Rdata")
load("intermediate_data/ca_sac_ortholog_physical_ppi_network.Rdata")

load("intermediate_data/ca_genes.Rdata")
load("intermediate_data/chromosome_features.Rdata")
load("intermediate_data/ca_go_annotations.Rdata")
load("intermediate_data/ca_go_annotations_by_subontology.Rdata")
load("intermediate_data/ca_go_annotations_by_evidence.Rdata")


go_pred <- list(
	ca_coexp_full_go_pred = EGAD::run_GBA(ca_coexp_network_full, ca_go_annotations),
	ca_blastp_go_pred = EGAD::run_GBA(ca_blastp_network, ca_go_annotations),
	ca_ppi_go_pred = EGAD::run_GBA(ca_ppi_network, ca_go_annotations),
	ca_silac_go_pred = EGAD::run_GBA(ca_silac_network, ca_go_annotations),
	ca_silac_pharmacological_unlabelled_go_pred = EGAD::run_GBA(silac_pharmacological_unlabelled_network, ca_go_annotations),
	ca_sac_ortholog_ppi_go_pred = EGAD::run_GBA(ca_sac_ortholog_ppi_network, ca_go_annotations),
	ca_sac_ortholog_genetic_ppi_go_pred = EGAD::run_GBA(ca_sac_ortholog_genetic_ppi_network, ca_go_annotations),
	ca_sac_ortholog_physical_ppi_go_pred = EGAD::run_GBA(ca_sac_ortholog_physical_ppi_network, ca_go_annotations),
	ca_yeast_net_go_pred = EGAD::run_GBA(yeast_net_network, ca_go_annotations),
	ca_yeast_net_CC_go_pred = EGAD::run_GBA(yeast_net_CC_network, ca_go_annotations),
	ca_yeast_net_CX_go_pred = EGAD::run_GBA(yeast_net_CX_network, ca_go_annotations),
	ca_yeast_net_DC_go_pred = EGAD::run_GBA(yeast_net_DC_network, ca_go_annotations),
	ca_yeast_net_GN_go_pred = EGAD::run_GBA(yeast_net_GN_network, ca_go_annotations),
	ca_yeast_net_GT_go_pred = EGAD::run_GBA(yeast_net_GT_network, ca_go_annotations),
	ca_yeast_net_HT_go_pred = EGAD::run_GBA(yeast_net_HT_network, ca_go_annotations),
	ca_yeast_net_LC_go_pred = EGAD::run_GBA(yeast_net_LC_network, ca_go_annotations),
	ca_yeast_net_PG_go_pred = EGAD::run_GBA(yeast_net_PG_network, ca_go_annotations),
	ca_yeast_net_TS_go_pred = EGAD::run_GBA(yeast_net_TS_network, ca_go_annotations))
save(go_pred, file="intermediate_data/go_pred.Rdata")


best_combo_abalation <- evaluate_network_gba(
	genes=ca_genes,
	annotation_sets=ca_go_annotations_by_subontology,
	networks=list(
		`Co-Exp` = ca_coexp_network_full,
#		BlastP = ca_blastp_network,
		SacPhys = ca_sac_ortholog_physical_ppi_network,
		SacGene = ca_sac_ortholog_genetic_ppi_network,
		SILAC = ca_silac_network),
	degrees=c(3,4))



gba_summary <- evaluate_network_gba(
	genes=ca_genes,
	annotation_sets=ca_go_annotations_by_subontology,
	networks = list(
		`Co-Exp` = ca_coexp_network_full,
		BlastP = ca_blastp_network,
		SacPhys = ca_sac_ortholog_physical_ppi_network,
		SacGene = ca_sac_ortholog_genetic_ppi_network))
save(gba_summary_full, file="intermediate_data/gba_summary_full.Rdata")



gba_summary <- evaluate_network_gba(
	genes=ca_genes,
	annotation_sets=ca_go_annotations_by_subontology,
	networks = list(
		`Co-Exp` = ca_coexp_network_full,
		BlastP = ca_blastp_network,
		SacPhys = ca_sac_ortholog_physical_ppi_network,
		SacGene = ca_sac_ortholog_genetic_ppi_network,
		YeastNet = yeast_net_network),
	nfold=10)
gba_summary %>% readr::write_tsv("product/gba_summary_10f_C-B-SP-SG-YN_180706.tsv")

gba_summary <- evaluate_network_gba(
	genes=ca_genes,
	annotation_sets=ca_go_annotations_by_evidence,
	networks = list(
		`Co-Exp` = ca_coexp_network_full,
		BlastP = ca_blastp_network,
		SacPhys = ca_sac_ortholog_physical_ppi_network,
		SacGene = ca_sac_ortholog_genetic_ppi_network,
		YeastNet = yeast_net_network),
	nfold=10)
gba_summary %>% readr::write_tsv("product/gba_summary_byE_10f_C-B-SP-SG-YN_180706.tsv")

