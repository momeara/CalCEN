# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(tidyr)
library(tibble)
library(EGAD)

source("scripts/evaluate_network_gba.R")


#load("product/CalCEN_network_180612.Rdata")
#load("product/CalCEN_network_full_180704.Rdata")
load("product/CalCEN_network_full_20201007.Rdata")
load("product/CalCEN_intra_study_coexp_network_spearman_20201216.Rdata")

load("intermediate_data/ca_blastp_network.Rdata")
load("intermediate_data/ca_blastp_rank_network.Rdata")
load("intermediate_data/ca_ppi_network.Rdata")
load("intermediate_data/ca_silac_network.Rdata")

load("intermediate_data/yeast_net_network.Rdata")
load("intermediate_data/yeast_net_CC_network.Rdata")
load("intermediate_data/yeast_net_CX_network.Rdata")
load("intermediate_data/yeast_net_DC_network.Rdata")
load("intermediate_data/yeast_net_GN_network.Rdata")
load("intermediate_data/yeast_net_GT_network.Rdata")
load("intermediate_data/yeast_net_HT_network.Rdata")
load("intermediate_data/yeast_net_LC_network.Rdata")
load("intermediate_data/yeast_net_PG_network.Rdata")
load("intermediate_data/yeast_net_TS_network.Rdata")

load("intermediate_data/CoCoCoNet_candida_20210104.Rdata")
load("intermediate_data/CoCoCoNet_meta_candida_20210104.Rdata")

load("intermediate_data/ca_sac_ortholog_ppi_network.Rdata")
load("intermediate_data/ca_sac_ortholog_genetic_ppi_network.Rdata")
load("intermediate_data/ca_sac_ortholog_physical_ppi_network.Rdata")

load("intermediate_data/ca_genes.Rdata")
load("intermediate_data/chromosome_features.Rdata")
load("intermediate_data/ca_go_annotations.Rdata")
load("intermediate_data/ca_go_annotations_by_subontology.Rdata")
load("intermediate_data/ca_go_annotations_by_evidence.Rdata")

# Walk through of run_GBA to refresh my memory about how the cross-validation works
# EGAD::run_GBA(CalCEN_network_full, ca_go_annotations)
#    after filtering: n_genes = 4480 n_terms = 361
#    neighbor_voting([n_genes,n_genes], [n_genes,n_terms])
#       l = 361, g = 4480, n = 26,448
#       test.genes.labels = [n_genes, 3*n_terms]
#       for each term j
#           clear the annotation 1/nFold annotations from each fold in test.genes.labels


go_pred <- list(
	CalCEN_full_go_pred = EGAD::run_GBA(CalCEN_network_full, ca_go_annotations),
	ca_blastp_go_pred = EGAD::run_GBA(ca_blastp_network, ca_go_annotations),
	ca_ppi_go_pred = EGAD::run_GBA(ca_ppi_network, ca_go_annotations),
	ca_silac_go_pred = EGAD::run_GBA(ca_silac_network, ca_go_annotations),
	# ca_silac_pharmacological_unlabelled_go_pred = EGAD::run_GBA(silac_pharmacological_unlabelled_network, ca_go_annotations),
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



best_combo_ablation <- evaluate_network_gba(
	genes=ca_genes,
	annotation_sets=ca_go_annotations_by_subontology,
	networks=list(
		CalCEN = CalCEN_network_full,
		BlastP = ca_blastp_network,
		SacPhys = ca_sac_ortholog_physical_ppi_network,
		SacGene = ca_sac_ortholog_genetic_ppi_network,
		SILAC = ca_silac_network),
	degrees=c(3,4))
save(best_combo_ablation, file="intermediate_data/gba_best_combo_ablation.Rdata")


gba_summary <- evaluate_network_gba(
	genes=ca_genes,
	annotation_sets=ca_go_annotations_by_subontology,
	networks = list(
		CalCEN = CalCEN_network_full,
		BlastP = log(ca_blastp_network + 1),
		SacPhys = ca_sac_ortholog_physical_ppi_network,
		SacGene = ca_sac_ortholog_genetic_ppi_network))
save(gba_summary, file="intermediate_data/gba_summary_full.Rdata")



gba_summary <- evaluate_network_gba(
	genes=ca_genes,
	annotation_sets=ca_go_annotations_by_subontology,
	networks = list(
		CalCENp = CalCEN_network_full,
		CalCENi = CalCEN_intra_study_coexp_network_spearman,
		BlastP = ca_blastp_rank_network,
		SacPhys = ca_sac_ortholog_physical_ppi_network,
		SacGene = ca_sac_ortholog_genetic_ppi_network,
		YeastNet = yeast_net_network),
	nfold=10)
gba_summary %>%
		readr::write_tsv("product/gba_summary_10f_Cp-Ci-Br-SP-SG-YN_20201216.tsv")

source("scripts/evaluate_network_gba.R")
gba_summary <- evaluate_network_gba(
	genes=ca_genes,
	annotation_sets=ca_go_annotations_by_evidence,
	networks = list(
		CalCEN = CalCEN_network_full,
		BlastP = ca_blastp_network,
		SacPhys = ca_sac_ortholog_physical_ppi_network,
		SacGene = ca_sac_ortholog_genetic_ppi_network,
		YeastNet = yeast_net_network),
	nfold=10)
gba_summary %>%
		readr::write_tsv("product/gba_summary_byE_10f_C-B-SP-SG-YN_20201113.tsv")


gba_summary <- evaluate_network_gba(
	genes=ca_genes,
	annotation_sets=ca_go_annotations_by_subontology,
	networks = list(
		CalCENi = CalCEN_intra_study_coexp_network_spearman,
		CoCoCoNet = CoCoCoNet_candida,
		CoCoCoNet_Meta = CoCoCoNet_meta_candida),
	degrees = c(1),
  nfold = 10)


load("intermediate_data/CoCoCoNet_candida_prio_agg_net_20210104.Rdata")
load("intermediate_data/CoCoCoNet_candida_meta_agg_net_20210104.Rdata")

common_genes <- data.frame(gene = ca_genes) %>%
		dplyr::filter(gene %in% rownames(candida_prio_agg_net)) %>%
		dplyr::filter(gene %in% rownames(candida_meta_agg_net))


gba_summary_common <- evaluate_network_gba(
	genes=common_genes$gene,
	annotation_sets=ca_go_annotations_by_subontology,
	networks = list(
		CalCEN = CalCEN_intra_study_coexp_network_spearman[common_genes$gene, common_genes$gene],
		CoCoCoNet_Prio = CoCoCoNet_candida[common_genes$gene, common_genes$gene],
		CoCoCoNet_Meta = CoCoCoNet_meta_candida[common_genes$gene, common_genes$gene]),
	degrees = c(1),
  nfold = 10)
