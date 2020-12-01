# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(stringr)
library(readr)
library(tidyr)
library(magrittr)

cat("Loading annotation datasets ... ")
load("intermediate_data/chromosome_features.Rdata")
load("intermediate_data/CalCEN_full.Rdata")
load("intermediate_data/ca_blastp.Rdata")
load("intermediate_data/ca_biogrid.Rdata")

load("intermediate_data/apms_tap.Rdata")
load("intermediate_data/apms_gfp.Rdata")
load("intermediate_data/ca_silac_hsp90_fold_change.Rdata")

load("intermediate_data/ca_ppi_network.Rdata")
load("intermediate_data/ca_sac_ortholog_physical_ppi_summary.Rdata")
load("intermediate_data/ca_sac_ortholog_physical_ppi.Rdata")
load("intermediate_data/ca_sac_ortholog_genetic_ppi.Rdata")
load("intermediate_data/co_annotation_count.Rdata")
cat(" DONE\n")

annotate_genes <- function(data){
	coexp_data <- data %>%
		dplyr::inner_join(
			CalCEN_full %>%
				dplyr::transmute(
					feature_name = feature_name_1,
					feature_name_2,
					source_organism="Ca",
					mode = "Co-Expression",
					method = "RNA-Seq",
					score) %>%
				dplyr::filter(
					feature_name != feature_name_2,
					score >= .995),
				by=c("feature_name"))

	blastp_data <- data %>%
		dplyr::inner_join(
			ca_blastp %>%
				dplyr::filter(
					feature_name_1 != feature_name_2,
					EValue < 1e-80) %>%
				dplyr::transmute(
					feature_name = feature_name_1,
					feature_name_2,
					source_organism="Ca",
					mode = "Sequence Similarity",
					method= "BlastP",
					score=EValue),
				by=c("feature_name"))

	apms_tap_data <- data %>%
		dplyr::inner_join(
			apms_tap %>%
				dplyr::filter(
					hit,
					!is.na(prey_feature_name)) %>%
				dplyr::transmute(
					feature_name = prey_feature_name,
					feature_name_2 = bait_feature_name,
					source_organism="Ca",
					mode="Physical PPI",
					method="AP-MS Bait-TAP",
					score=1),
			by=c("feature_name"))

	apms_gfp_data <- data %>%
		dplyr::inner_join(
			apms_gfp %>%
				dplyr::filter(hit) %>%
				dplyr::transmute(
					feature_name = prey_feature_name,
					feature_name_2 = bait_feature_name,
					source_organism="Ca",
					mode="Physical PPI",
					method="AP-MS Bait-GFP",
					score=1),
			by="feature_name")

	silac_data <- data %>%
		dplyr::inner_join(
			ca_silac_hsp90_fold_change %>%
				dplyr::filter(!is.na(feature_name)) %>%
				dplyr::group_by(feature_name) %>%
				dplyr::summarize(
					`Physical PPI HSP90-pharmacological` = sum(abs(pharmacological_log_fc)),
					`Physical PPI HSP90-genetic` = sum(abs(genetic_log_fc))) %>%
				tidyr::gather("mode","score", -feature_name) %>%
				dplyr::mutate(
					feature_name_2 = NA,
					source_organism="Ca",
					method="SILAC sum(abs(log(fc)))"),
			by="feature_name")

	ca_biogrid_data <- data %>%
		dplyr::inner_join(
			ca_biogrid %>%
				dplyr::filter(
					!((experimental_system_type == "genetic") & (feature_name_1 == feature_name_2))) %>%
				dplyr::transmute(
					feature_name = feature_name_1,
					feature_name_2,
					source_organism="Ca",
					mode=dplyr::case_when(
						experimental_system_type == "physical" ~ "Physical PPI",
						experimental_system_type == "genetic" ~ "Genetic PPI",
						TRUE ~ NA_character_),
					method = paste0(
						experimental_system_abbreviation, "|",
						throughput_abbreviation),
					score=1),
			by="feature_name")

	sac_biogrid_physical_data <- data %>%
		dplyr::inner_join(
			sac_biogrid %>%
				dplyr::filter(
					experimental_system_type == "physical",
					!is.na(ca_feature_name_1),
					!is.na(ca_feature_name_2)) %>%
				dplyr::select(
					feature_name = ca_feature_name_1,
					feature_name_2 = ca_feature_name_2) %>%
				dplyr::semi_join(data, by="feature_name") %>%
				dplyr::count(feature_name, feature_name_2) %>%
				dplyr::transmute(
					feature_name,
					feature_name_2,
					source_organism="Sc",
					mode="Physical PPI",
					method="BioGRID",
					score=n),
			by="feature_name")

	sac_biogrid_genetic_data <- data %>%
		dplyr::inner_join(
			ca_sac_ortholog_genetic_ppi %>%
				dplyr::transmute(
					feature_name = feature_name_1,
					feature_name_2,
					source_organism="Sc",
					mode="Genetic PPI",
					method="BioGRID",
					score = score) %>%
				dplyr::filter(
					feature_name != feature_name_2,
					score >= .999),
			by="feature_name")

	rbind(
		coexp_data,
		blastp_data,
		apms_tap_data,
		apms_gfp_data,
		silac_data,
		ca_biogrid_data,
		sac_biogrid_physical_data,
		sac_biogrid_genetic_data) %>%
		dplyr::left_join(
			chromosome_features %>%
				dplyr::select(
					feature_name = feature_name,
					gene_name = gene_name,
					sac_ortholog = sac_ortholog,
					description = description),
			by=c("feature_name")) %>%
		dplyr::left_join(
			chromosome_features %>%
				dplyr::select(
					feature_name_2 = feature_name,
					gene_name_2 = gene_name,
					sac_ortholog_2 = sac_ortholog,
					description_2 = description),
			by=c("feature_name_2"))

}



annotate_associations <- function(data){
	data %T>%

		(~cat("Add identifiers\n")) %>%
		dplyr::left_join(
			chromosome_features %>%
				dplyr::select(
					feature_name_1 = feature_name,
					gene_name_1 = gene_name,
					sac_ortholog_1 = sac_ortholog,
					description_1 = description),
			by=c("feature_name_1")) %>%
		dplyr::left_join(
			chromosome_features %>%
				dplyr::select(
					feature_name_2 = feature_name,
					gene_name_2 = gene_name,
					sac_ortholog_2 = sac_ortholog,
					description_2 = description),
			by=c("feature_name_2")) %T>%

		(~cat("Add Co-expression scores\n")) %>%
		dplyr::left_join(
			CalCEN_full %>%
				dplyr::select(
					feature_name_1,
					feature_name_2,
					coexp_score = score),
				by=c("feature_name_1", "feature_name_2")) %T>%

		(~cat("Add BlastP sequence similarity scores\n")) %>%
		dplyr::left_join(
			ca_blastp %>%
				dplyr::select(
					feature_name_1,
					feature_name_2,
					blastp_EValue=EValue),
				by=c("feature_name_1", "feature_name_2")) %T>%

		(~cat("Add AP-MS PPI hits with HSP90+Co-Chaperones baits labeled with TAP\n")) %>%
		dplyr::left_join(
			apms_tap %>%
				dplyr::filter(hit) %>%
				dplyr::transmute(
					feature_name_1 = bait_feature_name,
					feature_name_2 = prey_feature_name,
					apms_tap_hit = TRUE),
			by=c("feature_name_1", "feature_name_2")) %T>%

		(~cat("Add AP-MS PPI hits with HSP90 labeled with GFP\n")) %>%
		dplyr::left_join(
			apms_gfp %>%
				dplyr::filter(hit) %>%
				dplyr::transmute(
					feature_name_1 = prey_feature_name,
					apms_hsp90_gfp_interactor_1 = TRUE),
			by=c("feature_name_1")) %T>%
		dplyr::left_join(
			apms_gfp %>%
				dplyr::filter(hit) %>%
				dplyr::transmute(
					feature_name_2 = prey_feature_name,
					apms_hsp90_gfp_interactor_2 = TRUE),
			by=c("feature_name_2")) %T>%

		(~cat("Add SILAC sum(abs(log(fc))) with HSP90 depletion\n")) %>%
		dplyr::left_join(
			ca_silac_hsp90_fold_change %>%
				dplyr::rename(
					feature_name_1 = feature_name) %>%
				dplyr::group_by(feature_name_1) %>%
				dplyr::summarize(
					silac_pharmacological_log_fc_1 = sum(abs(pharmacological_log_fc)),
					silac_genetic_log_fc_1 = sum(abs(genetic_log_fc))),
			by=c("feature_name_1")) %>%
		dplyr::left_join(
			ca_silac_hsp90_fold_change %>%
				dplyr::rename(
					feature_name_2 = feature_name) %>%
				dplyr::group_by(feature_name_2) %>%
				dplyr::summarize(
					silac_pharmacological_log_fc_2 = sum(abs(pharmacological_log_fc)),
					silac_genetic_log_fc_2 = sum(abs(genetic_log_fc))),
			by=c("feature_name_2")) %T>%

		(~cat("Add known Ca PPI interactions from BioGRID\n")) %>%
		dplyr::left_join(
			ca_biogrid %>%
				dplyr::semi_join(
					data,
					by=c("feature_name_1", "feature_name_2")) %>%
				dplyr::transmute(
					feature_name_1,
					feature_name_2,
					ca_biogrid_ppi = paste(
						experimental_system_type,
						experimental_system_abbreviation,
						throughput_abbreviation, sep="|")),
			by=c("feature_name_1", "feature_name_2")) %T>%

		(~cat("Add orthologous physical interactions in sac from BioGRID\n")) %>%
		dplyr::left_join(
			ca_sac_ortholog_physical_ppi_summary %>%
				dplyr::select(
					feature_name_1,
					feature_name_2,
					sac_ppi_interactions=interactions),
			by=c("feature_name_1", "feature_name_2")) %T>%

		(~cat("Add orthologous genetic interaction in sac score from BioGRID\n")) %>%
		dplyr::left_join(
			ca_sac_ortholog_genetic_ppi %>%
				dplyr::select(
					feature_name_1,
					feature_name_2,
					sac_gi_score = score),
			by=c("feature_name_1", "feature_name_2")) %T>%

		(~cat("Add number of shared go-term annotations\n")) %>%
		dplyr::left_join(
			co_annotation_count %>%
				dplyr::select(
					feature_name_1,
					feature_name_2,
					n_shared_go_terms = score),
			by=c("feature_name_1", "feature_name_2"))
}
