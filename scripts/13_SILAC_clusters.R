# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(tidyr)
library(tibble)
library(EGAD)
library(apcluster)

load("intermediate_data/chromosome_features.Rdata")
load("intermediate_data/ca_silac_hsp90_intensities.Rdata")
load("intermediate_data/chromosome_features.Rdata")
load("intermediate_data/ca_silac_network.Rdata")
load("intermediate_data/apms_tap.Rdata")

" For each of the wildtype, and HSP90 pharmacological and genetic depletion conditions, the spearman rank correlation of intensities, averaging ties, was computed for each pair of genes using the EGAD R package (Ballouz, et al. 2017, doi: 10.1093/bioinformatics/btw695). Taking the correlations for each condition as a network, the genes were clustered using affinity propagation using the apcluster R package (Frey and Dueck, 2007, doi:10.1126/science.1136800)

"


silac_genes <- ca_silac_hsp90_intensities %>%
	dplyr::distinct(gene_label) %>%
	magrittr::extract2("gene_label")

find_silac_clusters <- function(input_condition){
	silac_network <- ca_silac_hsp90_intensities %>%
		dplyr::filter(condition==input_condition) %>%
		dplyr::select(gene_label, fraction_index, intensity) %>%
		tidyr::spread(fraction_index, intensity) %>%
		tibble::column_to_rownames("gene_label") %>%
		as.matrix() %>%
		EGAD::build_coexp_network(silac_genes)
	silac_network[is.na(silac_network)] <- 0

	silac_clusters <- silac_network %>%
		apcluster::apcluster()

	cluster_id <- 0
	silac_clusters@clusters %>%
		plyr::ldply(function(members){
			cluster_id <<- cluster_id + 1
			data.frame(
				cluster_id = cluster_id,
				exemplar = names(silac_clusters@exemplars)[cluster_id],
				gene_label = names(members))
		}) %>%
		dplyr::left_join(
			ca_silac_hsp90_intensities %>%
				dplyr::filter(condition==input_condition) %>%
				dplyr::distinct(gene_label) %>%
				dplyr::select(
					gene_label,
					feature_name,
					uniprot_entry,
					uniprot_accn,
					description),
			by="gene_label")
}

silac_clusters_wt <- find_silac_clusters("wildtype")
silac_clusters_ph <- find_silac_clusters("pharmacological")
silac_clusters_gt <- find_silac_clusters("genetic")

silac_clusters <- rbind(
	silac_clusters_wt,
	silac_clusters_ph,
	silac_clusters_gt)
save(silac_clusters, file="intermediate_data/silac_clusters.Rdata")

silac_clusters %>%
	dplyr::distinct(gene_label, condition, exemplar, .keep_all=TRUE) %>%
	dplyr::select(-fraction_index, -intensity) %>%
	readr::write_tsv("product/silac_clusters_180913.tsv")

silac_clusters_wt %>%
	dplyr::left_join(
		ca_silac_hsp90_intensities %>%
			dplyr::filter(condition=="wildtype") %>%
			dplyr::select(
				gene_label,
				fraction_index,
				intensity),
		by="gene_label") %>%
	tidyr::spread(fraction_index, intensity) %>%
	readr::write_tsv("product/silac_hsp90_depletion_180913.tsv")





#### Test if APMS HSP90 clients are enriched in HSP90 cluster
HSP90_cluster <- silac_clusters_long %>%
	dplyr::filter(gene_label == "HSP90") %>%
	magrittr::extract2("cluster_id")

HSP90_APMS_clients <- apms_tap %>%
	dplyr::filter(
		bait_gene == "HSP90",
		prey_gene != "HSP90",
		hit) %>%
	magrittr::extract2("prey_feature_name")

in_HSP90_cluster <- (silac_clusters_long$cluster_id == HSP90_cluster) & (silac_clusters_long$gene_label != "HSP90")
is_APMS_client <- silac_clusters_long$feature_name %in% HSP90_APMS_clients

z <- table(in_HSP90_cluster, is_APMS_client)

z %>% fisher.test
# p-value = 0.01502
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  1.101122 5.735977
# sample estimates:
# odds ratio
#   2.651136


overlap <- silac_clusters_long %>%
	dplyr::filter(
		in_HSP90_cluster,
		is_APMS_client) %>%
	magrittr::extract2("gene_label") %>%
	paste0(collapse=" ")
# ACC1 BMH1 CLU1 ECM17 FAS1 FAS2 SRV2 SUB2 TIF
