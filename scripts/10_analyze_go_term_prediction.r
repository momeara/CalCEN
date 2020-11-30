# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(GO.db)

load("intermediate_data/go_pred.Rdata")


# these are the current go-annotations but
# they may not be in sync with the terms provided by CGD
go_terms_current <- GO.db::GO_dbconn() %>%
	dplyr::tbl("go_term") %>%
	dplyr::collect(n=Inf) %>%
	dplyr::select(-`_id`)

load("intermediate_data/ca_go_propagated_filtered.Rdata")
go_terms <- ca_go_propagated_filtered %>%
		dplyr::distinct(go_id, .keep_all=TRUE) %>%
		dplyr::select(go_id, term, ontology, definition)

network_scores <- go_pred %>%
		purrr::map_dfr(
		.id = "network_name",
		.f = function(network_pred) {
			go_pred_auc_tbl <- network_pred[[1]] %>%
					as.data.frame() %>%
					tibble::rownames_to_column(var = "go_id")}) %>%
		dplyr::left_join(go_terms, by = c("go_id")) %>%
		dplyr::mutate(network_name = network_name %>% stringr::str_replace("_go_pred", ""))


network_scores %>%
	readr::write_tsv("product/go_pred_20201118.tsv")


p <- ggplot2::ggplot(
		data = network_scores %>%
				dplyr::filter(network_name == "ca_coexp_full")) +
		ggplot2::theme_bw() +
		ggplot2::theme(legend.position = c(.8, .2)) +
		ggplot2::geom_abline(
				slope = 1,
				color = "grey70") +
		ggplot2::geom_point(
				mapping = ggplot2::aes(
						x = degree_null_auc,
						y = auc,
						color = ontology),
				size = 1) +
		ggplot2::coord_fixed() +
		ggplot2::scale_x_continuous(
				"Degree-Null AUROC",
				limits = c(0, 1)) +
		ggplot2::scale_y_continuous(
				"Co-Expression Network AUROC",
				limits = c(0, 1)) +
		ggplot2::ggtitle(
				label = "GO term prediction : Co-Expession")

ggplot2::ggsave(
	file="product/figures/go_pred_coexp_scatter_20201117.pdf",
	height=4, width=4)
ggplot2::ggsave(
	file="product/figures/go_pred_coexp_scatter_20201117.png",
	height=4, width=4)



p <- ggplot2::ggplot(
		data = network_scores %>%
				dplyr::filter(
						network_name %in% c(
								"ca_coexp_full",
								"ca_blastp",
								"ca_sac_ortholog_genetic_ppi",
								"ca_sac_ortholog_physical_ppi",
								"ca_yeast_net"))) +
		ggplot2::theme_bw() +
		ggplot2::theme(legend.position = "bottom") +
		ggplot2::geom_abline(
				slope = 1,
				color = "grey70") +
		ggplot2::geom_point(
				mapping = ggplot2::aes(
						x = degree_null_auc,
						y = auc,
						color = ontology),
				size = 1) +
		ggplot2::coord_fixed() +
		ggplot2::scale_x_continuous(
				"Degree-Null AUROC",
				limits = c(0, 1)) +
		ggplot2::facet_wrap(facets = dplyr::vars(network_name)) +
		ggplot2::scale_y_continuous(
				"Network AUROC",
				limits = c(0, 1)) +
		ggplot2::ggtitle(label = "GO term prediction")

ggplot2::ggsave(
	file="product/figures/go_pred_scatter_20201117.pdf",
	height=10, width=10)
ggplot2::ggsave(
	file="product/figures/go_pred_scatter_20201117.png",
	height=10, width=10)






