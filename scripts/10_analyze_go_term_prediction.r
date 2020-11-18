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


go_pred_auc_tbl <- go_pred[[1]][[1]] %>%
		as.data.frame() %>%
		tibble::rownames_to_column(var = "go_id") %>%
		dplyr::left_join(
				go_terms,
				by = c("go_id"))

go_pred_auc_tbl %>%
	readr::write_tsv("product/go_pred_coexp_20201116.tsv")


p <- ggplot2::ggplot(data = go_pred_auc_tbl) +
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

ggsave(
	file="product/figures/go_pred_coex_scatter_p_20201117.pdf",
	height=4, width=4)
ggsave(
	file="product/figures/go_pred_coexp_scatter_20201117.png",
	height=4, width=4)



# Plot the ROC AUC for each term
# For each (gene, term) pair, the neighbor voting predictor sums the edge weight for all neighbors that have the term divided by weight of all edges.
# For each (gene, term) the degree null predictor assigns the weight of all edges to in the network to the gene, irrespective of the term.
p <- ggplot(data=go_pred_auc_tbl %>%
	dplyr::select(
		go_id,
		ontology,
		`Neighbor Voting`=auc,
		`Degree Null`=degree_null_auc) %>%
	tidyr::gather("prediction_method", "auc", -go_id, -ontology)) +
	theme_bw() +
	geom_density(
		mapping=aes(x=auc, fill=prediction_method),
		alpha=.8) +
	ggtitle("Predictability of GO terms by Co-Expression Network") +
	scale_x_continuous("Area Under the ROC curve") +
	scale_fill_discrete("Prediction Method")

ggsave(
	file="product/figures/go_pred_coexp_20201117.pdf",
	height=3, width=5)
ggsave(
	file="product/figures/go_pred_coexp_20201117.png",
	height=3, width=5)


p <- p + facet_wrap(~ontology, ncol=1)

ggsave(
	file="product/figures/go_pred_coexp_by_ontology_20201117.pdf",
	height=5, width=5)
ggsave(
	file="product/figures/go_pred_coexp_by_ontology_20201117.png",
	height=5, width=5)






