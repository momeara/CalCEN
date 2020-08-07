# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(GO.db)


load("intermediate_data/go_pred.Rdata")

go_terms <- GO.db::GO_dbconn() %>%
	dplyr::tbl("go_term") %>%
	dplyr::collect(n=Inf)


go_pred_auc_tbl <- go_pred[[1]] %>%
	as.data.frame() %>%

go_pred_auc_tbl %>%
	readr::write_tsv("product/go_pred_coexp_180618.tsv")


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
	file="product/figures/go_pred_coexp_180618.pdf",
	height=5, width=5)
ggsave(
	file="product/figures/go_pred_coexp_180618.png",
	height=5, width=5)


p <- p + facet_wrap(~ontology, ncol=1)

ggsave(
	file="product/figures/go_pred_coexp_by_ontology_180618.pdf",
	height=5, width=5)
ggsave(
	file="product/figures/go_pred_coexp_by_ontology_180618.png",
	height=5, width=5)






