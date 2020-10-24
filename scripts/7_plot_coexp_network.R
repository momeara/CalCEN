# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(stringr)
library(readr)
library(seriation)

source("scripts/seriation_heatmap.R")

load("intermediate_data/ca_runs_final.Rdata")

load("intermediate_data/estimated_expression.Rdata")
load("product/ca_coexp_network_full_20201007.Rdata")
date_code <- "20201022"

###########################
# Plot expression heatmap #
###########################
estimated_expression <- estimated_expression %>%
		dplyr::semi_join(ca_runs_final, by = c("run_accession"))

exprs <- reshape2::acast(
	data=estimated_expression,
	formula=gene_id ~ run_accession,
	value.var="FPKM")

exprs <- log10(exprs+1)

gene_slice <- 1:2000
run_slice <- 1:50
ca_coexp_dist <- as.dist(1-ca_coexp_network_full[,])
seriation_heatmap(
		x = exprs[,],
		#ref_x = ca_coexp_dist,
		color_scale = seriation::greys(100, power = 2),
		fname = paste0("product/figures/test_heatmap.pdf"),
		height = 8,
		width = 8,
		verbose = TRUE)

ca_coexp_dist <- as.dist(1-ca_coexp_network_full)
seriation_heatmap(
		x = exprs,
		ref_x = ca_coexp_dist,
		color_scale = seriation::greys(100, power = 2),
		fname = paste0("product/figures/expression_heatmap_", date_code, ".pdf"),
		height = 150,
		width = 30,
		verbose = TRUE)

#############################
# Plot co-expression matrix #
#############################

pdf(
heat_map <- seriation::hmap(ca_coexp_network_full)
dev.off()
