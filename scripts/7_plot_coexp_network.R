# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(stringr)
library(readr)
library(seriation)

load("intermediate_data/estimated_expression.Rdata")


exprs <- reshape2::acast(
	data=estimated_expression,
	formula=gene_id ~ run_accession,
	value.var="FPKM")

pdf("product/figures/expression_heatmap_180612.png")
heat_map <- seriation::hmap((exprs+1) %>% log())
dev.off()
