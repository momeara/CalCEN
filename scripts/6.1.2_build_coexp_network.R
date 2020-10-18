# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(magrittr)
library(stringr)
library(readr)
library(reshape2)
library(EGAD)

load("intermediate_data/chromosome_features.Rdata")
load("intermediate_data/estimated_expression.Rdata")
load("intermediate_data/ca_genes.Rdata")

exprs <- reshape2::acast(
	data=estimated_expression,
	formula=gene_id ~ run_accession,
	value.var="FPKM") %>%
	magrittr::extract(ca_genes,)

ca_coexp_network_full <- EGAD::build_coexp_network(
	exprs=exprs,
	gene.list=ca_genes)

save(ca_coexp_network_full, file="product/ca_coexp_network_full_20201007.Rdata")
ca_coexp_network_full %>%
	as.data.frame() %>%
	readr::write_tsv("product/ca_coexp_network_full_20201007.tsv")

ca_coexp_full <- ca_coexp_network_full %>%
	data.frame() %>%
	tibble::rownames_to_column("feature_name_1") %>%
	tidyr::gather(key="feature_name_2", value="score", -feature_name_1)
save(ca_coexp_full, file="intermediate_data/ca_coexp_full.Rdata")
