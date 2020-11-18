# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(magrittr)
library(stringr)
library(readr)
library(reshape2)
library(EGAD)
library(Hotelling)
library(huge)

load("intermediate_data/ca_runs_final.Rdata")
load("intermediate_data/chromosome_features.Rdata")
load("intermediate_data/estimated_expression.Rdata")
load("intermediate_data/ca_genes.Rdata")

exprs <- reshape2::acast(
		data=estimated_expression %>%
				dplyr::semi_join(ca_runs_final, by=c("run_accession")),
		formula=gene_id ~ run_accession,
		value.var="FPKM") %>%
		magrittr::extract(ca_genes,)

#############################
# spearman rank correlation #
#############################

ca_coexp_network_full_spearman <- EGAD::build_coexp_network(
	exprs=exprs,
	gene.list=ca_genes)

save(ca_coexp_network_full, file="product/ca_coexp_network_full_spearman_20201024.Rdata")
ca_coexp_network_full_spearman %>%
	as.data.frame() %>%
	readr::write_tsv("product/ca_coexp_network_full_spearman_20201024.tsv")

ca_coexp_full_spearman <- ca_coexp_network_full_spearman %>%
	data.frame() %>%
	tibble::rownames_to_column("feature_name_1") %>%
	tidyr::gather(key="feature_name_2", value="score", -feature_name_1)
save(ca_coexp_full, file="intermediate_data/ca_coexp_full_spearman.Rdata")



###################################
# Pearson Correlation coefficient #
###################################

ca_coexp_network_full_pearson <- EGAD::build_coexp_network(
    exprs = exprs,
    gene.list=ca_genes,
    method = "pearson",
    flag = FALSE)

save(ca_coexp_network_full_pearson, file="product/ca_coexp_network_full_pearson_20201024.Rdata")
ca_coexp_network_full_pearson %>%
	as.data.frame() %>%
	readr::write_tsv("product/ca_coexp_network_full_pearson_20201024.tsv")

ca_coexp_full_pearson <- ca_coexp_network_full_pearson %>%
    data.frame() %>%
    tibble::rownames_to_column("feature_name_1") %>%
    tidyr::gather(key="feature_name_2", value="score", -feature_name_1)
save(ca_coexp_full, file="intermediate_data/ca_coexp_full_pearson.Rdata")

######################
# Direct correlation #
######################

# more numerically stable version of Hotelling::clr
# https://stackoverflow.com/q/2602583/198401
clr <- function(data){
		log_gms <- apply(data, 1, function(x){mean(log(x))})
		log(data) - log_gms
}

nlambda <- 30
ca_coexp_full_direct <- (exprs+1) %>%
		t() %>%
		clr() %>%
		huge::huge.npn() %>%
		huge::huge(method="mb", nlambda = nlambda)
for (i in c(1:nlambda)) {
		rownames(ca_coexp_full_direct$path[[i]]) <- ca_genes
		colnames(ca_coexp_full_direct$path[[i]]) <- ca_genes
		rownames(ca_coexp_full_direct$beta[[i]]) <- ca_genes
		colnames(ca_coexp_full_direct$beta[[i]]) <- ca_genes
}

save(ca_coexp_full_direct, file="intermediate_data/ca_coexp_full_direct.Rdata")

graph_bootstrap <- graph %>%
		huge::huge.select(criterion = "stars", stars.thresh=0.05)
save(graph_bootstrap, file="intermediate_data/ca_coexp_full_direct_bootstrap.Rdata")
