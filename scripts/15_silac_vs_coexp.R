# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(readr)
library(reshape2)
library(EGAD)
library(ggplot2)


load("intermediate_data/ca_silac_hsp90_intensities.Rdata")
load("intermediate_data/ca_genes.Rdata")
load("intermediate_data/ca_silac_network.Rdata")
load("intermediate_data/ca_coexp_full.Rdata")


spearman_rank_correlation <- function(netA, netB, n){
	netA_ranks <- order(netA)
	induced_netB_ranks <- order( netB[netA_ranks] )
	z <- cor(induced_netB_ranks, 1:n, method="spearman")
}


silac_genes <- ca_silac_network %>% rownames()
coexp_genes <- ca_coexp_full %>%
	dplyr::distinct(feature_name_1) %>%
	magrittr::extract2("feature_name_1")
genes <- intersect(silac_genes, coexp_genes)

ca_silac_network <- ca_silac_network %>%
	embed_network(genes)

ca_coexp_network <- ca_coexp_full %>%
	dplyr::filter(feature_name_1 %in% genes) %>%
	dplyr::filter(feature_name_2 %in% genes) %>%
	reshape2::acast(feature_name_1 ~ feature_name_2, value.var="score") %>%
	embed_network(genes)

n <- length(genes) * length(genes)
spearman_rank_correlation(ca_coexp_network, ca_silac_network, n)




		

rankM <- function(network,n){
	network <- matrix(rank(network, na.last="keep",ties.method="average"), nrow=n, ncol=n)
	network <- Matrix(network, sparse=T)
}

### Ranks and unwinds sparse matrices
rankAndFlattenM <- function(network,n, fill=0){
	#diag(network) = 0
	diag(network) = fill
	temp = rankM(network,n)
	diag(temp) = NA
	flat_net = array(temp)
}






do_spearman_test <- function(
	netA, netB, n, B=999, plot=T, fname="spearman_rank_perm_test.pdf"){

	test <- function(y){
		suppressWarnings(
			cor.test(
				rankAndFlattenM(netB,n),
				y,
				method="spearman")$estimate)
	}

	fr_netB <- rankAndFlattenM(netB,n)

	rho <- test(fr_netB)

	nulldist <- data.frame(p=replicate(B, test( sample( fr_netB, length(fr_netB)))))
	if(plot) {
		p <- ggplot(data=nulldist) +
			theme_bw() +
			ggtitle("Spearman Rank Correlation Permutation Test") +
			geom_histogram(aes(x=p)) +
			geom_vline(xintercept=rho, color="blue", size=1.4) +
			scale_x_continuous("Spearman Rank Correlation Test Statistic") +
			scale_y_continuous("Density")
		ggsave(fname)
	}
}
