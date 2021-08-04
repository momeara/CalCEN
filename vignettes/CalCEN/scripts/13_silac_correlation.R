# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(tidyr)
library(tibble)
library(EGAD)
library(seriation)
library(gplots)

load("intermediate_data/ca_silac_hsp90_intensities.Rdata")

silac_genes <- ca_silac_hsp90_intensities %>%
	dplyr::distinct(gene_label) %>%
	magrittr::extract2("gene_label")


input_condition <- "wildtype"
fname <- "product/figures/silac_hsp90_correlation_heatmap_190121.pdf"
height=30
width=30
silac_network <- ca_silac_hsp90_intensities %>%
	dplyr::filter(condition==input_condition) %>%
	dplyr::select(gene_label, fraction_index, intensity) %>%
	tidyr::spread(fraction_index, intensity) %>%
	tibble::column_to_rownames("gene_label") %>%
	as.matrix() %>%
	EGAD::build_coexp_network(silac_genes)

silac_network[is.na(silac_network)] <- 0

o_network <- seriate(as.dist(silac_network), method = "OLO", control = NULL)[[1]]
args <- list()
args$trace <- "none"
args$density.info <- "none"
args$cexRow <- 1
args$cexCol <- 1
args$dendrogram="none"
args$key = FALSE
args$keysize = 0.03
args$colsep = seq(0, ncol(silac_network), by=5000000)
args$rowsep = seq(0, nrow(silac_network), by=5000000)
#args$sepwidth = c(0.02, 0.02)
args$sepwidth = c(0, 0)
args$margins = c(3,7)
#args <- c(list(x = silac_network, Colv = FALSE, Rowv = o_network), args)
args <- c(list(x = silac_network, Colv = FALSE, Rowv = as.dendrogram(o_network)), args)


cat("plotting '", fname, "'...\n", sep="")
pdf(fname, height=height, width=width)
	suppressWarnings(ret <- do.call(gplots::heatmap.2, args))
dev.off()
