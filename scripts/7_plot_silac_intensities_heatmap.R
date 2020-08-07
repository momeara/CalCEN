# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(tidyr)
library(tibble)
library(seriation)
library(gplots)

load("intermediate_data/chromosome_features.Rdata")
load("intermediate_data/ca_silac_hsp90_intensities.Rdata")
load("intermediate_data/ca_silac_hsp90_fold_change.Rdata")

load("intermediate_data/genes_of_interest.Rdata")

#########
# given a matrix x,
# This makes a dendrogram where the columns are sorted to minimize the row-to-row difference
# the sorting is done relative to the ref_x, if it is supplied
silac_heatmap <- function (
	x,
	ref_x=x,
	fname=NULL,
	height=180,
	width=15)
{
	if(!is.null(ref_x)){
		dist_row <- dist(ref_x)
		o_row <- seriate(dist_row, method = "OLO", control = NULL)[[1]]
		Rowv=as.dendrogram(o_row)
	} else {
		Rowv = FALSE
	}
	args <- list()
  if (any(x < 0, na.rm = TRUE)){
		args$col <- seriation::bluered(n=100, bias=1)
	} else {
		args$col <- seriation::greys(n=100, power=1)
	}
	args$trace <- "none"
	args$density.info <- "none"
	args$cexRow <- 1
	args$cexCol <- 1
	args$dendrogram="none"
	args$key = FALSE
	args$keysize = 0.03
	args$colsep = seq(0, ncol(x), by=5000000)
	args$rowsep = seq(0, nrow(x), by=5000000)
	#args$sepwidth = c(0.02, 0.02)
	args$sepwidth = c(0, 0)
	args$margins = c(3,7)
	args <- c(list(x = x, Colv = FALSE, Rowv = Rowv), args)
	cat("plotting '", fname, "'...\n", sep="")
	pdf(fname, height=height, width=width)
		suppressWarnings(ret <- do.call(gplots::heatmap.2, args))
	dev.off()
}


### WT Intensities ###
# log intensity
wt_data <- ca_silac_hsp90_intensities %>%
	dplyr::filter(condition=="wildtype") %>%
	dplyr::mutate(log_intensity=log(intensity + 1)) %>%
	dplyr::select(gene_label, fraction_index, log_intensity) %>%
	tidyr::spread(fraction_index, log_intensity) %>%
	tibble::column_to_rownames("gene_label") %>%
	as.matrix()

silac_heatmap(
	x=wt_data,
	fname="product/figures/silac_hsp90_intensities_wildtype_heatmap_190219.pdf")

# normed intensity
wt_data <- ca_silac_hsp90_intensities %>%
	dplyr::filter(condition=="wildtype") %>%
	dplyr::group_by(gene_label) %>%
	dplyr::mutate(normed_intensity= intensity / max(intensity)) %>%
	dplyr::ungroup() %>%
	dplyr::mutate(
		normed_intensity = ifelse(is.na(normed_intensity), 0, normed_intensity)) %>%
	dplyr::filter(!is.na(normed_intensity)) %>%
	dplyr::select(gene_label, fraction_index, normed_intensity) %>%
	tidyr::spread(fraction_index, normed_intensity) %>%
	tibble::column_to_rownames("gene_label") %>%
	as.matrix()

silac_heatmap(
	x=wt_data,
	fname="product/figures/silac_hsp90_normed_intensities_wildtype_heatmap_190219.pdf")



### HSP90 depletion fold change ###
# log_fc
ph_fc_data <- ca_silac_hsp90_fold_change %>%
	dplyr::select(gene_label, fraction_index, pharmacological_log_fc) %>%
	tidyr::spread(fraction_index, pharmacological_log_fc) %>%
	tibble::column_to_rownames("gene_label") %>%
	as.matrix()

silac_heatmap(
	x=ph_fc_data,
	ref_x=wt_data,
	fname="product/figures/silac_hsp90_fold_change_pharmacolical_heatmap_190219.pdf")

gt_fc_data <- ca_silac_hsp90_fold_change %>%
	dplyr::select(gene_label, fraction_index, genetic_log_fc) %>%
	tidyr::spread(fraction_index, genetic_log_fc) %>%
	tibble::column_to_rownames("gene_label") %>%
	as.matrix()

silac_heatmap(
	x=gt_fc_data,
	ref_x=wt_data,
	fname="product/figures/silac_hsp90_fold_change_genetic_heatmap_190219.pdf")

# normed log_fc
ph_fc_data <- ca_silac_hsp90_fold_change %>%
	dplyr::select(gene_label, fraction_index, pharmacological_log_fc) %>%
	dplyr::group_by(gene_label) %>%
	dplyr::mutate(
    max_log_fc = max(abs(pharmacological_log_fc)),
		normed_pharmacological_log_fc =	pharmacological_log_fc / max_log_fc) %>%
	dplyr::ungroup() %>%
	dplyr::mutate(normed_pharmacological_log_fc = ifelse(
		is.na(normed_pharmacological_log_fc), 0,
		normed_pharmacological_log_fc)) %>%
	dplyr::select(-pharmacological_log_fc, -max_log_fc) %>%
	tidyr::spread(fraction_index, normed_pharmacological_log_fc) %>%
	tibble::column_to_rownames("gene_label") %>%
	as.matrix()

gt_fc_data <- ca_silac_hsp90_fold_change %>%
	dplyr::select(gene_label, fraction_index, genetic_log_fc) %>%
	dplyr::group_by(gene_label) %>%
	dplyr::mutate(
    max_log_fc = max(abs(genetic_log_fc)),
		normed_genetic_log_fc =	genetic_log_fc / max_log_fc) %>%
	dplyr::ungroup() %>%
	dplyr::mutate(normed_genetic_log_fc = ifelse(
		is.na(normed_genetic_log_fc), 0,
		normed_genetic_log_fc)) %>%
	dplyr::select(-genetic_log_fc, -max_log_fc) %>%
	tidyr::spread(fraction_index, normed_genetic_log_fc) %>%
	tibble::column_to_rownames("gene_label") %>%
	as.matrix()



silac_heatmap(
	x=wt_data,
	width=15,
	height=180,
	fname="product/figures/silac_hsp90_normed_intensities_wildtype_heatmap_190219.pdf")

silac_heatmap(
	x=ph_fc_data,
	ref_x=wt_data,
	width=15,
	height=180,
	fname="product/figures/silac_hsp90_normed_fold_change_pharmacolical_heatmap_190219.pdf")


silac_heatmap(
	x=gt_fc_data,
	ref_x=wt_data,
	width=15,
	height=180,
	fname="product/figures/silac_hsp90_normed_fold_change_genetic_heatmap_190219.pdf")



### Combined WT/HSP90 depletion fold changes for figure 1 ###

silac_heatmap(
	x=wt_data,
	width=15,
	height=40,
	fname="product/figures/silac_hsp90_fig1_wildtype_intensities_heatmap_180808.pdf")


data <- cbind(ph_fc_data, gt_fc_data)
data <- log(1 + abs(data))*sign(data)
silac_heatmap(
	x=data,
	ref_x=wt_data,
	width=30,
	height=40,
	fname="product/figures/silac_hsp90_fig1_genetic_pharmacolical_fold_change_heatmap_180808.pdf")

silac_heatmap(
		x=matrix(c(
			-1, -.75, -.5, -.25, 0, .25, .5, .75, 1,
			-1, -.75, -.5, -.25, 0, .25, .5, .75, 1,
			-1, -.75, -.5, -.25, 0, .25, .5, .75, 1,
			-1, -.75, -.5, -.25, 0, .25, .5, .75, 1,
			-1, -.75, -.5, -.25, 0, .25, .5, .75, 1,
			-1, -.75, -.5, -.25, 0, .25, .5, .75, 1,
			-1, -.75, -.5, -.25, 0, .25, .5, .75, 1,
			-1, -.75, -.5, -.25, 0, .25, .5, .75, 1,
			-1, -.75, -.5, -.25, 0, .25, .5, .75, 1,
			-1, -.75, -.5, -.25, 0, .25, .5, .75, 1,
			-1, -.75, -.5, -.25, 0, .25, .5, .75, 1,
			-1, -.75, -.5, -.25, 0, .25, .5, .75, 1,
			-1, -.75, -.5, -.25, 0, .25, .5, .75, 1,
			-1, -.75, -.5, -.25, 0, .25, .5, .75, 1,
			-1, -.75, -.5, -.25, 0, .25, .5, .75, 1), nrow=15, byrow=TRUE),
		width=15,
		height=2,
		fname="product/figures/silac_normed_fc_legend.pdf")


silac_heatmap(
		x=matrix(c(
			0, 0.125, .25, 0.375, .5, 0.625, .75, 0.875, 1,
			0, 0.125, .25, 0.375, .5, 0.625, .75, 0.875, 1,
			0, 0.125, .25, 0.375, .5, 0.625, .75, 0.875, 1,
			0, 0.125, .25, 0.375, .5, 0.625, .75, 0.875, 1,
			0, 0.125, .25, 0.375, .5, 0.625, .75, 0.875, 1,
			0, 0.125, .25, 0.375, .5, 0.625, .75, 0.875, 1,
			0, 0.125, .25, 0.375, .5, 0.625, .75, 0.875, 1,
			0, 0.125, .25, 0.375, .5, 0.625, .75, 0.875, 1,
			0, 0.125, .25, 0.375, .5, 0.625, .75, 0.875, 1,
			0, 0.125, .25, 0.375, .5, 0.625, .75, 0.875, 1,
			0, 0.125, .25, 0.375, .5, 0.625, .75, 0.875, 1,
			0, 0.125, .25, 0.375, .5, 0.625, .75, 0.875, 1,
			0, 0.125, .25, 0.375, .5, 0.625, .75, 0.875, 1,
			0, 0.125, .25, 0.375, .5, 0.625, .75, 0.875, 1,
			0, 0.125, .25, 0.375, .5, 0.625, .75, 0.875, 1), nrow=15, byrow=TRUE),
		width=15,
		height=2,
		fname="product/figures/silac_normed_wt_legend.pdf")



#############

genes_of_interest %>%
	dplyr::filter(set %in% c("R2TP Complex")) %>%
	dplyr::filter(include) %>%
	plyr::d_ply(c("set"), function(genes){
	set_id <- genes$set[1]
	cat("plotting change parmacological log fold change heatmaps for ", set_id, " ...\n", sep="")
	ph_fc_data <- genes %>%
		dplyr::select(feature_name, gene_name) %>%
		dplyr::left_join(ca_silac_hsp90_fold_change, by=c("feature_name")) %>%
		dplyr::mutate(gene_label = ifelse(!is.na(gene_label), gene_label, gene_name)) %>%
		# normalize intensity per-gene
		dplyr::group_by(gene_label) %>%
		dplyr::mutate(
  	  max_log_fc = max(abs(pharmacological_log_fc)),
			normed_pharmacological_log_fc =	pharmacological_log_fc / max_log_fc) %>%
		dplyr::ungroup() %>%
		dplyr::mutate(
			normed_pharmacological_log_fc = ifelse(
				is.na(normed_pharmacological_log_fc), 0,
				normed_pharmacological_log_fc)) %>%
		dplyr::select(-pharmacological_log_fc, -max_log_fc) %>%
		#
		dplyr::select(gene_label, fraction_index, normed_pharmacological_log_fc) %>%
		dplyr::filter(!is.na(fraction_index)) %>%
		tidyr::spread(fraction_index, normed_pharmacological_log_fc) %>%
		tibble::column_to_rownames("gene_label") %>%
		as.matrix()
	print(ph_fc_data)
	silac_heatmap(
		x=ph_fc_data,
		ref_x=NULL,
		fname=paste0("product/figures/silac_hsp90_fold_change_pharmacolical_heatmap_", set_id, "_190304.pdf"),
		width=15,
		height=nrow(ph_fc_data))
	cat("plotting genetic log fold change heatmaps for ", set_id, " ...\n", sep="")
	gt_fc_data <- genes %>%
		dplyr::select(feature_name, gene_name) %>%
		dplyr::left_join(ca_silac_hsp90_fold_change, by=c("feature_name")) %>%
		dplyr::mutate(gene_label = ifelse(!is.na(gene_label), gene_label, gene_name)) %>%
		# normalize intensity per-gene
		dplyr::group_by(gene_label) %>%
		dplyr::mutate(
  	  max_log_fc = max(abs(genetic_log_fc)),
			normed_genetic_log_fc =	genetic_log_fc / max_log_fc) %>%
		dplyr::ungroup() %>%
		dplyr::mutate(normed_genetic_log_fc = ifelse(
			is.na(normed_genetic_log_fc), 0,
			normed_genetic_log_fc)) %>%
		dplyr::select(-genetic_log_fc, -max_log_fc) %>%
		#
		dplyr::select(gene_label, fraction_index, normed_genetic_log_fc) %>%
		dplyr::filter(!is.na(fraction_index)) %>%
		tidyr::spread(fraction_index, normed_genetic_log_fc) %>%
		tibble::column_to_rownames("gene_label") %>%
		as.matrix()
	silac_heatmap(
		x=gt_fc_data,
		ref_x=NULL,
		fname=paste0("product/figures/silac_hsp90_fold_change_genetic_heatmap_", set_id, "_190304.pdf"),
		width=15,
		height=nrow(gt_fc_data))
})
