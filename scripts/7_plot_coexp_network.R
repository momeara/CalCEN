# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(stringr)
library(readr)
library(seriation)
library(EGAD)

source("scripts/seriation_heatmap.R")

load("intermediate_data/ca_runs_final.Rdata")

load("intermediate_data/estimated_expression.Rdata")
load("product/CalCEN_network_full_20201007.Rdata")

date_code <- "20201130"

#############################
# Plot co-expression matrix #
#############################

load("intermediate_data/CalCEN.Rdata")

### OLO ##
gene_order <- seriation::seriate(
		x = dist(1-CalCEN_network_full),
		method = "OLO") %>%
		seriation::get_order()

heatmap <- as.raster(CalCEN_network_full[gene_order, gene_order])
png(
		"product/figures/CalCEN_heatmap_raster1_20201130.png",
		width = length(gene_order),
		height = length(gene_order))
plot.new()
rasterImage(
		image = heatmap,
		xleft = 0,
		ybottom = 0,
		xright = 1,
		ytop = 1,
		interpolate = FALSE)
dev.off()

### TSP ##
gene_order <- seriation::seriate(
		x = dist(1-CalCEN_network_full),
		method = "TSP") %>%
		seriation::get_order()

heatmap <- as.raster(CalCEN_network_full[gene_order, gene_order])
png(
		"product/figures/CalCEN_heatmap_raster_TSP_20201112.png",
		width = length(gene_order),
		height = length(gene_order))
plot.new()
rasterImage(
		image = heatmap,
		xleft = 0,
		ybottom = 0,
		xright = 1,
		ytop = 1,
		interpolate = FALSE)
dev.off()

###########################
# Plot expression heatmap #
###########################


estimated_expression <- estimated_expression %>%
		dplyr::semi_join(ca_runs_final, by = c("run_accession"))

exprs <- reshape2::acast(
	data=estimated_expression,
	formula=gene_id ~ run_accession,
	value.var="FPKM")

exprs <- log(1+exprs)


run_network <- EGAD::build_coexp_network(
		exprs = exprs %>% t(),
		gene.list = colnames(exprs))

run_order <- seriation::seriate(
		x = dist(1-run_network),
		method = "OLO") %>%
		seriation::get_order()

exprs_log1p <- log(1+exprs)
heatmap <- as.raster(exprs_log1p[gene_order, run_order]/max(exprs_log1p))

png(
		"product/figures/ca_expression_heatmap_raster_20201130.png",
		width = length(run_order),
		height = length(gene_order))
plot.new()
rasterImage(
		image = heatmap,
		xleft = 0,
		ybottom = 0,
		xright = 1,
		ytop = 1,
		interpolate = FALSE)
dev.off()

