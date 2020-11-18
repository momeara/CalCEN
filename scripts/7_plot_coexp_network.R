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
date_code <- "20201110"

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

gene_slice <- 1:6226
run_slice <- 1:853
#ca_coexp_dist <- as.dist(exprs[gene_slice,run_slice])
for(seriate_method in c("TSP", "OLO", "VAT")){
  seriation_heatmap(
  		x = exprs[gene_slice,run_slice],
  		ref_x = exprs[gene_slice,run_slice],
  		color_scale = seriation::greys(100, power = 1),
  		fname = paste0(
					"product/figures/test_heatmap_seriation_g=", length(gene_slice), ",r=",length(run_slice), "_method=", seriate_method,"_h=40.pdf"),
  		height = 40,
  		width = 6,
			seriate_method = seriate_method,
  		verbose = TRUE)
	}

ca_coexp_dist <- as.dist(1-ca_coexp_network_full[1:200, 1:200])

seriation_heatmap(
		x = 1-ca_coexp_network_full[1:200, 1:200],
		ref_x = 1-ca_coexp_network_full[1:200, 1:200],
		color_scale = seriation::greys(100, power = 2),
		fname = paste0("product/figures/ca_coexp_heatmap_", date_code, ".pdf"),
		height = 10,
		width = 10,
		verbose = TRUE)

seriation_heatmap(
		x = 1-ca_coexp_network_full[1:200, 1:200],
		ref_x = 1-ca_coexp_network_full[1:200, 1:200],
		color_scale = seriation::greys(100, power = 2),
		fname = paste0("product/figures/ca_coexp_heatmap_raster_", date_code, ".png"),
		height = 10,
		width = 10,
		verbose = TRUE)



#############################
# Plot co-expression matrix #
#############################

load("intermediate_data/ca_coexp.Rdata")

### OLO ##
gene_order <- seriation::seriate(
		x = dist(1-ca_coexp_network_full),
		method = "OLO") %>%
		seriation::get_order()
heatmap <- as.raster(ca_coexp_network[gene_order, gene_order])
png(
		"product/figures/ca_coexp_heatmap_raster1_20201112.png",
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
		x = dist(1-ca_coexp_network_full),
		method = "TSP") %>%
		seriation::get_order()

heatmap <- as.raster(ca_coexp_network_full[gene_order, gene_order])
png(
		"product/figures/ca_coexp_heatmap_raster_TSP_20201112.png",
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

