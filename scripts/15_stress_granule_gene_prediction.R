# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(stringr)
library(readr)
library(magrittr)
library(ggplot2)
library(ggrepel)
library(UpSetR)
load("intermediate_data/genes_of_interest.Rdata")


make_plot <- function(genes_of_interest){
	genes <- genes_of_interest %>%
		dplyr::distinct(feature_name, set, .keep_all=TRUE)
	list(
		`Candida Stress Granule` =
			genes %>% dplyr::filter(set=="Stress Granule") %>% magrittr::extract2("gene_name"),
		`Sac GO:Stress Granule` =
			genes %>% dplyr::filter(set=="Sac Ortholog Stress Granule") %>% magrittr::extract2("gene_name"),
		`Sac GO:P-Body` =
			genes %>% dplyr::filter(set=="Sac Ortholog P-body") %>% magrittr::extract2("gene_name"),
		`Human Stress Granule` =
			genes %>% dplyr::filter(set=="Predicted Stress Granule") %>% magrittr::extract2("gene_name")) %>%
		UpSetR::fromList() %>%
		UpSetR::upset(
			keep.order=TRUE,
			point.size=3.5,
			mainbar.y.label="Candida Albicans Genes",
			sets.x.label="Evidence Source")
}

pdf(
   file="product/figures/predicted_stress_granules_180822.pdf",
   width=4, height=4, useDingbats=FALSE, onefile=FALSE)
make_plot(genes_of_interest)
dev.off()

png(
   file="product/figures/predicted_stress_granules_180822.png",
   width=360*4, height=360*4, res=72*4)
make_plot(genes_of_interest)
dev.off()


stress_granule_table <- genes_of_interest %>%
	dplyr::distinct(feature_name, set, .keep_all=TRUE) %>%
	dplyr::select(-notes, -reference) %>%
	dplyr::mutate(value=TRUE) %>%
	tidyr::spread(
		value="value",
		key="set",
		fill=FALSE)

stress_granule_table %>%
	readr::write_tsv("product/predicted_stress_granules_180822.tsv")
