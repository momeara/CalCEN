g# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(UpSetR)

load("intermediate_data/chromosome_features.Rdata")
load("intermediate_data/estimated_expression.Rdata")
load("intermediate_data/ca_coexp.Rdata")
load("intermediate_data/ca_blastp.Rdata")
load("intermediate_data/ca_biogrid.Rdata")
load("intermediate_data/sac_biogrid.Rdata")


load("intermediate_data/yeast_net.Rdata")
load("intermediate_data/yeast_net_CC.Rdata")
load("intermediate_data/yeast_net_CX.Rdata")
load("intermediate_data/yeast_net_DC.Rdata")
load("intermediate_data/yeast_net_GN.Rdata")
load("intermediate_data/yeast_net_GT.Rdata")
load("intermediate_data/yeast_net_HT.Rdata")
load("intermediate_data/yeast_net_LC.Rdata")
load("intermediate_data/yeast_net_PG.Rdata")
load("intermediate_data/yeast_net_TS.Rdata")

load("intermediate_data/ca_go_propagated_filtered.Rdata")




ca_coexp_1p <- ca_coexp %>%
		dplyr::filter(
				feature_name_1 != feature_name_2,
				score >= .99)

ca_blastp_5E <- ca_blastp %>%
		dplyr::filter(
				feature_name_1 != feature_name_2,
				EValue <= 1e-5)

sac_ortholog_ppi <- sac_biogrid %>%
	dplyr::select(
		sac_ortholog_1=gene_symbol_1,
		sac_ortholog_2=gene_symbol_2,
		experimental_system_type) %>%
	dplyr::inner_join(
		chromosome_features %>%
			dplyr::filter(feature_class=="ORF") %>%
			dplyr::select(feature_name_1=feature_name, sac_ortholog_1 = sac_ortholog),
		by="sac_ortholog_1") %>%
	dplyr::inner_join(
		chromosome_features %>%
			dplyr::filter(feature_class=="ORF") %>%
			dplyr::select(feature_name_2=feature_name, sac_ortholog_2 = sac_ortholog),
		by="sac_ortholog_2")

network_overlap_plot <- function() {
	list(
		`Co-Exp` = c(
			ca_coexp_1p$feature_name_1 %>% unique,
			ca_coexp_1p$feature_name_2 %>% unique) %>% unique,
		`BlastP` = c(
			ca_blastp_5E$feature_name_1,
			ca_blastp_5E$feature_name_2) %>% unique,
		`PPI` = c(
			ca_biogrid$feature_name_1,
			ca_biogrid$feature_name_2) %>% unique,
		`Sac Genetic PPI` = c(
			sac_ortholog_ppi %>%
				dplyr::filter(experimental_system_type == "genetic") %>%
				magrittr::extract2("feature_name_1") %>%
				unique(),
			sac_ortholog_ppi %>%
				dplyr::filter(experimental_system_type == "genetic") %>%
				magrittr::extract2("feature_name_2") %>%
				unique) %>% unique,
		`Sac Physical PPI` = c(
			sac_ortholog_ppi %>%
				dplyr::filter(experimental_system_type == "physical") %>%
				magrittr::extract2("feature_name_1") %>%
				unique(),
			sac_ortholog_ppi %>%
				dplyr::filter(experimental_system_type == "physical") %>%
				magrittr::extract2("feature_name_2") %>%
				unique) %>% unique,
		`YeastNet` = c(
			yeast_net$feature_name_1,
			yeast_net$feature_name_2) %>% unique,
		`YeastNet-CC` = c(
			yeast_net_CC$feature_name_1,
			yeast_net_CC$feature_name_2) %>% unique,
		`YeastNet-CX` = c(
			yeast_net_CX$feature_name_1,
			yeast_net_CX$feature_name_2) %>% unique,
		`YeastNet-DC` = c(
			yeast_net_DC$feature_name_1,
			yeast_net_DC$feature_name_2) %>% unique,
		`YeastNet-GN` = c(
			yeast_net_GN$feature_name_1,
			yeast_net_GN$feature_name_2) %>% unique,
		`YeastNet-GT` = c(
			yeast_net_GT$feature_name_1,
			yeast_net_GT$feature_name_2) %>% unique,
		`YeastNet-HT` = c(
			yeast_net_HT$feature_name_1,
			yeast_net_HT$feature_name_2) %>% unique,
		`YeastNet-LC` = c(
			yeast_net_LC$feature_name_1,
			yeast_net_LC$feature_name_2) %>% unique,
		`YeastNet-PG` = c(
			yeast_net_PG$feature_name_1,
			yeast_net_PG$feature_name_2) %>% unique,
		`YeastNet-TS` = c(
			yeast_net_TS$feature_name_1,
			yeast_net_TS$feature_name_2) %>% unique,
		`CGD GO` = ca_go_propagated_filtered$feature_name %>% unique) %>%
		UpSetR::fromList() %>%
		UpSetR::upset(
			keep.order=TRUE,
			order.by="freq",
			point.size=3.5,
			mainbar.y.label="Candida Albicans Genes",
			sets.x.label="Network Nodes")
}

pdf(
	file="product/figures/network_overlap_20201118.pdf",
	width=9, height=5, useDingbats=FALSE, onefile=FALSE)
network_overlap_plot()
dev.off()

png(
	file="product/figures/network_overlap_20201118.png",
	width=360*4, height=360*4, res=72*4)
network_overlap_plot()
dev.off()
