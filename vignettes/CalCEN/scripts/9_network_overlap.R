g# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(UpSetR)

load("intermediate_data/chromosome_features.Rdata")
load("intermediate_data/estimated_expression.Rdata")
load("intermediate_data/CalCEN.Rdata")
load("intermediate_data/ca_blastp.Rdata")
load("intermediate_data/ca_biogrid.Rdata")
load("intermediate_data/sac_biogrid.Rdata")


load("intermediate_data/yeast_net_network.Rdata")

load("intermediate_data/ca_go_propagated_filtered.Rdata")




CalCEN_1p <- CalCEN %>%
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

yeast_net_long <- yeast_net_network %>%
		as.data.frame() %>%
		tibble::rownames_to_column(var = "feature_name_1") %>%
		tidyr::pivot_longer(
				cols=-feature_name_1,
				names_to = "feature_name_2",
				values_to = "score") %>%
		dplyr::filter(score > 0)


network_overlap_plot <- function() {
	list(
		`Co-Exp` = c(
			CalCEN_1p$feature_name_1 %>% unique,
			CalCEN_1p$feature_name_2 %>% unique) %>% unique,
		`BlastP` = c(
			ca_blastp_5E$feature_name_1,
			ca_blastp_5E$feature_name_2) %>% unique,
		`PPI` = c(
			ca_biogrid$feature_name_1,
			ca_biogrid$feature_name_2) %>% unique,
		`SacGene` = c(
			sac_ortholog_ppi %>%
				dplyr::filter(experimental_system_type == "genetic") %>%
				magrittr::extract2("feature_name_1") %>%
				unique(),
			sac_ortholog_ppi %>%
				dplyr::filter(experimental_system_type == "genetic") %>%
				magrittr::extract2("feature_name_2") %>%
				unique) %>% unique,
		`SacPhys` = c(
			sac_ortholog_ppi %>%
				dplyr::filter(experimental_system_type == "physical") %>%
				magrittr::extract2("feature_name_1") %>%
				unique(),
			sac_ortholog_ppi %>%
				dplyr::filter(experimental_system_type == "physical") %>%
				magrittr::extract2("feature_name_2") %>%
				unique) %>% unique,
		`YeastNet` = c(
			yeast_net_long$feature_name_1,
			yeast_net_long$feature_name_2) %>% unique) %>%
#		`CGD GO` = ca_go_propagated_filtered$feature_name %>% unique) %>%
		UpSetR::fromList() %>%
			UpSetR::upset(
					sets = c(
							"YeastNet",
							"SacPhys",
							"SacGene",
							"Co-Exp",
							"BlastP"),
			keep.order = TRUE,
			order.by = "freq",
			point.size = 3.5,
			set_size.show = FALSE,
			mainbar.y.label = "Candida albicans genes",
			sets.x.label = "Network Nodes")
}

pdf(
	file="product/figures/network_overlap_20201124.pdf",
	width=6, height=5, useDingbats=FALSE, onefile=FALSE)
network_overlap_plot()
dev.off()

png(
	file="product/figures/network_overlap_20201124.png",
	width=2400*4, height=360*4, res=72*4)
network_overlap_plot()
dev.off()
