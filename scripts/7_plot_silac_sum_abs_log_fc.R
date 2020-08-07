# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(stringr)
library(readr)
library(magrittr)
library(ggplot2)
library(ggrepel)
load("intermediate_data/ca_silac_hsp90_fold_change.Rdata")


data <-	rbind(
	ca_silac_hsp90_fold_change %>%
		dplyr::filter(!is.na(gene_label)) %>%
		dplyr::group_by(gene_label) %>%
			dplyr::summarize(
				uniprot_accn=uniprot_accn[1],
				feature_name=feature_name[1],
				expression = sum(wildtype),
				delta = sum(abs(pharmacological_log_fc))) %>%
		dplyr::mutate(condition="pharmacological"),
	ca_silac_hsp90_fold_change %>%
		dplyr::filter(!is.na(gene_label)) %>%
		dplyr::group_by(gene_label) %>%
			dplyr::summarize(
				uniprot_accn=uniprot_accn[1],
				feature_name=feature_name[1],
				expression = sum(wildtype),
				delta = sum(abs(genetic_log_fc))) %>%
		dplyr::mutate(condition="genetic"))

silac_genes_of_interest <- data %>%
	dplyr::select(-expression) %>%
	dplyr::semi_join(
		data %>%
			dplyr::filter(delta > 25),
		by="gene_label") %>%
	tidyr::spread("condition", "delta") %>%
	dplyr::mutate(note=paste0("G",signif(genetic,3), "|P", signif(pharmacological,3)))

silac_genes_of_interest %>%
	readr::write_tsv("product/silac_genes_of_interest.tsv")


p <- ggplot2::ggplot() +
	theme_bw() +
	geom_point(
 		data=data,
		mapping=aes(
			x=expression,
			y=delta)) +
	geom_label_repel(
		data=data %>% dplyr::filter(delta > 25),
		mapping=aes(
			x=expression,
			y=delta,
			label=gene_label),
		force=2,
		max.iter=20000,
#		ylim=c(30, NA),
#		xlim=c(1.2e11,1e14),
		size=2) +
	facet_wrap(~condition, nrow=1) +
#	coord_fixed(ratio = 1) +
	scale_x_log10("Total Intensity") +
	scale_y_continuous("Sum of the Absolute Log Fold Change") +
	scale_color_discrete(guide=FALSE)

ggsave(
	"product/figures/silac_sum_abs_log_fc.pdf",
	height=4, width=8,
	useDingbats=FALSE)
ggsave(
	"product/figures/silac_sum_abs_log_fc.png",
	height=4, width=8)

