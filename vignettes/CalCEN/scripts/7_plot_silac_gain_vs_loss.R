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


fc_sig <- 1

data <-	rbind(
	ca_silac_hsp90_fold_change %>%
		dplyr::group_by(gene_label) %>%
			dplyr::summarize(
				gain = sum(max(pharmacological_log_fc, 0)),
				loss = -sum(min(pharmacological_log_fc, 0))) %>%
		dplyr::mutate(condition="pharmacological"),
	ca_silac_hsp90_fold_change %>%
		dplyr::group_by(gene_label) %>%
			dplyr::summarize(
				gain = sum(max(genetic_log_fc, 0)),
				loss = -sum(min(genetic_log_fc, 0))) %>%
		dplyr::mutate(condition="genetic")) %>%
	dplyr::mutate(
		gain_sig = gain > fc_sig,
		loss_sig = loss > fc_sig)

p <- ggplot2::ggplot() +
	theme_bw() +
	geom_vline(xintercept=fc_sig, size=.5) +
	geom_hline(yintercept=fc_sig, size=.5) +

	geom_point(
 		data=data,
		mapping=aes(
			x=loss,
			y=gain,
			color=interaction(gain_sig, loss_sig))) +
	geom_label_repel(
		data=data %>% dplyr::filter(
			gain > 2.5*fc_sig |
			loss > 2.5*fc_sig |
			(gain > 1.2*fc_sig & loss > 1.2*fc_sig)),
		mapping=aes(
			x=loss,
			y=gain,
			label=gene_label,
			color=interaction(gain_sig, loss_sig)),
		size=2) +
	facet_wrap(~condition, nrow=1) +
	coord_fixed(ratio = 1) +
	scale_x_continuous("Protein Lost with HSP90 Deplection") +
	scale_y_continuous("Protein Gained with HSP90 Depletion") +
	scale_color_discrete(guide=FALSE)

ggsave(
	"product/figures/silac_gain_vs_loss.pdf",
	height=4, width=8,
	useDingbats=FALSE)
ggsave(
	"product/figures/silac_gain_vs_loss.png",
	height=4, width=8)

