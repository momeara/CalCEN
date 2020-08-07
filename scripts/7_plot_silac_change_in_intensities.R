# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(stringr)
library(readr)
library(magrittr)
library(ggplot2)
library(ggrepel)
load("intermediate_data/ca_silac_hsp90_intensities.Rdata")

load("intermediate_data/genes_of_interest.Rdata")

data <- ca_silac_hsp90_intensities %>%
	dplyr::filter(gene_label %in% genes_of_interest) %>%
	dplyr::mutate(gene_label = as.factor(gene_label)) %>%
	dplyr::select(gene_label, condition, fraction_index, intensity)




p <- ggplot2::ggplot(data=data) +
	theme_bw() +
	geom_line(
		mapping=aes(x=fraction_index, y=intensity, color=condition)) +
	facet_wrap(
		~gene_label,
		scales="free_y",
		ncol=1,
		strip.position = "right") +
	scale_y_continuous("Intensity") +
	scale_x_discrete("Elution Fraction") +
	scale_color_manual(
		"Condition",
		values=c(
			"wildtype"="black",
			"pharmacological"="blue",
			"genetic"="orange")) +
	ggtitle("SILAC Elution After HSP90 Depletion")

ggsave(
	file="product/figures/silac_change_in_intensities.pdf",
	height=7,
	width=10)
ggsave(
	file="product/figures/silac_change_in_intensities.png",
	height=7,
	width=10)


#####
normed_intensity_threshold <- .01
vertical_spread <- .66
sig_threshold <- 1.5
min_frac <- .05
genes_of_interest %>%
	dplyr::filter(include) %>%
	plyr::d_ply(c("set"), function(genes){
	set_id <- genes$set[1]
	cat("plotting change in intensities for ", set_id, " ...\n", sep="")
	data <- genes %>%
		dplyr::select(feature_name, gene_name) %>%
		dplyr::left_join(ca_silac_hsp90_intensities, by=c("feature_name")) %>%
		dplyr::mutate(gene_label = ifelse(!is.na(gene_label), gene_label, gene_name)) %>%
		dplyr::select(gene_label, condition, fraction_index, intensity) %>%
		dplyr::group_by(gene_label) %>%
		#dplyr::mutate(normed_intensity = log(intensity+1)) %>%
		dplyr::mutate(normed_intensity = intensity / max(ifelse(condition == "wildtype", intensity, 0))) %>%
		dplyr::ungroup()
	sig_candidate <- data %>%
		dplyr::select(-intensity) %>%
		tidyr::spread(condition, normed_intensity) %>%
		dplyr::mutate(
			wildtype_sig_frac = wildtype > min_frac,
			gene_sig_frac = genetic > min_frac,
			pharm_sig_frac = pharmacological > min_frac) %>%
		dplyr::select(gene_label, fraction_index, wildtype_sig_frac, gene_sig_frac, pharm_sig_frac)
	sig_diff <- data %>%
		dplyr::select(-normed_intensity) %>%
		tidyr::spread(condition, intensity) %>%
		dplyr::left_join(sig_candidate, by=c("gene_label", "fraction_index")) %>%
		dplyr::mutate(
			pharm_diff = (
				((pharm_sig_frac) &
					((wildtype == 0) | ((wildtype > 0) & (pharmacological/wildtype >= sig_threshold)))) |
				((wildtype_sig_frac) &
					((pharmacological == 0) | ((pharmacological > 0) & (wildtype/pharmacological >= sig_threshold))))),
			gene_diff = (
				((gene_sig_frac) &
					((wildtype == 0) | ((wildtype > 0) & (genetic/wildtype >= sig_threshold)))) |
				((wildtype_sig_frac) &
					((genetic == 0) | ((genetic > 0) & (wildtype/genetic >= sig_threshold)))))) %>%
#		dplyr::mutate(
#			pharm_diff = (
#				((wildtype == 0) & (pharmacological > 0)) |
#				((wildtype > 0) & (pharmacological/wildtype >= sig_threshold)) |
#				((pharmacological == 0) & (wildtype > 0)) |
#				((pharmacological > 0) & (wildtype/pharmacological >= sig_threshold))),
#			gene_diff = (
#				((wildtype == 0) & (genetic > 0)) |
#				((wildtype > 0) & (genetic/wildtype >= sig_threshold)) |
#				((genetic == 0) & (wildtype > 0)) |
#				((genetic > 0) & (wildtype/genetic >= sig_threshold)))) %>%
		dplyr::select(
			gene_label, fraction_index, wildtype, genetic, pharmacological, pharm_diff, gene_diff)
	first_fraction_index <- data %>%
		dplyr::filter(normed_intensity > normed_intensity_threshold) %>%
		dplyr::arrange(fraction_index) %>%
		dplyr::slice(1) %>%
		magrittr::extract2("fraction_index")
	last_fraction_index <- data %>%
		dplyr::filter(normed_intensity > normed_intensity_threshold) %>%
		dplyr::arrange(desc(fraction_index)) %>%
		dplyr::slice(1) %>%
		magrittr::extract2("fraction_index")
	data <- data %>%
		dplyr::filter(
			first_fraction_index <= fraction_index,
			fraction_index <= last_fraction_index)
	sig_diff <- sig_diff %>%
		dplyr::filter(
			first_fraction_index <= fraction_index,
			fraction_index <= last_fraction_index)
	gene_labels <- data %>%
		dplyr::distinct(gene_label) %>%
		dplyr::mutate(
			gene_index = n() - dplyr::row_number() + 1,
			gene_y = gene_index*vertical_spread)
	data <- data %>% dplyr::left_join(gene_labels, by=c("gene_label"))
	sig_diff <- sig_diff %>% dplyr::left_join(gene_labels, by=c("gene_label"))
	p <- ggplot2::ggplot(data=data) +
		theme_bw() +
		geom_rect(
			data=sig_diff,
			mapping=aes(
				xmin=fraction_index-.5, xmax=fraction_index+.5,
				ymin=gene_y+vertical_spread/2, ymax=gene_y+vertical_spread,
				alpha=pharm_diff*.3,
				fill="pharmacological")) +
		geom_rect(
			data=sig_diff,
			mapping=aes(
				xmin=fraction_index-.5, xmax=fraction_index+.5,
				ymin=gene_y, ymax=gene_y+vertical_spread/2,
				alpha=gene_diff*.3,
				fill="genetic")) +
		geom_line(
			mapping=aes(
				x=fraction_index,
				y=normed_intensity+gene_y,
				color=condition,
				group=interaction(condition, gene_label))) +
		scale_alpha(guide='none') +
		scale_y_continuous(
			"Normalized Intensity",
			breaks=gene_labels$gene_y,
			labels=gene_labels$gene_label) +
		scale_x_discrete("Elution Fraction") +
		scale_fill_manual(
			paste0(sig_threshold, "-fold change"),
			values=c(
				"pharmacological"="lightblue",
				"genetic"="orange")) +
		scale_color_manual(
			"Condition",
			values=c(
				"wildtype"="black",
				"pharmacological"="blue",
				"genetic"="red")) +
		ggtitle(paste0("SILAC elution of ", set_id, " proteins after HSP90 Depletion")) +
		theme(
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			panel.background = element_blank(),
			axis.line = element_line(colour = "black"),
			axis.ticks = element_blank(),
			axis.text.x = element_blank())
	ggsave(
		file=paste0("product/figures/silac_change_in_intensities_overlap_", set_id, "_sig_threshold_",sig_threshold,".pdf"),
		height=length(gene_labels),
		width=10)
	ggsave(
		file=paste0("product/figures/silac_change_in_intensities_overlap_", set_id, "_sig_threshold_",sig_threshold,".png"),
		height=length(gene_labels),
		width=10)
})
