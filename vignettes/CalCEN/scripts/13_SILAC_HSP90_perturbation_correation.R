# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(ggplot2)
library(ggrepel)

load("intermediate_data/ca_silac_hsp90_fold_change.Rdata")


data <- ca_silac_hsp90_fold_change %>%
	dplyr::group_by(gene_label) %>%
	dplyr::summarize(
		sum_abs_ph_lfc = sum(abs(pharmacological_log_fc)),
		sum_abs_gt_lfc = sum(abs(genetic_log_fc)))

r_squared <- summary(lm(data$sum_abs_ph_lfc ~ data$sum_abs_gt_lfc))$r.squared


p <- ggplot() +
	theme_bw() +
	geom_point(
		data=data,
		mapping=aes(
			x=sum_abs_ph_lfc,
			y=sum_abs_gt_lfc),
		size=.8) +
	stat_smooth(
		data=data,
		mapping=aes(
			x=sum_abs_ph_lfc,
			y=sum_abs_gt_lfc),
		method = lm,
		color="blue") +
	 geom_label_repel(
		data=data %>% dplyr::filter(
			(sum_abs_ph_lfc > 30) |
			(sum_abs_gt_lfc > 30)),
		mapping=aes(
			x=sum_abs_ph_lfc,
			y=sum_abs_gt_lfc,
			label=gene_label),
		size=2) +
	geom_text(
			data=data.frame(sum_abs_ph_lfc=55, sum_abs_gt_lfc=5, label=paste0("R^2: ", signif(r_squared, 2))),
			aes(
				x=sum_abs_ph_lfc,
				y=sum_abs_gt_lfc,
				label=label),
			parse=TRUE) +
	ggtitle("Correlation of SILAC-fc after perturbation of HSP90") +
	scale_x_continuous("Pharmacological HSP90 inactivation\nsum(abs(log(intensity fold-change)") +
	scale_y_continuous("Genetic HSP90 depletion\nsum(abs(log(intensity fold-change)")

ggsave(
	"product/figures/silac_hsp90_perturb_correlation_180914.pdf",
	height=5, width=5,
	useDingbats=FALSE)
ggsave(
	"product/figures/silac_hsp90_perturb_correlation_180914.png",
	height=5, width=5)

