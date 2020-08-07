# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(stringr)
library(readr)
library(magrittr)
library(ggplot2)
library(ggrepel)


load("intermediate_data/silac.Rdata")

p <- ggplot2::ggplot(data=silac) +
	theme_bw() +
	geom_density(
		mapping=aes(x=log(fold_change), color=depletion)) +
	scale_x_continuous("Log(fold change)") +
	ggtitle("SILAC protein-cluster scores")

ggsave(
	"product/figures/silac_fold_change.pdf",
	height=3, width=5)
ggsave(
	"product/figures/silac_fold_change.png",
	height=3, width=5)


fc_sig <- 1.5

mean_fc <- silac %>%
	dplyr::filter(!is.na(fold_change)) %>%
	dplyr::filter(!is.na(feature_name)) %>%
	dplyr::mutate(fold_change = 1/fold_change) %>%
	dplyr::group_by(name, depletion) %>%
	dplyr::summarize(mean_fc = mean(fold_change)) %>%
	tidyr::spread(depletion, mean_fc) %>%
	dplyr::mutate(
		genetic_sig = genetic > fc_sig | genetic < 1/fc_sig,
		pharmacological_sig = pharmacological > fc_sig | pharmacological < 1/fc_sig)

p <- ggplot2::ggplot() +
	theme_bw() +
	geom_vline(xintercept=log(fc_sig), size=.5) +
	geom_vline(xintercept=-log(fc_sig), size=.5) +
	geom_hline(yintercept=log(fc_sig), size=.5) +
	geom_hline(yintercept=-log(fc_sig), size=.5) +

	geom_point(
  		data=mean_fc,
		mapping=aes(
			x=log(genetic),
			y=log(pharmacological),
			color=interaction(genetic_sig, pharmacological_sig))) +
	geom_label_repel(
		data=mean_fc %>% dplyr::filter(
				genetic > 3*fc_sig | genetic < 1/(3*fc_sig) |
				pharmacological > 3*fc_sig | pharmacological < 1/(3*fc_sig)),
		mapping=aes(
			x=log(genetic), y=log(pharmacological),
			label=name,
			color=interaction(genetic_sig, pharmacological_sig)),
		size=2) +
	scale_x_continuous("log(Genetic /  WT) Spectral Counts") +
	scale_y_continuous("log(Pharamcological / WT) Spectral Counts") +
	scale_color_discrete(guide=FALSE) +
	ggtitle("SILAC mean fold change for HSP90 depletion vs. WT")

ggsave(
	"product/figures/silac_genetic_vs_pharmacolgocal_mean_fold_change.pdf",
	height=7, width=7,
	useDingbats=FALSE)
ggsave(
	"product/figures/silac_genetic_vs_pharmacolgocal_mean_fold_change.png",
	height=7, width=7)


#####################33
# version 2

load("intermediate_data/ca_silac_hsp90_fold_change.Rdata")

lfc <- ca_silac_hsp90_fold_change %>%
		dplyr::group_by(gene_label) %>%
		dplyr::summarize(
				total_wt = sum(wildtype),
				total_ph = sum(pharmacological),
				total_gt = sum(genetic),
				ph_lfc = log(total_ph/total_wt),
				gt_lfc = log(total_gt/total_wt)) %>%
		dplyr::filter(
				total_wt > 0) %>%
		dplyr::filter(
				is.finite(ph_lfc),
				is.finite(gt_lfc))

gt_lm <- lm(
	data= lfc %>% dplyr::transmute(x=log(total_wt),y=log(total_gt)),
	formula = y ~ x)
gt_r2 <- summary(gt_lm)$r.squared

ph_lm <- lm(
	data= lfc %>% dplyr::transmute(x=log(total_wt),y=log(total_ph)),
	formula = y ~ x)
ph_r2 <- summary(ph_lm)$r.squared



p <- ggplot2::ggplot() +
	theme_bw() +
	geom_smooth(
		data=lfc %>%
			dplyr::transmute(
				x=log(total_wt),
				y=log(total_gt),
				treatment="Genetic"),
		mapping=aes(x=x, y=y),
		method='lm',
		se = TRUE,
		formula=y~x) +
	geom_smooth(
		data=lfc %>%
			dplyr::transmute(
				x=log(total_wt),
				y=log(total_ph),
				treatment="Pharmacological"),
		mapping=aes(x=x, y=y),
		method='lm',
		se = TRUE,
		formula=y~x) +
	geom_point(
 		data=lfc %>% dplyr::mutate(treatment="Genetic"),
		mapping=aes(
			x=log(total_wt),
			y=log(total_gt))) +
	geom_point(
 		data=lfc %>% dplyr::mutate(treatment="Pharmacological"),
		mapping=aes(
			x=log(total_wt),
			y=log(total_ph))) +
	geom_label_repel(
		data=lfc %>%
			dplyr::mutate(treatment="Genetic") %>%
			dplyr::filter(abs(gt_lfc) > 2),
		mapping=aes(
			x=log(total_wt), y=log(total_gt),
			label=gene_label),
		size=2) +
	geom_label_repel(
		data=lfc %>%
			dplyr::mutate(treatment="Pharmacological") %>%
			dplyr::filter(abs(ph_lfc) > 2),
		mapping=aes(
			x=log(total_wt), y=log(total_ph),
			label=gene_label),
		size=2) +
	facet_wrap(~treatment) +
	geom_label(
		data=data.frame(
			treatment=c("Genetic", "Pharmacological"),
			x = c(25, 25),
			y = c(14, 14),
			label=paste0("R2 ", c(signif(gt_r2, 2), signif(ph_r2, 2)))),
		mapping=aes(x=x, y=y, label=label)) +
	scale_x_continuous("Wildtype log(Intensity)") +
	scale_y_continuous("Treatment Log(Intensity)") +
	ggtitle("SILAC Treatment vs Wildtype intensity")

ggsave(
	"product/figures/silac_treatment_vs_wilttype_log_intensity_181110.pdf",
	height=5, width=7,
	useDingbats=FALSE)
ggsave(
	"product/figures/silac_treatment_vs_wilttype_log_intensity_181110.png",
	height=5, width=7)



##########################

robust_fit <- MASS::rlm(lfc$gt_lfc ~ lfc$ph_lfc)
r_squared <- summary(lm(lfc$gt_lfc ~ lfc$ph_lfc, weights=robust_fit$w))$r.squared

SSe <- sum((robust_fit$resid)^2)
observed <- robust_fit$resid+robust_fit$fitted
SSt <- sum((robust_fit$w*observed-mean(robust_fit$w*observed))^2)
r_squared <- 1 - SSe/SSt


p <- ggplot2::ggplot() +
	theme_bw() +
	geom_vline(xintercept=0, size=.5) +
	geom_hline(yintercept=0, size=.5) +
	geom_point(
 		data=lfc,
		mapping=aes(
			x=ph_lfc,
			y=gt_lfc)) +
	stat_smooth(
		data=lfc,
		mapping=aes(
			x=ph_lfc,
			y=gt_lfc),
		method = rlm,
		color="blue") +
	geom_text(
			data=data.frame(ph_lfc=1.5, gt_lfc=6, label=paste0("R^2: ", signif(r_squared, 2))),
			aes(
				x=ph_lfc,
				y=gt_lfc,
				label=label),
			parse=TRUE) +
	geom_label_repel(
		data=lfc %>% dplyr::filter(
				(abs(ph_lfc) > 2) |
				(abs(gt_lfc) > 2)),
		mapping=aes(
			x=ph_lfc, y=gt_lfc,
			label=gene_label),
		size=2) +
	scale_x_continuous("log(Pharamcological / WT) Intensity") +
	scale_y_continuous("log(Genetic /  WT) Intensity") +
	ggtitle("SILAC mean fold change for HSP90 depletion vs. WT")

ggsave(
	"product/figures/silac_genetic_vs_pharmacolgocal_mean_fold_change_180914.pdf",
	height=5, width=5,
	useDingbats=FALSE)
ggsave(
	"product/figures/silac_genetic_vs_pharmacolgocal_mean_fold_change_180914.png",
	height=5, width=5)
