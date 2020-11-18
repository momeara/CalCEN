# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(ggplot2)

gba_summary <- readr::read_tsv(
	file="product/gba_summary_10f_C-B-SP-SG-YN_20201113.tsv",
	col_types=readr::cols(
	  anno_id = readr::col_character(),
	  network_id = readr::col_character(),
	  auroc_mean = readr::col_double(),
	  auroc_std = readr::col_double()))

auroc_data <- gba_summary %>%
	dplyr::mutate(networks = stringr::str_split(network_id, "[|]")) %>%
	tidyr::unnest(networks) %>%
	dplyr::mutate(value=1) %>%
	tidyr::spread("networks", value, 0) %>%
	dplyr::mutate(
		degree = dplyr::select(., -anno_id, -network_id, -auroc_mean, -auroc_std) %>% rowSums()) %>%
	dplyr::arrange(degree, network_id) %>%
	dplyr::group_by(anno_id) %>%
	dplyr::mutate(set_id = row_number()) %>%
	dplyr::ungroup()


sets_1hot <- auroc_data %>%
	dplyr::distinct(set_id, .keep_all=TRUE) %>%
	dplyr::select(-anno_id, -auroc_mean, -auroc_std, -network_id, -degree) %>%
	as.data.frame() %>%
	tibble::column_to_rownames("set_id") %>%
	t

	# this adapts the bottom panel from UpSetR

	name_size_scale <- 2
	shade_alpha <- 0.2
	shade_color <- "gray80"
	point_size <- 5
	line_size <- 2
	dot_color <- "gray10"
	dot_alpha <- 1
	n_sets <- ncol(sets_1hot)
	n_items <- nrow(sets_1hot)

	sets_data <- expand.grid(
		y = 1:nrow(sets_1hot),
		x = 1:ncol(sets_1hot)) %>%
		dplyr::mutate(
			value = as.vector(sets_1hot),
			color = ifelse(value > 0L, dot_color, "gray83"),
			alpha = ifelse(value > 0L, 1, dot_alpha),
			intersection = ifelse(
				value > 0L,
				paste0(x, "yes"),
				paste0(row_number(), "no")))


	upper_sets_shading_data <- sets_data %>%
		dplyr::distinct(x) %>%
		dplyr::filter(x %% 2 != 0) %>%
		dplyr::mutate(
			ymin=0.5,
			ymax=1,
			xmin=x-0.5,
			xmax=x+0.5,
			color=shade_color)

	lower_sets_shading_data <- sets_data %>%
		dplyr::distinct(x) %>%
		dplyr::filter(x %% 2 != 0) %>%
		dplyr::mutate(
			ymin=0 + 0.5,
			ymax=n_items+0.5,
			xmin=x-0.5,
			xmax=x+0.5,
			color=shade_color)

	item_shading_data <- sets_data %>%
		dplyr::distinct(y) %>%
		dplyr::filter(y %% 2 != 0) %>%
		dplyr::mutate(
			xmin=0 + 0.5,
			xmax=max(sets_data$x) + 0.5,
			ymin=y-0.5,
			ymax=y+0.5,
			color=shade_color)

	dots_plot <- ggplot() +
		theme(
			legend.position = c(0.5, 0.9),
			legend.direction = "horizontal",
			panel.background = element_rect(fill = "white"),
	    plot.margin=unit(c(-0.2,0.5,-.5,0.5), "lines"),
	    panel.grid.major.x = element_blank(),
	    panel.grid.minor.x = element_blank(),
	    panel.grid.major.y = element_line(colour="grey80", size=0.5),
	    panel.grid.minor.y = element_line(colour="grey90", size=0.5),
	    axis.text.x = element_blank(),
	    axis.ticks.x = element_blank(),
	    axis.ticks.y = element_blank(),
	    axis.text.y = element_text(
				colour = "gray0", size = 7*name_size_scale, hjust = 0.4),
			axis.title.y = element_text(
			  colour = "gray0", size = 8*name_size_scale,
				margin = margin(t = 0, r = -40, b = 0, l = 0))) +
			xlab(NULL) +
			scale_y_continuous(
				"Area Under the ROC Curve",
				limits=c(0.5, 1)) +
	    scale_x_continuous(
				limits=c(0, n_sets+1),
				expand = c(0, 0)) +
			scale_color_discrete("GO Ontology") +
	    geom_rect(
				data = upper_sets_shading_data,
				aes(
					xmin = xmin, xmax = xmax,
	        ymin = ymin, ymax = ymax),
	        fill = shade_color,
					alpha = shade_alpha) +
#		geom_line(
#			data=auroc_data %>% dplyr::filter(anno_id != 'all'),
#			aes(x=set_id, y=auroc_mean, color=anno_id, group=degree),
#			size=1) +
##		geom_errorbar(
##			data=auroc_data %>% dplyr::filter(anno_id != 'all'),
##			aes(
##				x=set_id,
##				ymin=auroc_mean-auroc_std,
##				ymax=auroc_mean+auroc_std,
##				color=anno_id),
##			size=.5) +
		geom_point(
			data=auroc_data %>% dplyr::filter(anno_id != 'all'),
			aes(x=set_id, y=auroc_mean, color=anno_id),
			size=2) +
#		geom_line(
#			data=auroc_data %>% dplyr::filter(anno_id == 'all'),
#			aes(x=set_id, y=auroc_mean, group=degree),
#			size=1.5) +
		geom_errorbar(
			data=auroc_data %>% dplyr::filter(anno_id == 'all'),
			aes(x=set_id, ymin=auroc_mean-auroc_std/sqrt(10), ymax=auroc_mean+auroc_std/sqrt(10)),
			width = .05,
			size=.2) +
		geom_point(
			data=auroc_data %>% dplyr::filter(anno_id == 'all'),
			aes(x=set_id, y=auroc_mean),
			size=3)


	sets_plot <- ggplot() +
		theme(
			panel.background = element_rect(fill = "white"),
	    plot.margin=unit(c(-0.5,0.5,0.5,0.5), "lines"),
	    axis.text.x = element_blank(),
	    axis.ticks.x = element_blank(),
	    axis.ticks.y = element_blank(),
	    axis.text.y = element_text(
				colour = "gray0", size = 7*name_size_scale, hjust = 1),
	    panel.grid.major = element_blank(),
	    panel.grid.minor = element_blank()) +
			xlab(NULL) + ylab("   ") +
	    scale_y_continuous(
				breaks = c(1:n_items),
	      limits = c(0.5, n_items+0.5),
	      labels = rownames(sets_1hot),
				expand = c(0,0)) +
	    scale_x_continuous(
				limits = c(0, n_sets+1),
				expand = c(0, 0)) +
	    geom_rect(
				data = item_shading_data,
				aes(
					xmin = xmin, xmax = xmax,
	        ymin = ymin, ymax = ymax),
	        fill = shade_color,
					alpha = shade_alpha) +
	    geom_rect(
				data = lower_sets_shading_data,
				aes(
					xmin = xmin, xmax = xmax,
	        ymin = ymin, ymax = ymax),
	        fill = shade_color,
					alpha = shade_alpha) +
	    geom_point(
				data=sets_data,
				aes(x=x, y=y),
				colour = sets_data$color,
	      size= point_size,
				alpha = sets_data$alpha,
				shape=16) +
			geom_line(
				data= sets_data,
				aes(x=x, y=y,
					group = intersection,
	        colour=color),
				size = line_size) +
	    scale_color_identity()

dots_plot <- dots_plot %>% ggplotGrob
sets_plot <- sets_plot %>% ggplotGrob

dots_plot$widths[2:5] <- sets_plot$widths[2:5]

p <- gridExtra::grid.arrange(
	dots_plot,
	sets_plot,
	heights=c(4.5, 1.5),
	ncol=1)

ggsave(
	file="product/figures/go_pred_SubO_10f_C-B-SP-SG_YN_20201117.pdf",
	plot=p,
	width=7.5, height=5,
	useDingbats=FALSE)

ggsave(
	file="product/figures/go_pred_SubO_10f_C-B-SP-SG_Y20201117.png",
	plot=p,
	width=7.5, height=5)
