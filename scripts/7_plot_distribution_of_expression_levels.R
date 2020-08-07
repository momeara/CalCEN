# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(magrittr)
library(stringr)
library(readr)
library(ggplot2)
library(scales)

load("intermediate_data/estimated_expression.Rdata")

p <- ggplot() +
	theme_bw() +
	geom_density(
		data=estimated_expression,
		mapping=aes(x=log(FPKM), group=run_accession, color=study_accession),
		size=.5) +
	scale_color_discrete("Study Accession") +
	scale_x_continuous("log(FPKM)") +
	scale_y_continuous("Density")

ggplot2::ggsave(
	filename="product/figures/distribution_of_expression_levels.pdf",
	height=5, width=5)
ggplot2::ggsave(
	filename="product/figures/distribution_of_expression_levels.png",
	height=5, width=5)