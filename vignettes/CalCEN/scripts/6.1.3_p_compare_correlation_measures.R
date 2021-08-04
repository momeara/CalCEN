# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(tidyverse)

load("intermediate_data/CalCEN_full_spearman.Rdata")
load("intermediate_data/CalCEN_full_pearson.Rdata")
load("intermediate_data/CalCEN_full_direct.Rdata")

########################################################
# compare spearman vs pearson correlation coefficients #
########################################################

data <- dplyr::bind_cols(
    CalCEN_full_spearman %>% dplyr::rename(score_spearman = score),
    CalCEN_full_pearson %>% dplyr::rename(score_pearson = score))

ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::geom_point(
        data = data %>% dplyr::sample_n(100000),
        mapping = ggplot2::aes(
            x = score_spearman,
            y = score_pearson),
        size = .1,
        shape=16,
        alpha = .4) +
    ggplot2::ggtitle(
        label = "Pearson vs. Spearman expression correlation coefficients") +
    ggplot2::scale_x_continuous(
        "Spearman Rank Correlation") +
    ggplot2::scale_y_continuous(
        "Pearson correlation Coefficient")

ggplot2::ggsave(
		filename="product/figures/pearson_vs_spearman_20201024.pdf",
		height=8, width=8,
		useDingbats=FALSE)

ggplot2::ggsave(
		filename="product/figures/pearson_vs_spearman_20201024.png",
		height=8, width=8)


################################
# direct correlation threshold #
################################

load("intermediate_data/CalCEN_full_direct.Rdata")

# following broom:::tidy.dgTMatrix which is depricated, but I'm not sure why
# assume rownames and column names are defined
tidy_sparse_matrix <- function(matrix){
		s <- Matrix::summary(matrix)
		tibble::tibble(
				row = rownames(matrix)[s$i],
				column = colnames(matrix)[s$j],
				value = s$x)
}

for (i in 1:30){
		z <- CalCEN_full_direct$beta[[i]] %>% tidy_sparse_matrix()
		cat("for path value '", i, "', n edges: ", nrow(z), "\n", sep = "")
}

load("intermediate_data/CalCEN_full_direct_bootstrap.Rdata")
