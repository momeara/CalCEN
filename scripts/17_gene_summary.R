
library(plyr)
library(dplyr)
library(stringr)
library(readr)

load("intermediate_data/chromosome_features.Rdata")
load("intermediate_data/CalCEN_full.Rdata")


ca_summary <- function(query_feature_name, coexp_cutoff=.95){
  CalCEN_full %>%
    dplyr::select(
      feature_name = feature_name_1,
      coexp_feature_name = feature_name_2,
      coexp_score = score) %>%
    dplyr::filter(
      feature_name %in% query_feature_name,
      coexp_score >= coexp_cutoff,
      !is.na(coexp_feature_name),
      coexp_feature_name != query_feature_name) %>%
    dplyr::left_join(
      chromosome_features %>%
        dplyr::select(
          coexp_feature_name = feature_name,
          coexp_primary_cgd_id = primary_cgd_id,
          coexp_description = description),
        by=c("coexp_feature_name")) %>%
    dplyr::arrange(feature_name, desc(coexp_score))
}



b <- ca_summary(c("C7_00310C_A", "C1_11670W_A"))

b %>% readr::write_tsv("product/coexp_summary_191107.tsv")
