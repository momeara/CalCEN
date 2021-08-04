

library(plyr)
library(dplyr)
library(EGAD)

load("intermediate_data/ca_genes.Rdata")
load("intermediate_data/chromosome_features.Rdata")
#load("intermediate_data/sac_chromosome_features.Rdata")


load("intermediate_data/CoCoCoNet_candida_prio_agg_net_20210104.Rdata")
load("intermediate_data/CoCoCoNet_candida_meta_agg_net_20210104.Rdata")

CoCoCoNet_candida <- expand.grid(
    feature_name_1 = ca_genes,
    feature_name_2 = ca_genes) %>%
    dplyr::left_join(
        candida_prio_agg_net %>%
        as.data.frame() %>%
        dplyr::add_rownames(var = "feature_name_1") %>%
        tidyr::pivot_longer(
            cols = -"feature_name_1",
            names_to = "feature_name_2",
            values_to = "score"),
        by = c("feature_name_1", "feature_name_2")) %>%
    dplyr::mutate(score = ifelse(!is.na(score), score, 0)) %>%
    tidyr::pivot_wider(
        id_cols = "feature_name_1",
        names_from = "feature_name_2",
        values_from = "score")

CoCoCoNet_candida <- data.frame(CoCoCoNet_candida)
rownames(CoCoCoNet_candida) <- CoCoCoNet_candida$feature_name_1
CoCoCoNet_candida <- CoCoCoNet_candida %>%
    dplyr::select(-feature_name_1) %>%
    as.matrix()

save(CoCoCoNet_candida, file = "intermediate_data/CoCoCoNet_candida_20210104.Rdata")


##########
CoCoCoNet_meta_candida <- expand.grid(
    feature_name_1 = ca_genes,
    feature_name_2 = ca_genes) %>%
    dplyr::left_join(
        candida_meta_agg_net %>%
        as.data.frame() %>%
        dplyr::add_rownames(var = "feature_name_1") %>%
        tidyr::pivot_longer(
            cols = -"feature_name_1",
            names_to = "feature_name_2",
            values_to = "score"),
        by = c("feature_name_1", "feature_name_2")) %>%
    dplyr::mutate(score = ifelse(!is.na(score), score, 0)) %>%
    tidyr::pivot_wider(
        id_cols = "feature_name_1",
        names_from = "feature_name_2",
        values_from = "score")

CoCoCoNet_meta_candida <- data.frame(CoCoCoNet_meta_candida)
rownames(CoCoCoNet_meta_candida) <- CoCoCoNet_meta_candida$feature_name_1
CoCoCoNet_meta_candida <- CoCoCoNet_meta_candida %>%
    dplyr::select(-feature_name_1) %>%
    as.matrix()

save(CoCoCoNet_meta_candida, file = "intermediate_data/CoCoCoNet_meta_candida_20210104.Rdata")
