

library(plyr)
library(dplyr)
library(hdf5r)

raw_path <- "raw_data/CoCoCoNet_20210104"

system(paste0("mkdir ", raw_path))

###############################
# candida Albicans CoCoCo Net #
###############################
system(
    paste0(
        "cd ", raw_path, " && ",
        "wget ftp://milton.cshl.edu/data/networks/candida_MetaAggNet.hdf5"))
system(
    paste0(
        "cd ", raw_path, " && ",
        "wget ftp://milton.cshl.edu/data/networks/candida_prioAggNet.hdf5"))

dataset <- hdf5r::H5File$new(paste0(raw_path, "/candida_prioAggNet.hdf5"))
dataset$ls(recursive = TRUE)

candida_prio_agg_genes <- dataset[["row"]][] %>%
    tibble::tibble(gene_name = .) %>%
    dplyr::mutate(
        feature_name = paste(
            gene_name %>% stringr::str_sub(8, 9),
            gene_name %>% stringr::str_sub(10, 15),
            gene_name %>% stringr::str_sub(16, 16),
            sep = "_"))

candida_prio_agg_net <- dataset[["agg"]][, ]

colnames(candida_prio_agg_net) <- candida_prio_agg_genes$feature_name
rownames(candida_prio_agg_net) <- candida_prio_agg_genes$feature_name
save(candida_prio_agg_net, file = "intermediate_data/CoCoCoNet_candida_prio_agg_net_20210104.Rdata")


system(
    paste0(
        "cd ", raw_path, " && ",
        "wget ftp://milton.cshl.edu/data/networks/yeast_MetaAggNet.hdf5"))
system(
    paste0(
        "cd ", raw_path, " && ",
        "wget ftp://milton.cshl.edu/data/networks/yeast_prioAggNet.hdf5"))

meta <- hdf5r::H5File$new(paste0(raw_path, "/candida_metaAggNet.hdf5"))
meta$ls(recursive = TRUE)

candida_meta_agg_genes <- meta[["row"]][] %>%
    tibble::tibble(gene_name = .) %>%
    dplyr::mutate(
        feature_name = paste(
            gene_name %>% stringr::str_sub(8, 9),
            gene_name %>% stringr::str_sub(10, 15),
            gene_name %>% stringr::str_sub(16, 16),
            sep = "_"))
candida_meta_agg_net <- meta[["agg"]][, ]
colnames(candida_meta_agg_net) <- candida_meta_agg_genes$feature_name
rownames(candida_meta_agg_net) <- candida_meta_agg_genes$feature_name
save(candida_meta_agg_net, file = "intermediate_data/CoCoCoNet_candida_meta_agg_net_20210104.Rdata")

