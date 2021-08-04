

library(plyr)
library(tidyverse, quietly = TRUE, warn.conflicts = FALSE)
library(monocle3, quietly = TRUE, warn.conflicts = FALSE)


load("intermediate_data/ca_runs_final.Rdata")
load("intermediate_data/estimated_expression.Rdata")
load("product/CalCEN_network_full_20201007.Rdata")
date_code <- "20201113"
date_code <- "20210310"
date_code <- "20210319"

estimated_expression <- estimated_expression %>%
		dplyr::semi_join(ca_runs_final, by = c("run_accession"))

estimated_expression <- estimated_expression %>%
    dplyr::mutate(FPKM_normed = log(FPKM + 1))


if (! dir.exists(paths = "product/figures/CalCEN_embedding")) {
    cat("Creating 'product/figures/CalCEN_embedding'\n", sep = "")
    dir.create(path = "product/figures/CalCEN_embedding", recursive = TRUE)
}



#####################
# Embed experiments #
#####################

# translating to monocle terminology
# "gene" <- feature_name a.k.a CGD gene id
# "cell" <- run_accession

exprs <- reshape2::acast(
    data=estimated_expression,
    formula=gene_id ~ run_accession,
    value.var="FPKM")

gene_metadata <- data.frame(gene_short_name = exprs %>% colnames())

cell_metadata <-  data.frame(run_accession = exprs %>% colnames) %>%
    dplyr::left_join(
        ca_runs_final %>%
            dplyr::rename(experient_sample_name = sample_name),
        by = "run_accession") %>%
    dplyr::mutate(
        study_accession = as.factor(study_accession))

row.names(gene_metadata) <- gene_metadata$gene_short_name
row.names(exprs) <- gene_metadata$gene_short_name
row.names(cell_metadata) <- cell_metadata$run_accession

coexp_cds <- monocle3::new_cell_data_set(
    expression_data = exprs,
    gene_metadata = gene_metadata,
    cell_metadata = cell_metadata)

coexp_cds <- coexp_cds %>% monocle3::preprocess_cds(
    num_dims = 100)

coexp_cds <- coexp_cds %>%
    monocle3::reduce_dimension(
        preprocess_method = "PCA",
        umap.min_dis = .5,
        umap.n_neighbors = 30L,
        verbose = TRUE)


coexp_cds <- coexp_cds %>%
    monocle3::cluster_cells(
        k = 30,
        num_iter = 10,
        resolution = .1,
        verbose = TRUE)

# id like to label the graph by study accession...
clusters <- cell_metadata$study_accession
names(clusters) <- cell_metadata$run_accession
coexp_cds@clusters[["UMAP"]]$clusters <- clusters

coexp_cds %>%
    monocle3::plot_cells(
        show_trajectory_graph = FALSE,
        cell_size = 0.8)

if (! dir.exists(paths = "product/figures/CalCEN_embedding")) {
    cat("Creating 'product/figures/CalCEN_embedding'\n", sep = "")
    dir.create(path = "product/figures/CalCEN_embedding", recursive = TRUE)
}

ggplot2::ggsave(
    filename = "product/figures/CalCEN_embedding/UMAP_experiments_20201113.pdf",
    width = 6,
    height = 6)


clusters <- coexp_cds@clusters[["UMAP"]]$clusters %>%
    data.frame(
        run_accession = exprs %>% colnames(),
        cluster_label = .)

ca_runs_final <- ca_runs_final %>%
    dplyr::left_join(
        clusters,
        by = c("run_accession"))

ca_runs_final %>%
    dplyr::count(study_accession, cluster_label) %>%
    tidyr::pivot_wider(id_cols = "study_accession", names_from = "cluster_label", values_from = "n")


###############
# embed genes #
###############
# translating to monocle terminology
# "gene" expression_data columns <- rna-seq run accession
# "cell" expression_data rows    <- gene (chromosome feature)

exprs <- reshape2::acast(
    data = estimated_expression,
    formula = gene_id ~ run_accession,
    value.var = "FPKM_normed") %>%
    t()

gene_metadata <- data.frame(run_accession = exprs %>% rownames) %>%
    dplyr::left_join(
        ca_runs_final %>%
        dplyr::mutate(
            gene_short_name = run_accession),
        by = "run_accession") %>%
    dplyr::mutate(
        study_accession = as.factor(study_accession))
row.names(gene_metadata) <- gene_metadata$gene_short_name

cell_metadata <- data.frame(feature_name = exprs %>% colnames())
row.names(exprs) <- gene_metadata$gene_short_name
row.names(cell_metadata) <- cell_metadata$feature_name

coexp_cds <- monocle3::new_cell_data_set(
    expression_data = exprs,
    gene_metadata = gene_metadata,
    cell_metadata = cell_metadata)

coexp_cds <- coexp_cds %>% monocle3::preprocess_cds(
    num_dims = 500)

coexp_cds <- coexp_cds %>%
    monocle3::reduce_dimension(
        preprocess_method = "PCA",
        #umap.min_dis = 0.01,
        #spread = .3,
        a = 50,
        b = 0.5,
        umap.n_neighbors = 30L,
        verbose = TRUE,
        n_epochs = 2000,
        negative_sample_rate = 50,
        repulsion_strength = 3)

# march 2019
coexp_cds <- coexp_cds %>%
    monocle3::reduce_dimension(
        preprocess_method = "PCA",
        #umap.min_dis = 0.01,
        #spread = .3,
        a = 20,
        b = 0.5,
        umap.n_neighbors = 30L,
        verbose = TRUE,
        n_epochs = 2000,
        negative_sample_rate = 50,
        repulsion_strength = 3)


coexp_cds <- coexp_cds %>%
    monocle3::cluster_cells(
        k = 30,
        num_iter = 10,
        resolution = .00001,
        verbose = TRUE)

save(coexp_cds, file = paste0("intermediate_data/coexp_cds_by_gene_", date_code, ".Rdata"))


coexp_cds %>%
    monocle3::plot_cells(
        show_trajectory_graph = FALSE,
        cell_size = 0.8)

ggplot2::ggsave(
    filename = paste0("product/figures/CalCEN_embedding/UMAP_genes_", date_code, ".pdf"),
    width = 6,
    height = 6)


# write out gene embedding and clusters
gene_clusters <- coexp_cds@clusters[["UMAP"]]$clusters %>%
    data.frame(
        feature_name = exprs %>% colnames(),
        cluster_label = .)

load("intermediate_data/chromosome_features.Rdata")
gene_clusters <- gene_clusters %>%
    dplyr::left_join(
        chromosome_features,
        by = "feature_name")

umap_coordinates <- coexp_cds %>% reducedDim("UMAP") %>%
    as.data.frame() %>%
    dplyr::rename(
        UMAP_1 = 1,
        UMAP_2 = 2)

gene_clusters <- gene_clusters %>%
    dplyr::bind_cols(umap_coordinates)

gene_clusters %>% readr::write_tsv(
    paste0("product/figures/CalCEN_embedding/UMAP_genes_cluster_labels_", date_code, ".tsv"))




### plot with old clusters ###
old_gene_clusters <- readr::read_tsv(
    "product/figures/CalCEN_embedding/UMAP_genes_cluster_labels_20201113.tsv")

data.frame(
    old_clustering = as.character(old_gene_clusters$cluster_label),
    new_clustering = as.character(gene_clusters$cluster_label)) %>%
    dplyr::count(old_clustering, new_clustering) %>%
    tidyr::pivot_wider(
        names_from = new_clustering,
        values_from = n)


#                      new_clustering
#     old_clustering   `1`  `10`  `12`  `13`  `14`  `15`   `2`   `3`   `4`   `5`   `6`   `7`   `8`   `9`  `18`  `11`  `20`  `19`  `21`  `17`  `16`
#   1 1                654     3    77     1     2     2     1     1     3     2     3    18     1     3    NA    NA    NA    NA    NA    NA    NA
#   3 11                 1   138     1     2     2    NA    31    15    NA    NA    NA     2    43    NA     3    NA    NA     1    29    NA    NA
#  13 4                  6   128   195    NA   240     1     1    NA    15    NA    12     1     4     3    NA    NA    NA    NA    NA    NA    NA
#  17 8                  4     1    NA   256    NA    NA    NA     1    NA     2    NA    11    NA    NA     1     1    NA    74    NA    NA    NA
#  11 2                  1     2    NA     1    NA    NA   532     9    NA    NA    NA    NA     4    NA    NA    NA    69     1     1    NA    NA
#  12 3                 NA     7    NA    NA     1   195    NA     1   299     4     2     2     1    NA    NA     2    NA     1    NA    97    NA
#   2 10                NA    NA    NA    NA    NA    NA     1   289    NA    NA    NA    NA    NA    NA     1    NA     3    NA    NA    NA    NA
#   5 13                NA    NA    NA    NA    NA    NA     2   175    NA     1    NA    NA     1    NA    NA    NA    NA    NA    NA    NA    NA
#  14 5                 27     2     6    NA    11     5    NA    NA   183     1   339     3    NA     4     1     1    NA    NA    NA     1    NA
#  15 6                  1    NA    NA     1    NA    NA     1    NA     1   425     1    NA    NA    NA    NA     3    NA    NA    NA    NA    NA
#  16 7                  2     4     1     3     3     7    NA    NA     2     5     7   302     2     5    NA    54    NA     1    NA    NA    NA
#   6 14                 1     5    NA    NA    NA    NA     5    NA    NA    NA    NA    NA   135    NA     1    NA    NA    NA     4    NA    NA
#   7 15                NA    NA    NA    NA     3    NA    NA     1    NA    NA    NA    NA   131    NA    NA    NA    NA    NA    NA    NA    NA
#  18 9                 28    NA     1    NA    NA    NA    NA    NA     2     1     4     1    NA   297    NA    NA    NA    NA    NA    NA    NA
#   9 17                NA     3     3    NA    NA    NA    NA     5    NA    NA    NA    NA     2    NA    80    NA    NA    NA     1    NA     5
#   4 12                NA     2    NA     4    NA     3    NA    NA    NA     8     4    12    NA    NA    NA   224    NA    NA    NA     1    NA
#  10 18                NA    NA    NA    NA    NA    NA    NA    31    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA
#   8 16                 1    NA    NA     2    NA    NA    NA     2    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    97


coexp_cds@clusters[["UMAP"]]$clusters <- old_gene_clusters$cluster_label %>% factor()

coexp_cds %>%
    monocle3::plot_cells(
        show_trajectory_graph = FALSE,
        cell_size = 0.8)

ggplot2::ggsave(
    filename = paste0("product/figures/CalCEN_embedding/UMAP_genes_old_clustering", date_code, ".pdf"),
    width = 6,
    height = 6)



##### Consensus Clustering #####

n_clusterings <- 2

clusterings <- data.frame(clustering = 1:n_clusterings) %>%
    dplyr::rowwise() %>%
    dplyr::do({
        data <- .
        cat("Computing embedding ", data$clustering, "\n", sep = "")
        coexp_cds <- coexp_cds %>%
            monocle3::reduce_dimension(
                preprocess_method = "PCA",
                init = "random",
                #umap.min_dis = 0.01,
                #spread = .3,
                a = 50,
                b = 0.5,
                umap.n_neighbors = 30L,
                verbose = TRUE,
                n_epochs = 200, #2000,
                negative_sample_rate = 50,
                repulsion_strength = 3)
        coexp_cds <- coexp_cds %>%
            monocle3::cluster_cells(
                k = 30,
                num_iter = 10,
                resolution = .001,
                verbose = TRUE)
        gene_clusters <- coexp_cds@clusters[["UMAP"]]$clusters %>%
            data.frame(
                clustering = data$clustering,
                feature_name = exprs %>% colnames(),
                cluster_label = .)
        gene_clusters
    })

clustering <- clusterings %>%
    dplyr::ungroup() %>%
    dplyr::mutate(cluster_label = as.character(cluster_label)) %>%
    tidyr::pivot_wider(
        id_cols = feature_name,
        names_from = clustering,
        values_from = cluster_label)
# these are all the same


library(Dune)
merger <- Dune::Dune(
    clusMat = clusterings,
    verbose = TRUE)
