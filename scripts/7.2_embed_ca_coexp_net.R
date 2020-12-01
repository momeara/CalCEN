

library(plyr)
library(tidyverse, quietly = TRUE, warn.conflicts = FALSE)
library(monocle3, quietly = TRUE, warn.conflicts = FALSE)


load("intermediate_data/ca_runs_final.Rdata")
load("intermediate_data/estimated_expression.Rdata")
load("product/CalCEN_network_full_20201007.Rdata")
date_code <- "20201113"


estimated_expression <- estimated_expression %>%
		dplyr::semi_join(ca_runs_final, by = c("run_accession"))

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
    data=estimated_expression,
    formula=gene_id ~ run_accession,
    value.var="FPKM") %>%
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


coexp_cds <- coexp_cds %>%
    monocle3::cluster_cells(
        k = 30,
        num_iter = 10,
        resolution = .001,
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
    filename = "product/figures/CalCEN_embedding/UMAP_genes_20201113.pdf",
    width = 6,
    height = 6)


gene_clusters <- coexp_cds@clusters[["UMAP"]]$clusters %>%
    data.frame(
        feature_name = exprs %>% colnames(),
        cluster_label = .)

load("intermediate_data/chromosome_features.Rdata")
gene_clusters <- gene_clusters %>%
    dplyr::left_join(
        chromosome_features,
        by = "feature_name")

gene_clusters %>% readr::write_tsv("product/figures/CalCEN_embedding/UMAP_genes_cluster_labels_20201113.tsv")

