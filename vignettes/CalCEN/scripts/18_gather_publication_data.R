
library(plyr)
library(dplyr)
library(readr)

publication_data_path <- "product/publication_data_CalCEN_v1.0.0_20201230"

if (!dir.exists(paths = publication_data_path)) {
    cat("Creating pulication data path: ", publication_data_path, "\n", sep = "")
    dir.create(publication_data_path)
}

load("intermediate_data/ca_genes.Rdata")
ca_genes %>%
    data.frame(feature_name = .) %>%
    readr::write_tsv(paste0(publication_data_path, "/genes.tsv"))

load("intermediate_data/chromosome_features.Rdata")
chromosome_features %>%
    readr::write_tsv(paste0(publication_data_path, "/chromosome_features.tsv"))

load("intermediate_data/ca_runs_final.Rdata")
ca_runs_final %>%
    readr::write_tsv(paste0(publication_data_path, "/expression_runs.tsv"))

load("intermediate_data/estimated_expression.Rdata")
estimated_expression %>%
    readr::write_tsv(paste0(publication_data_path, "/estimated_expression.tsv"))


load("intermediate_data/CalCEN_full.Rdata")
CalCEN_full %>%
    as.data.frame() %>%
    readr::write_tsv(paste0(publication_data_path, "/CalCEN_network.tsv"))

load("intermediate_data/ca_blastp_rank_network_full.Rdata")
ca_blastp_rank_network_full %>%
    as.data.frame() %>%
    dplyr::filter(score != 0) %>%
    readr::write_tsv(paste0(publication_data_path, "/blastp_network.tsv"))

load("intermediate_data/ca_sac_ortholog_genetic_ppi_network.Rdata")
ca_sac_ortholog_genetic_ppi_network %>%
    data.frame() %>%
    tibble::rownames_to_column("feature_name_1") %>%
    tidyr::gather(key="feature_name_2", value="score", -feature_name_1) %>%
    readr::write_tsv(paste0(publication_data_path, "/sac_gene_network.tsv"))

load("intermediate_data/ca_sac_ortholog_physical_ppi_network.Rdata")
ca_sac_ortholog_physical_ppi_network %>%
    data.frame() %>%
    tibble::rownames_to_column("feature_name_1") %>%
    tidyr::gather(key="feature_name_2", value="score", -feature_name_1) %>%
    dplyr::filter(score != 0) %>%
    readr::write_tsv(paste0(publication_data_path, "/sac_phys_network.tsv"))

load("intermediate_data/yeast_net_network.Rdata")
z <- yeast_net_network %>%
    data.frame() %>%
    tibble::rownames_to_column("feature_name_1") %>%
    tidyr::gather(key="feature_name_2", value="score", -feature_name_1) %>%
    dplyr::filter(score != 0) %>%
    readr::write_tsv(paste0(publication_data_path, "/yeast_net_network.tsv"))

load("intermediate_data/ca_go_annotations.Rdata")
ca_go_annotations %>%
    data.frame() %>%
    tibble::rownames_to_column("feature_name") %>%
    tidyr::gather(key="go_id", value="annotation", -feature_name) %>%
    dplyr::filter(annotation != 0) %>%
    dplyr::select(-annotation) %>%
    readr::write_tsv(paste0(publication_data_path, "/go_annotations.tsv"))

file.copy(
    from = "product/gba_summary_10f_C-B-SP-SG-YN_20201124.tsv",
    to = paste0(publication_data_path, "/gba_summary.tsv"))


tar(
    tarfile = paste0(publication_data_path, ".tar.gz"),
    files = publication_data_path,
    compression = "gzip",
    compression_level = 9)
