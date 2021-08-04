
#' parse BioGRID tab2 format
#'
#' Biogrid is a database curating protein-protein interactions
#' from hight- and low-throughput physical and genetic interactions
#'
#' @param fname path to .tab2 file
#' @param taxon filter for interactions where both partners are from ncbi <taxon>
#' @return tibble::tibble with columns
#'    biogrid_interaction_id
#'    gene_id_1
#'    gene_id_2
#'    biogrid_id_1
#'    biogrid_id_2
#'    feature_name_1
#'    feature_name_2
#'    gene_symbol_1
#'    gene_symbol_2
#'    synonyms_1
#'    synonyms_2
#'    experimental_system
#'    experimental_system_type
#'    author
#'    pubmed_id
#'    taxon_1
#'    taxon_2
#'    throughput
#'    score
#'    modification
#'    phenotypes
#'    qualifications
#'    tags
#'    source_database
#'
#' usage:
#' From the command line, download the by-organism biogrid dataset for the desired release and filter for the desired organism
#' 
#'     BIOGRID_RELEASE=3.4.161
#'     ORGANISM_NAME=Candida_albicans_SC5314
#'     pushd raw_data
#'     wget https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-${BIOGRID_RELEASE}/BIOGRID-ORGANISM-${BIOGRID_RELEASE}.tab2.zip
#'     unzip BIOGRID-${BIOGRID_RELEASE}/BIOGRID-ORGANISM-${BIOGRID_RELEASE}.tab2.zip
#'     ls | grep -v -e "${ORANISM_NAME} | xargs rm
#'     popd
#'
#' from R, parse biogrid data
#'
#'     biogrid_release <- "3.4.161"
#'     organism_name <- "Candida_albicans_SC5314"
#'     taxon <- "237561"
#'     biogrid_data <- CalCEN::read_biogrid_tab2(
#'        fname = paste0("raw_data/biogrid/BIOGRID-ORGANISM-", organism_name, "-", biogrid_release, ".tab2.txt"),
#'        taxon = taxon)
#' 
#' @export
read_biogrid_tab2 <- function(fname, taxon) {
  readr::read_tsv(
    file = fname,
    skip = 1,
    col_names = c(
      "biogrid_interaction_id",
      "gene_id_1",
      "gene_id_2",
      "biogrid_id_1",
      "biogrid_id_2",
      "feature_name_1",
      "feature_name_2",
      "gene_symbol_1",
      "gene_symbol_2",
      "synonyms_1",
      "synonyms_2",
      "experimental_system",
      "experimental_system_type",
      "author",
      "pubmed_id",
      "taxon_1",
      "taxon_2",
      "throughput",
      "score",
      "modification",
      "phenotypes",
      "qualifications",
      "tags",
      "source_database"),
    col_types=readr::cols(
      biogrid_interaction_id = readr::col_double(),
      gene_id_1 = readr::col_double(),
      gene_id_2 = readr::col_double(),
      biogrid_id_1 = readr::col_double(),
      biogrid_id_2 = readr::col_double(),
      feature_name_1 = readr::col_character(),
      feature_name_2 = readr::col_character(),
      gene_symbol_1 = readr::col_character(),
      gene_symbol_2 = readr::col_character(),
      synonyms_1 = readr::col_character(),
      synonyms_2 = readr::col_character(),
      experimental_system = readr::col_character(),
      experimental_system_type = readr::col_character(),
      author = readr::col_character(),
      pubmed_id = readr::col_double(),
      taxon_1 = readr::col_double(),
      taxon_2 = readr::col_double(),
      throughput = readr::col_character(),
      score = readr::col_character(),
      modification = readr::col_character(),
      phenotypes = readr::col_character(),
      qualifications = readr::col_character(),
      tags = readr::col_character(),
      source_database = readr::col_character())) %>%
    dplyr::filter(taxon_1 == taxon, taxon_2 == taxon) %>%
    dplyr::mutate(
      experimental_system_abbreviation = dplyr::case_when(
        experimental_system == "Dosage Lethality" ~ "DL",
        experimental_system == "Dosage Growth Defect" ~ "DGD",
        experimental_system == "Dosage Rescue" ~ "DR",
        experimental_system == "Synthetic Rescue" ~ "SR",
        experimental_system == "Phenotypic Suppression" ~ "PS",
        experimental_system == "Phenotypic Enhancement" ~ "PE",
        experimental_system == "Synthetic Lethality" ~ "SL",
        experimental_system == "Synthetic Growth Defect" ~ "SGD",
        experimental_system == "Positive Genetic" ~ "PG",
        experimental_system == "Negative Genetic" ~ "NG",
        experimental_system == "Synthetic Haploinsufficiency" ~ "HIP",
        experimental_system == "Affinity Capture-Luminescence" ~ "APL",
        experimental_system == "Proximity Label-MS" ~ "PL-MS",
        experimental_system == "Far Western" ~ "FW",
        experimental_system == "FRET" ~ "FRET",
        experimental_system == "Protein-RNA" ~ "P-RNA",
        experimental_system == "Co-localization" ~ "CoL",
        experimental_system == "Protein-peptide" ~ "Pro-Pep",
        experimental_system == "Co-crystal Structure" ~ "CoX",
        experimental_system == "Co-fractionation" ~ "CoF",
        experimental_system == "Co-purification" ~ "CoP",
        experimental_system == "Biochemical Activity" ~ "BA",
        experimental_system == "PCA" ~ "PCA",
        experimental_system == "Reconstituted Complex" ~ "ReC",
        experimental_system == "Two-hybrid" ~ "2H",
        experimental_system == "Affinity Capture-Western" ~ "AP-W",
        experimental_system == "Affinity Capture-RNA" ~ "AP-RNA",
        experimental_system == "Affinity Capture-MS" ~ "AP-MS",
        TRUE ~ experimental_system),
      throughput_abbreviation = ifelse(
          throughput == "High Throughput", "ht", "lt"))
}
