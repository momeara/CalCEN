
source("parameters.R")

if (!dir.exists(paths = "intermediate_data/crossmap")) {
    cat("Creating 'intermediate_data/crossmap' directory ...\n")
    dir.create(path = "intermediate_data/crossmap")
}

# convert .gff to .bed
# https://bedops.readthedocs.io/en/latest/content/reference/file-management/conversion/gff2bed.html
# BEDOPS v2.4.39
command <- paste0(
    "gff2bed < raw_data/C_albicans_SC5314_A22_current_features.gff > intermediate_data/crossmap/c_albicans_exon.bed")
cat(command, "\n", sep = "")
system(command)

command <- paste0(
    "gff2bed < raw_data/C_albicans_SC5314_A22_current_intergenic.gff > intermediate_data/crossmap/c_albicans_utr.bed")
cat(command, "\n", sep = "")
system(command)


# process annotation data
command <- paste0(
    "Rscript ",
    crossmap_path, "/gtf_to_txt.R ",
    "--gtf raw_data/C_albicans_SC5314_A22_current_features.gtf ",
    "-f exon,UTR ",
    "-o intermediate_data/crossmap/c_albicans_exon_utr.txt ")
cat(command, "\n", sep = "")
system(command)

# generate gene-mappability for k-mer mappability
command <- paste0(
    "Rscript ",
    crossmap_path, "/compute_mappability.R ",
    "--annot intermediate_data/crossmap/c_albicans_exon_utr.txt ",
    "--kmap_exon intermediate_data/crossmap/c_albicans_exon.bed ",
    "--kmap_utr intermediate_data/crossmap/c_albicans_utr.bed ",
    "--verbose 1 ",
    "-o intermediate_data/crossmap/gene_mappability.txt")
cat(command, "\n", sep = "")
system(command)
