# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(stringr)
library(readr)


dir.create("raw_data/yeast_net", showWarnings = FALSE)

get_yeast_net <- function(tag=NULL, type, name=NULL){
	fname <- paste0("YeastNet.v3", ifelse(is.null(tag), "", paste0(".", tag)), ".txt")
	name_query <- ifelse(is.null(name), "", paste0("&name=", name))
	system(paste0("\\
	  cd raw_data/yeast_net && \\
	  wget -O ", fname, " 'https://www.inetbio.org/yeastnet/download.php?type=", type, name_query, "'"))

	readr::read_tsv(
		file=paste0("raw_data/yeast_net/", fname),
		col_names=c("gene1", "gene2", "score"),
		col_types=readr::cols(
			gene1=readr::col_character(),
			gene2=readr::col_character(),
			score=readr::col_double()))
}

yeast_net <- get_yeast_net("", 1)
save(yeast_net, file="intermediate_data/yeast_net.Rdata")

yeast_net_benchmark <- get_yeast_net("benchmark", 2) %>% dplyr::select(-score)
save(yeast_net_benchmark, file="intermediate_data/yeast_net_benchmark.Rdata")

yeast_net_CC <- get_yeast_net("CC", 3, "INT.CC.YeastNet.v3.4345gene.82319link.txt")
save(yeast_net_CC, file="intermediate_data/yeast_net_CC.Rdata")

yeast_net_CX <- get_yeast_net("CX", 3, "INT.CX.YeastNet.v3.5730gene.242504link.txt")
save(yeast_net_CX, file="intermediate_data/yeast_net_CX.Rdata")

yeast_net_DC <- get_yeast_net("DC", 3, "INT.DC.YeastNet.v3.3679gene.29880link.txt")
save(yeast_net_DC, file="intermediate_data/yeast_net_DC.Rdata")

yeast_net_GN <- get_yeast_net("GN", 3, "INT.GN.YeastNet.v3.1863gene.29475link.txt")
save(yeast_net_GN, file="intermediate_data/yeast_net_GN.Rdata")

yeast_net_GT <- get_yeast_net("GT", 3, "INT.GT.YeastNet.v3.4365gene.149498link.txt")
save(yeast_net_GT, file="intermediate_data/yeast_net_GT.Rdata")

yeast_net_HT <- get_yeast_net("HT", 3, "INT.HT.YeastNet.v3.5487gene.141347link.txt")
save(yeast_net_HT, file="intermediate_data/yeast_net_HT.Rdata")

yeast_net_LC <- get_yeast_net("LC", 3, "INT.LC.YeastNet.v3.5293gene.54421link.txt")
save(yeast_net_LC, file="intermediate_data/yeast_net_LC.Rdata")

yeast_net_PG <- get_yeast_net("PG", 3, "INT.PG.YeastNet.v3.2463gene.54496link.txt")
save(yeast_net_PG, file="intermediate_data/yeast_net_PG.Rdata")

yeast_net_TS <- get_yeast_net("TS", 3, "INT.TS.YeastNet.v3.1101gene.3510link.txt")
save(yeast_net_TS, file="intermediate_data/yeast_net_TS.Rdata")







