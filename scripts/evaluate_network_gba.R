# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(EGAD)
library(purrr)
library(glue)


evaluate_network_gba <- function(
	genes,
	annotation_sets,
	networks,
	nfold=3,
	degrees=1:length(networks)){

	if(is.data.frame(networks) | is.matrix(networks)){
		networks=list("network"=networks)
		degrees = c(1)
	} else if(is.list(networks)){
		if(!all(degrees <= length(networks) & degrees >= 0)){
			stop("Degrees must be a list of numbers in [1, ", length(networks), "] but instead it is ", paste0(degrees, collapse=", "), "\n", sep="")
		}
	}

	for(i in 1:length(networks)){
		if(!all(genes == rownames(networks[i]))){
			stop("The row names for the ", i, "th network, ", names(networks)[i], ", does not match the names in 'genes'\n", sep="")
		}
		if(!all(genes == colnames(networks[i]))){
			stop("The column names for the ", i, "th network, ", names(networks)[i], ", does not match the names in 'genes'\n", sep="")
		}
	}

	if(is.data.frame(annotation_sets) | is.matrix(annotation_sets)){
		annotation_sets=list("annotation_set"=annotation_sets)
	}

	for(i in 1:length(annotation_sets)){
		if(!all(genes == rownames(annotation_sets[i]))){
			stop("The row names for the ", i, "th annotation set, ", names(annotation_sets)[i], ", does not match the names in 'genes'\n", sep="")
		}
	}


	cat("Computing Guilt-by-Association gene function prediction for:\n")
	cat("  annotation sets:", annotation_sets %>% names %>% paste(collapse=", "), "\n")
	cat("  networks:", networks %>% names %>% paste(collapse=", "), "\n")
	cat("  combinations of k networks:", paste0(degrees, collapse=", "), "\n")

	networks %>% purrr::map2(names(.), ., function(network_id, network){
		if(!all(rownames(network) == genes)){
			stop("The row names of '{network_id}' network do not match the supplied genes\n" %>% glue::glue())
		}
		if(!all(colnames(network) == genes)){
			stop("The column names of {network_id} do not match the supplied genes\n" %>% glue::glue())
		}
	})

	z <- plyr::ldply(names(annotation_sets), function(anno_id){
		annotation_set <- annotation_sets[[anno_id]]

		plyr::ldply(degrees, function(network_degree){
		plyr::ldply(networks %>% names %>% combn(network_degree, simplify=FALSE),
			function(network_ids){
				network_id <- paste0(network_ids, collapse="|")

				cat("Computing GBA for '", anno_id, "' annotations, using the ", network_id, " network: ", sep="")

				net <- networks[network_ids] %>% purrr::reduce(`+`)

				rownames(net) <- genes
				colnames(net) <- genes
				z <- tryCatch({
					gba <- EGAD::run_GBA(net, annotation_set, nfold=nfold)
					cat(
						"mean: ", gba[[1]][,1] %>% mean(na.rm=TRUE),
						"std: ", gba[[1]][,1] %>% sd(na.rm=TRUE), "\n", sep=" ")
					tibble::tibble(
						anno_id=anno_id,
						network_id=network_id,
						auroc_mean = gba[[1]][,1] %>% mean(na.rm=TRUE),
						auroc_std = gba[[1]][,1] %>% sd(na.rm=TRUE))
				}, error=function(e){
					cat("ERROR: ", e$message, "\n", sep = "")
					tibble::tibble()
				})
		})
		})
	})
}
