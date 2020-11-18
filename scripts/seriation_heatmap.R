

library(seriation)
library(gplots)

#########
# given a matrix x,
# This makes a dendrogram where the columns are sorted to minimize the row-to-row difference
# the sorting is done relative to the ref_x, if it is supplied
seriation_heatmap <- function (
    x,
    ref_x = x,
    fname = NULL,
    height = 180,
    width = 15,
    color_scale = NULL,
    seriate_method = "OLO",
    seriate_control = NULL,
    verbose = FALSE){
    if(!is.null(ref_x)){
        if(verbose){
            cat("Using reference data.frame to compute row order ...\n")
        }

        if(any(class(ref_x) == "dist")){
            if(verbose){
                cat("ref_x is dist, so using it as is ...\n")
            }
            if(attr(ref_x, 'Size') != nrow(x)){
                cat("ERROR: the size of the provided ref_x is ", attr(ref_x, "Size"), " while nrow(x)=", nrow(x), "\n", sep = "")
            }
            dist_row <- ref_x
            dist_col <- t(ref_x)
        } else {
            if(verbose){
                cat("Computing euclidian distances between ref_x rows ...\n")
            }
            if(nrow(ref_x) != nrow(x)){
                cat("ERROR: nrow(ref_x) = ", nrow(ref_x), " !=  nrow(x) = ", nrow(x), "\n", sep = "")
            }
            dist_row <- dist(ref_x)
            dist_col <- dist(t(ref_x))
        }
        if(verbose){
            cat("Computing seriation with method '", seriate_method, "'\n", sep="")
        }
        o_row <- seriate(dist_row, method = seriate_method, control = seriate_control)[[1]]
        Rowv <- o_row %>% seriation::get_order()
        o_col <- seriate(dist_col, method = seriate_method, control = seriate_control)[[1]]
        Colv <- o_col %>% seriation::get_order()
    } else {
        Rowv <- FALSE
        Colv <- FALSE 
    }
    args <- list()
    if(is.null(color_scale)){
        if (any(x < 0, na.rm = TRUE)){
            if(verbose){
                cat("Scores have values that are less than zero, so using blue-red colors ...\n")
            }
            args$col <- seriation::bluered(n=100, bias=1)
        } else {
            if(verbose){
                cat("Scores are all greater than zero using white-black colors ...\n")
            }
            args$col <- seriation::greys(n=100, power=1)
        }
    } else {
        args$col <- color_scale
    }
    args$trace <- "none"
    args$density.info <- "none"
    args$cexRow <- 1
    args$cexCol <- 1
    args$dendrogram="none"
    args$key = FALSE
    args$keysize = 0.03
    args$colsep = seq(0, ncol(x), by=5000000)
    args$rowsep = seq(0, nrow(x), by=5000000)
    #args$sepwidth = c(0.02, 0.02)
    args$sepwidth = c(0, 0)
    args$margins = c(3,7)
    args$useRaster = TRUE
    args <- c(list(x = x, Colv = Colv, Rowv = Rowv), args)
    cat("plotting '", fname, "'...\n", sep="")
    pdf(fname, height=height, width=width)
        suppressWarnings(ret <- do.call(gplots::heatmap.2, args))
    dev.off()
}
