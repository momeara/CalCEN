# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

# https://github.com/mbojan/alluvial/blob/master/R/alluvial.R

alluvial <- function(
	...,
	freq,
	col="gray", border=0, layer, hide=FALSE, alpha=0.5,
	gap.width=0.05, xw=0.1, cw=0.1,
	blocks = TRUE,
	ordering=NULL,
	axis_labels=NULL,
	mar = c(2, 1, 1, 1),
	cex=par("cex"),
	xlim_offset= c(0, 0),
	ylim_offset= c(0, 0),
	cex.axis=par("cex.axis"),
	axes=TRUE,
	ann=TRUE,
	title = NULL){

	# Data and graphical parameters
	p <- data.frame( ..., freq=freq, col, alpha, border, hide, stringsAsFactors=FALSE)
	np <- ncol(p) - 5                    # Number of dimensions
	# check if 'ordering' is of proper form
	if( !is.null(ordering) ){
		stopifnot(is.list(ordering))
		if( length(ordering) != np ){
			stop("'ordering' argument should have ", np, " components, has ", length(ordering))
    }
	}
	n <- nrow(p)
	# Layers determine plotting order
	if(missing(layer)){
		layer <- 1:n
	}
	p$layer <- layer
	d <- p[ , 1:np, drop=FALSE]          # Dimensions dframe
	p <- p[ , -c(1:np), drop=FALSE]      # Parameteres dframe
	p$freq <- with(p, freq/sum(freq))    # Frequencies (weights)
	# Converting colors to hexcodes
	col <- col2rgb(p$col, alpha=TRUE)
	if(!identical(alpha, FALSE)) {
		col["alpha", ] <- p$alpha*256
	}
	p$col <- apply(col, 2, function(x) do.call(rgb, c(as.list(x), maxColorValue = 256)))
	# convert character vectors in data to factors
	isch <- sapply(d, is.character)
	d[isch] <- lapply(d[isch], as.factor)
	# Convert blocks to vector
	if (length(blocks) == 1){
		blocks <- if (!is.na(as.logical(blocks))){
			rep(blocks, np)
		} else if (blocks == "bookends") {
			c(TRUE, rep(FALSE, np - 2), TRUE)
		}
	}
	# Axis labels
	if(is.null(axis_labels)) {
		axis_labels <- names(d)
	} else {
		if(length(axis_labels) != ncol(d)){
			stop("`axis_labels` should have length ", names(d), ", has ", length(axis_labels))
		}
	}
	# Compute endpoints of flows (polygons)
	# i = dimension id
	# d = data frame of dimensions
	# f = weights
	# w = gap between categories
	cat("data matrix:\n")
	print(d)
	cat("\n\n")
	getp <- function(i, d, f, w=gap.width) {
		# Ordering dimension ids for lexicographic sorting
		a <- c(i, (1:ncol(d))[-i])
		# Order of rows of d starting from i-th dimension
		if( is.null(ordering[[i]]) ){
			o <- do.call(order, d[a])
		} else {
			d2 <- d[a]
			d2[1] <- ordering[[i]]
			o <- do.call(order, d2)
		}
		cat("dimension lexicographic sorting: ", paste0(a, collapse=" "), "\n", sep="")
		cat("requested order for dimension ", i, ": ", paste0(ordering[[i]], collapse=" "), "\n", sep="")
		cat("Sort order for column ", i, ":\n", sep="")
		print(d[o,])
		cat("\n\n")
		# Breakpoints on a dimension
		x <- c(0, cumsum(f[o])) * (1-w)
		# Stripe coordinates on a dimension
		x <- cbind(x[-length(x)], x[-1])
		# By how much stripes need to be shifted upwards (gap/max(gap))
		gap <- cumsum( c(0L, diff(as.numeric(d[o,i])) != 0) )
		mx <- max(gap)
		if (mx == 0) mx <- 1
		# shifts
		gap <- gap / mx * w
		# add gap-related shifts to stripe coordinates on dimension i
		(x + gap)[order(o),]
	}
	# Calculate stripe locations on dimensions: list of data frames. A component
	# for a dimension. Data frame contains 'y' locations of stripes.
	dd <- lapply(seq_along(d), getp, d=d, f=p$freq)

	# Plotting
	op <- par(mar=mar)
	plot(NULL, type="n", xlim=c(1-cw, np+cw) + xlim_offset, ylim=c(0, 1) + ylim_offset, xaxt="n", yaxt="n",
			 xaxs="i", yaxs="i", xlab='', ylab='', main=title, frame=FALSE)
	# For every stripe
	ind <- which(!p$hide)[rev(order(p[!p$hide, ]$layer))]
	for(i in ind ){
		# For every inter-dimensional segment
		for(j in 1:(np-1) ){
			# Draw stripe
			xspline( c(j, j, j+xw, j+1-xw, j+1, j+1, j+1-xw, j+xw, j) + rep(c(cw, -cw, cw), c(3, 4, 2)),
							c( dd[[j]][i, c(1, 2, 2)], rev(dd[[j+1]][i, c(1, 1, 2, 2)]), dd[[j]][i,c(1, 1)]),
							shape = c(0,0,1,1,0,0,1,1,0, 0),
							open=FALSE,
							col=p$col[i], border=p$border[i])
		}
	}

	# Category blocks with labels
	for(j in seq_along(dd)){
		ax <- lapply(split(dd[[j]], d[,j]), range)
		if (blocks[j]){
			for(k in seq_along(ax)){
				rect( j-cw, ax[[k]][1], j+cw, ax[[k]][2] )
			}
		} else {
			for (i in ind) {
				x <- j + c(-1, 1) * cw
				y <- t(dd[[j]][c(i, i), ])
				w <- xw * (x[2] - x[1])
				xspline(x = c(x[1], x[1], x[1] + w, x[2] - w,
										x[2], x[2], x[2] - w, x[1] + w, x[1]),
								y = c(y[c(1, 2, 2), 1], y[c(2, 2, 1, 1), 2], y[c(1, 1), 1]),
								shape = c(0, 0, 1, 1, 0, 0, 1, 1, 0, 0),
								open = FALSE, col = p$col[i], border = p$border[i])
			}
		}
		for(k in seq_along(ax)) {
			if(ann) text( j, mean(ax[[k]]), labels=names(ax)[k], cex=cex)
		}
	}
	# X axis
	if(axes) {
		axis(1, at= rep(c(-cw, cw), ncol(d)) + rep(seq_along(d), each=2),
				 line=0.5, col="white", col.ticks="black", labels=FALSE)
		axis(1, at=seq_along(d), tick=FALSE, labels=axis_labels, cex.axis=cex.axis)
	}
	par(op)
	rval <- list(
		# Endpoints of alluvia
		endpoints = do.call(
			"rbind",
			lapply(
				seq(along=dd),
				function(i) {
					df <- as.data.frame(
						structure(dd[[i]], dimnames=list(NULL,c(".bottom", ".top"))),
						stringsAsFactors=FALSE)
					df$.axis <- i
					cbind(d, df)
				})))
		# Category midpoints
		rval$category_midpoints <- structure(
			lapply(
				seq(1, ncol(d)),
				function(i) {
					mi <- with(
						rval$endpoints[rval$endpoints$.axis == i , ],
						tapply(.bottom, d[[i]], min))
					ma <- with(
						rval$endpoints[ rval$endpoints$.axis == i , ],
						tapply(.top, d[[i]], max))
					(mi + ma)/2
				}),
			names = names(d))
	# alluvium midpoints
	rval$alluvium_midpoints <- rval$endpoints %>%
		tidyr::gather_(".endpoint", ".value", c(".bottom", ".top")) %>%
		dplyr::group_by_( .dots=c(names(d), ".axis")) %>%
		dplyr::summarise_(.dots = list(m = ~mean(.value))) %>%
		dplyr::arrange_(.dots=c(names(d), ".axis")) %>%
		dplyr::group_by_( .dots = names(d) ) %>%
		dplyr::mutate_(.dots=list(
			.axis_from = ~dplyr::lag(.axis),
			.axis_to = ".axis",
			.x = ~(.axis_from + .axis_to)/2,
			.y = ~(m + dplyr::lag(m))/2,
			.slope = ~ (m - dplyr::lag(m)) / (.axis_to - .axis_from - cw))) %>%
		dplyr::ungroup() %>%
		dplyr::filter_(~!is.na(.axis_from)) %>%
		dplyr::select_(.dots=c(names(d), ".axis_from", ".axis_to", ".x", ".y", ".slope")) %>%
		as.data.frame(stringsAsFactors=FALSE)

	invisible(rval)
}
