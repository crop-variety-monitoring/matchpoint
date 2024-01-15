

group_dend <- function(x, add=NULL, col=viridis::turbo, cex=0.5, ...) {
	x <- as.matrix(x)
	xlabs <- colnames(x)
	colnames(x) <- rownames(x) <- 1:ncol(x)
	hc <- stats::hclust(stats::as.dist(x))
	j <- as.integer(hc$order)
	lbs <- xlabs[j]
	id <- as.integer(as.factor(lbs))
	n <- length(unique(id))
	if (inherits(col, "function")) {
		cols <- col(n)
	} else {
		cols <- rep_len(col, n)
	}
	cols <- cols[id]
	if (!is.null(add)) {
		lbs <- paste0(lbs, add[j])	
	}
	hc <- stats::as.dendrogram(hc)
	hc <- dendextend::set(hc, "labels", lbs)
	hc <- dendextend::set(hc, "labels_cex", cex)
	hc <- dendextend::set(hc, "labels_col", cols)
	hc
}



