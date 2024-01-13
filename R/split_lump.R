
split_groups <- function(x, threshold, keep_nngb=TRUE, verbose=FALSE) {
	keep <- x
	x <- as.matrix(x)
	stopifnot(nrow(x) == ncol(x))
#	stopifnot(all(colnames(x) == rownames(x)))
	stopifnot(length(unique(colnames(x))) < ncol(x))
	adj <- (x < threshold)
	vars <- colnames(x)
	uvars <- unique(vars)
	
	for (v in uvars) {
		i <- which(vars == v)
		a <- adj[i,i]
		if (all(a)) next
		
		# break far apart clusters 
		g <- igraph::graph_from_adjacency_matrix(a)
		n <- igraph::count_components(g)
		if (n == 1) next 
		m <- igraph::clusters(g)$membership

		if (keep_nngb) {
			for (j in 1:n) {
				ids <- i[m == j]
				if (any(is.na(ids))) next
				d <- x[ids, ,drop=FALSE]
				d[cbind(1:length(ids), ids)] <- Inf
				nearest <- apply(d, 1, \(k) which(k == min(k))) |> unlist() |> unique()
				nearest <- nearest[nearest %in% i[m != j]]
				if (length(nearest) > 0) {
					# this cluster 
					i[m == j] <- NA
					# near cluster(s)
					k <- which(i %in% nearest)
					i[m %in% unique(m[k])] <- NA
				}
			}
			m <- m[!is.na(i)]
			i <- i[!is.na(i)]
		}

		if (length(i) > 0) {
			vars[i] <- paste0(vars[i], "__", letters[m])
			if (verbose) print(table(vars[i]))
			if (inherits(keep, "matrix")) {
				dimnames(keep) <- list(vars, vars)
			} else {
				attr(keep, "Labels") <- vars 
			}
		}
	}
	keep
}



lump_similar <- function(x, threshold, verbose=FALSE) {

	keep <- x
	x <- as.matrix(x)
	stopifnot(nrow(x) == ncol(x))
#	stopifnot(all(colnames(x) == rownames(x)))
	stopifnot(length(unique(colnames(x))) < ncol(x))
	oldnms <- nms <- colnames(x)
	dimnames(x) <- list(1:ncol(x), 1:ncol(x))

	adj <- (x < threshold)
	g <- igraph::graph_from_adjacency_matrix(adj)
	n <- igraph::count_components(g)
	if (n == ncol(adj)) {
		return(keep)
	}
	gg <- igraph::decompose(g)
	for (i in 1:n) {
		j <- igraph::clusters(gg[[i]])$membership |> names() |> as.integer()
		inms <- nms[j]
		if (length(unique(inms)) > 1) {
			nms[j] <- paste0(sort(unique(inms)), collapse="-#-")
		}
	}

	d <- data.frame(id=1:length(nms), from=oldnms, comb=nms, to=nms)
	d$changed <- d$from != d$to
	dc <- d[d$changed, ]
	if (nrow(dc) > 0) {
		ag <- stats::aggregate(dc$from, dc["to"], \(i) {
				tab <- table(i) / length(i)
				if (length(i) < 4) {
					j <- FALSE
				} else if (length(tab) == 2) {
					j <- tab > 0.66
				} else {
					j <- tab > 0.55		
				}
				if (any(j)) {
					return(paste0(names(j[j]), " *"))
				} else {
					NA
				}})
		ag <- ag[!is.na(ag$x), ]
		if (nrow(ag) > 1) {
			nms <- colnames(d)
			d <- merge(d, ag, by="to", all.x=TRUE)
			d$to[!is.na(d$x)] <- d$x[!is.na(d$x)]
			d <- d[,nms]
		}
		d$to <- gsub("-#-", " # ", d$to)
#		i <- gsub(" \\*", "", d$to) == d$from
#		d$to[i] <- d$from[i]
		d$to <- gsub("\\* \\*$", "*", d$to)
		d <- d[order(d$id), ]
	}
	if (verbose) {
		print(unique(d[, c("from", "to")]))
	}
	if (inherits(keep, "matrix")) {
		dimnames(keep) <- list(d$to, d$to)
	} else {
		attr(keep, "Labels") <- d$to 
	}

	keep
}
	


