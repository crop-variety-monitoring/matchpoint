
var_dist <- function(d, varieties) {
	d <- as.matrix(d)
	diag(d) <- NA
	u <- unique(varieties)
	v <- sapply(1:length(u), \(i) {
		j <- which(u[i] == varieties)
		dv <- d[j, j]
		mean(dv, na.rm=TRUE)
	})
	v[!is.finite(v)] <- 0
	v
}

break_far_relatives <- function(dstm, vars, threshold, verbose=FALSE) {
	adj <- (dstm < threshold)
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

		# do not break nearest neighbors
		for (j in 1:n) {
			ids <- i[m == j]
			if (any(is.na(ids))) next
			d <- dstm[ids, ,drop=FALSE]
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
		if (length(i) > 0) {
			vars[i] <- paste0(vars[i], "__", letters[m])
			if (verbose) print(vars[i])
		}
	}
	vars
}



clean_reference <- function(dst, vars, threshold, nmax=100, verbose=FALSE) {
	sampids <- colnames(dst)
	colnames(dst) <- rownames(dst) <- vars
	dst <- (dst < threshold)
	diag(dst) <- TRUE
	n <- 1
	while(TRUE) {
		x <- apply(dst, 1, which)
		y <- sapply(x, \(i) { length(unique(names(i))) > 1 })
		idx <- which(y)
		if (length(idx) == 0) break
		ids <- x[[ idx[1] ]]
		nms <- sort(unique(unlist(strsplit(names(ids), "-#-"))))
		nms <- paste0(sort(nms), collapse="-#-")
		colnames(dst)[ids] <- rownames(dst)[ids] <- nms 
		n = n + 1
		if (n > nmax) {
			warning("clean ref unsuccessful")
			break
		}
	}
	if (verbose) print(n)

	d <- data.frame(id=1:length(sampids), sample=sampids, from=vars, comb=colnames(dst), to=colnames(dst))

	d$changed <- d$from != d$to
	dc <- d[d$changed, ]
	if (nrow(dc) > 0) {
		ag <- aggregate(dc$from, dc["to"], \(i) {
				tab <- table(i) / length(i)
				j <- tab > 0.5
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
	d$id <- NULL
	rownames(d) <- NULL
	d
}

