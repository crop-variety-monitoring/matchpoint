
punity <- function(x, thresholds) {
	x <- as.matrix(x)
	stopifnot(nrow(x) == ncol(x))
	stopifnot(length(unique(colnames(x))) < ncol(x))
	stopifnot(all(thresholds >= 0))
	pure <- rep(NA, length(thresholds))
	unis <- rep(NA, length(thresholds))
	for (i in 1:length(thresholds)) {
		a <- (x < thresholds[i])
		g <- igraph::graph_from_adjacency_matrix(a, "undirected", diag=FALSE)	
		n <- igraph::count_components(g)
		if (n == 1) {
			pure[i] <- 0
			unis[i] <- 1
			next
		}
		gg <- igraph::decompose(g)
		z <- lapply(1:length(gg), \(i) 
			cbind(id=i, names=igraph::components(gg[[i]])$membership |> names() |> unique()))
		z <- do.call(rbind, z)
		pure[i] <- length(unique(z[, "id"])) / nrow(z)
		tab <- table(z[, "names"])
		unis[i] <- length(tab[tab==1]) / length(tab)
	}
	cbind(threshold=thresholds, purity=pure, unison=unis, mean=(pure+unis)/2)
}



split_groups <- function(x, threshold, keep_nngb=TRUE, verbose=FALSE) {
	keep <- x
	x <- as.matrix(x)
	stopifnot(nrow(x) == ncol(x))
#	stopifnot(all(colnames(x) == rownames(x)))
	stopifnot(length(unique(colnames(x))) < ncol(x))
	stopifnot(threshold > 0)
	if (threshold >= max(x, na.rm=TRUE)) return(x)

	adj <- (x < threshold)
	vars <- colnames(x)
	uvars <- unique(vars)
	
	for (v in uvars) {
		i <- which(vars == v)
		a <- adj[i,i]
		if (all(a)) next
		
		# break far apart clusters 
		g <- igraph::graph_from_adjacency_matrix(a, "undirected", diag=FALSE)
		
		n <- igraph::count_components(g)
		if (n == 1) next 
		m <- igraph::components(g)$membership

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
			vars[i] <- paste0(vars[i], "_", letters[m])
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
	if (length(unique(colnames(x))) == ncol(x)) {
		stop("no duplicates to work with")
	}
	oldnms <- nms <- colnames(x)
	dimnames(x) <- list(1:ncol(x), 1:ncol(x))

	adj <- (x < threshold)
	g <- igraph::graph_from_adjacency_matrix(adj, "undirected", diag=FALSE)
	
	n <- igraph::count_components(g)
	if (n == ncol(adj)) {
		return(keep)
	}
	gg <- igraph::decompose(g)
	for (i in 1:n) {
		j <- igraph::components(gg[[i]])$membership |> names() |> as.integer()
		inms <- unique(nms[j])
		if (length(inms) > 1) {
			nms[j] <- paste0(sort(inms), collapse="_#_")
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
					names(j[j])
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
		d$to <- sapply(strsplit(d$to, "_#_"), \(n) {
					paste(sort(unique(n)), collapse="_#_")
				})
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


split_lump <- function(mdst, lump=0.05, split=0.2) {

	oldnms <- colnames(mdst)

	mdst <- matchpoint::split_groups(mdst, split, TRUE)
	mdst <- matchpoint::lump_similar(mdst, lump) 
	newnms <- colnames(mdst)

	i <- !grepl("_#_", newnms)
	j <- (newnms[i] != oldnms[i]) & (!grepl("_.$", newnms[i]))
	newnms[i][j] <- paste0(newnms[i][j], "*")

	i <- grep("_#_", newnms)
	newnms[i] <- gsub("_.$", "", newnms[i])
	newnms[i] <- gsub("_._#_", "_#_", newnms[i])

	simplify <- function(names, pattern) {
		nms <- names[!grepl("_#_", names)] |> unique()
		p <- unique(grep(pattern, nms, value=TRUE))
		if (length(p) > 0) {
			for (i in 1:length(p)) {
				n <- grep(gsub(pattern, "", p[i]), nms, value=TRUE) |> unique() |> length() 
				if (n == 1) {
					names[names==p[i]] <- gsub(pattern, "", p[i])
				}
			}
		}
		names
	}
	newnms <- simplify(newnms, "_a$")
	newnms <- simplify(newnms, "_b$")
	
	newnms <- collapse_names(newnms)
	newnms <- gsub("_#_", " # ", newnms)

	list(old=oldnms, new=newnms, split=split, lump=lump)
}

split_lump_old <- function(mdst, maxlump=0.05, minsplit=0.2) {

	oldnms <- colnames(mdst)
	v <- matchpoint::self_dist(mdst)
	qth <- stats::quantile(v$value)
	#outh <- max(minsplit, qth[4] + 1.5 * (qth[4] - qth[2]))
	outh <- minsplit
	splitln <- Inf
	if (outh < qth[5]) {
		splitln <- outh
		mdst <- matchpoint::split_groups(mdst, outh, TRUE)
		v <- matchpoint::self_dist(mdst)
		qth <- stats::quantile(v$value)
#		outh <- max(minsplit, qth[4] + 1.5 * (qth[4] - qth[2]))
#		if ((outh < qth[5]) && (outh < splitln) && (outh > minsplit)) {
#			splitln <- outh
#			mdst <- matchpoint::split_groups(mdst, outh, TRUE)
#			v <- matchpoint::self_dist(mdst)
#			qth <- quantile(v$value)
#		}
	}

	mdst <- matchpoint::lump_similar(mdst, max(.01, qth[2])) 

	v <- matchpoint::self_dist(mdst)
	cut <- min(maxlump, max(0.01, stats::quantile(v$value)[3]))
	mdst <- matchpoint::lump_similar(mdst, cut) 

	v <- matchpoint::self_dist(mdst)
	cut <- min(maxlump, max(0.01, stats::quantile(v$value)[4]))
	mdst <- matchpoint::lump_similar(mdst, cut) 
	
	newnms <- colnames(mdst)

	i <- !grepl("_#_", newnms)
	j <- (newnms[i] != oldnms[i]) & (!grepl("_.$", newnms[i]))
	newnms[i][j] <- paste0(newnms[i][j], "*")

	i <- grep("_#_", newnms)
	newnms[i] <- gsub("_.$", "", newnms[i])
	newnms[i] <- gsub("_._#_", "_#_", newnms[i])

	simplify <- function(names, pattern) {
		nms <- names[!grepl("_#_", names)] |> unique()
		p <- unique(grep(pattern, nms, value=TRUE))
		if (length(p) > 0) {
			for (i in 1:length(p)) {
				n <- grep(gsub(pattern, "", p[i]), nms, value=TRUE) |> unique() |> length() 
				if (n == 1) {
					names[names==p[i]] <- gsub(pattern, "", p[i])
				}
			}
		}
		names
	}
	newnms <- simplify(newnms, "_a$")
	newnms <- simplify(newnms, "_b$")
	
	newnms <- collapse_names(newnms)
	newnms <- gsub("_#_", " # ", newnms)

	list(old=oldnms, new=newnms, split=splitln, lump=cut)
}

