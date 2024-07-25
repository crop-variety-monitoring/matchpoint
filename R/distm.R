
cut.dm <- function(x, d) {
	keep <- x
	x <- as.matrix(x)
	cns <- colnames(x)
	colnames(x) <- 1:ncol(x)
	a <- (x <= d)
	g <- igraph::graph_from_adjacency_matrix(a)
	n <- igraph::count_components(g)
	if (n == 1) {
		return(list(dm=x, key=data.frame(groups=1:length(cns), cases=1, labels=cns)))
	}
	gg <- igraph::decompose(g)
	z <- lapply(1:length(gg), \(i) 
		cbind(group=i, cases=igraph::components(gg[[i]])$membership |> names() |> as.integer()))
	z <- do.call(rbind, z) |> data.frame()
	z <- z[order(z$cases), ]
	a <- aggregate(x, z["group"], mean, na.rm=TRUE)
	a <- aggregate(t(a[,-1]), z["group"], mean, na.rm=TRUE)
	nms <- a[,1]
	a <- as.matrix(a[,-1])
	rownames(a) <- colnames(a) <- nms
	rownames(z) <- NULL
	z$labels <- cns
	list(dm=a, key=z)
}


ngroups.dm <- function(x, d) {
	keep <- x
	x <- as.matrix(x)
	cns <- colnames(x)
	colnames(x) <- 1:ncol(x)
	a <- (x <= d)
	g <- igraph::graph_from_adjacency_matrix(a)
	igraph::count_components(g)
}


#ff <- function() {
#	r <- sapply(1:n, \(i) {
#		e <- sort(cns[z$cases[z$group==i]])
#		e <- sort(e[e != ""])
#		if (length(unique(e)) > 1) {
#			tab <- as.data.frame(table(e))
#			e <- apply(tab, 1, \(j) paste0(j[1], " (", j[2], ")"))
#			e <- paste0(e, collapse="; ")
#		} else {
#			e <- paste0(e, " (", length(e), ")")
#		}
#		e
#	})
#	r <- data.frame(group=1:n, cases=r)
#}

old_var_groups <- function(x, d, ref.id) {
	a <- (x >= d)
	g <- igraph::graph_from_adjacency_matrix(a)
	n <- igraph::count_components(g)
	gg <- igraph::decompose(g)
	z <- lapply(1:length(gg), \(i) 
		cbind(group=i, cases=igraph::components(gg[[i]])$membership |> names()))
	z <- do.call(rbind, z) |> data.frame()
	z$ref <- z$cases %in% ref.id
	gs <- aggregate(z["ref"], z["group"], \(i) c(length(i), sum(i)))
	gs <- cbind(gs[1], gs[[2]])
	colnames(gs)[2:3] <- c("size", "nrefs")
	z <- merge(z, gs, by="group")

	i <- z$size == z$nrefs

	zi <- z[!i, ]
	zr <- unique(zi[zi$ref, 1:2])
	zr <- aggregate(zr["cases"], zr["group"], \(i) paste(i, collapse=";"))
	names(zr)[2] <- "ref_sample"
	
	m <- merge(z[, c("group", "cases", "size")], zr, by="group", all.x=TRUE)
	m <- m[order(is.na(m$ref_sample), -m$size, m$group), ]
	m$group <- m$group |> factor(levels=unique(m$group)) |> as.integer()

	j <- is.na(m$ref_sample)
	m$cases[j] <- paste0("nr_", m$group[j] |> factor(levels=unique(m$group[j])) |> as.integer())
	colnames(m)[2] <- "field_sample"
	m
}


var_groups_x <- function(x, d, ref.id, varnames) {

	fap <- function(i, d, nms) {
		#i <- i[i >= (1-2*(1-d))]
		j <- which(i >= d)
		if (length(j) == 0) return("")
		if (length(j) > 1) {
			j <- j[i[j] >= (max(i) - .25*(1-d))]
		}
		paste0(unique(nms[j]), collapse=";")
	}


	varfact <- as.factor(varnames)
	varint <- as.integer(varfact)
#	varcode <- paste0("var_", varint)
	i <- colnames(x) %in% ref.id
	z <- x[i, !i]
	m <- match(rownames(z), ref.id)
	rownames(z) <- varint[m]

	s1 <- apply(z, 2, fap, d=d, nms=rownames(z))
	
	#s1b <- apply(z >= d, 2, \(i) unique(rownames(z)[i]))
	#unique(cbind(s1, s1b))

	s1 <- sapply(s1, \(i) paste(i, collapse=";"))
	n1 <- s1[s1 != ""]

	ref.id2 <- c(ref.id, names(n1))
	varint2 <- c(varint, n1)
	i <- colnames(x) %in% ref.id2
	z2 <- x[i, !i] # >= d
	m2 <- match(rownames(z2), ref.id2)
	rownames(z2) <- varint2[m2]

	s1 <- apply(z2, 2, fap, d=d, nms=rownames(z2))

	#s2 <- apply(z2, 2, \(i) unique(rownames(z2)[i]))
	s2 <- sapply(s2, \(i) paste(i, collapse=";"))
	n2 <- s2[s2 != ""]

	nn <- c(n1, n2)
	nn <- data.frame(ref=as.integer(names(nn)), id=nn)

	a <- (x >= d)
	g <- igraph::graph_from_adjacency_matrix(a)
	n <- igraph::count_components(g)
	gg <- igraph::decompose(g)
	z <- lapply(1:length(gg), \(i) 
		cbind(group=i, cases=igraph::components(gg[[i]])$membership |> names()))
	z <- do.call(rbind, z) |> data.frame()
	z$ref <- z$cases %in% ref.id
	gs <- aggregate(z["ref"], z["group"], \(i) c(length(i), sum(i)))
	gs <- cbind(gs[1], gs[[2]])
	colnames(gs)[2:3] <- c("size", "nrefs")
	z <- merge(z, gs, by="group")

	i <- z$size == z$nrefs

	zi <- z[!i, ]
	zr <- unique(zi[zi$ref, 1:2])
	zr <- aggregate(zr["cases"], zr["group"], \(i) paste(i, collapse=";"))
	names(zr)[2] <- "ref_sample"
	
	m <- merge(z[, c("group", "cases", "size")], zr, by="group", all.x=TRUE)
	m <- m[order(is.na(m$ref_sample), -m$size, m$group), ]
	m$group <- m$group |> factor(levels=unique(m$group)) |> as.integer()

	j <- is.na(m$ref_sample)
	m$cases[j] <- paste0("nr_", m$group[j] |> factor(levels=unique(m$group[j])) |> as.integer())
	colnames(m)[2] <- "field_sample"
	m
}



var_groups <- function(x, d, ref.id, varnames) {
	a <- (x >= d)
	g <- igraph::graph_from_adjacency_matrix(a, "undirected", diag=FALSE)
	n <- igraph::count_components(g)
	gg <- igraph::decompose(g)
	z <- lapply(1:length(gg), \(i) 
		cbind(group=i, cases=igraph::components(gg[[i]])$membership |> names()))
	z <- do.call(rbind, z) |> data.frame()
	z$ref <- z$cases %in% ref.id

	zz <- lapply(1:n, \(i) {
		k <- z[z$group==i, ]
		if (sum(k$ref) == 0) {
			k$varname <- ""
			return(k)
		}
		if (all(k$ref)) {
			k$varname <- "ref"
			return(k)
		}

		m1 <- match(k$cases, colnames(x))
		m2 <- match(k$cases[k$ref], colnames(x))
		s <- sapply(m2, \(i) x[cbind(m1, i)])
		s[s < (1 - 1.5 * (1-d))] <- NA
		
		k$varname <- apply(s, 1, \(i) paste0(unique(varnames[match(colnames(x)[m2[!is.na(i)]], ref.id)]), collapse=";"))
		k$group <- paste0(k$group, letters[as.integer(as.factor(k$varname))])
		k
	})
	z <- do.call(rbind, zz)
	z <- z[!z$ref, ]
	z$ref <- NULL
	
	gs <- aggregate(z["varname"], z["group"], \(i) c(length(i), length(unique(i[i!=""]))))
	gs <- cbind(gs[1], gs[[2]])
	colnames(gs)[2:3] <- c("size", "nvars")

	z <- merge(z, gs, by="group")
	zr <- unique(z[zi$ref, 1:2])
	zr <- aggregate(zr["cases"], zr["group"], \(i) paste(unique(i), collapse=";"))
	names(zr)[2] <- "ref_sample"
	
	
	m <- merge(z[, c("group", "cases", "size")], zr, by="group", all.x=TRUE)
	m <- m[order(is.na(m$ref_sample), -m$size, m$group), ]
	m$group <- m$group |> factor(levels=unique(m$group)) |> as.integer()

	j <- is.na(m$ref_sample)
	m$cases[j] <- paste0("nr_", m$group[j] |> factor(levels=unique(m$group[j])) |> as.integer())
	colnames(m)[2] <- "field_sample"
	m
}

