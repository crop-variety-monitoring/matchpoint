
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

var_groups <- function(x, d, ref.id) {
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

	m <- merge(z[, c("group", "size")], zr, by="group", all.x=TRUE)
	m <- m[order(m$cases, -m$size, m$group), ]
	j <- is.na(m$cases)
	m$cases[j] <- paste0("nr_", m$group[j] |> factor(levels=unique(m$group[j])) |> as.integer())
	colnames(m)[3] <- "references"
	m
}

