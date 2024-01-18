
cut.dm <- function(x, d) {
	keep <- x
	x <- as.matrix(x)
	cns <- colnames(x)
	colnames(x) <- 1:ncol(x)
	a <- (x <= d)
	g <- igraph::graph_from_adjacency_matrix(a)
	n <- igraph::count_components(g)
	if (n == 1) return(x)
	gg <- igraph::decompose(g)
	z <- lapply(1:length(gg), \(i) 
		cbind(group=i, cases=igraph::clusters(gg[[i]])$membership |> names() |> as.integer()))
	z <- do.call(rbind, z) |> data.frame()
	z <- z[order(z$cases), ]
	a <- aggregate(x, z["group"], mean, na.rm=TRUE)
	a <- aggregate(t(a[,-1]), z["group"], mean, na.rm=TRUE)
	nms <- a[,1]
	a <- as.matrix(a[,-1])
	rownames(a) <- colnames(a) <- nms
	z$labels <- cns
	list(dm=a, key=z)
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


