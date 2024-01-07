
purity1 <- function(x) {

	stopifnot(x$type == "counts")

	p <- as.matrix(x$snp[,-(1:2)])
	i <- colSums(is.na(p)) < (0.2 * nrow(p))
	p <- p[,i]
	# not hetero
	p[p==0] <- NA

	i <- seq(1, nrow(p), 2)
	p1 <- p[i,]
	p2 <- p[i+1,]
	nv <- ncol(p1) 
	out <- rep(NA, nv)
	for (i in 1:nv) {
		p_i <- na.omit(cbind(p1[,i], p2[,i]))
		if (!(isTRUE(nrow(p_i) > 25))) next
		p_i <- p_i[rowSums(p_i)> 4, ]
		if (!(isTRUE(nrow(p_i) > 25))) next
		p_i <- apply(p_i, 1, sort) |> t()
		s <- p_i[,1] / p_i[,2]
		# if s if very small it cannot be meaningful
		out[i] <- median(s[s > 0.2])
	}
	out
}



purity <- function(xc, mr=0.2, minhet=0.2) {

	stopifnot(xc$type == "counts")
	p <- as.matrix(xc$snp[,-(1:2)])
	# not hetero
	p[p==0] <- NA

	i <- seq(1, nrow(p), 2)
	p1 <- p[i,]
	p2 <- p[i+1,]
	
	s <- p1 / p2
	j <- which(p1 > p2)
	s[j] <- p2[j] / p1[j]
	
	s[s < minhet] <- NA
	s[, colSums(is.na(s)) < ((1-mr) * nrow(s))] <- NA

	apply(s, 2, median, na.rm=TRUE)
}

