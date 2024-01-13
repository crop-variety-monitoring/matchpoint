

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
	s[, colSums(is.na(s)) < (mr * nrow(s))] <- NA
	apply(s, 2, stats::median, na.rm=TRUE)
}

