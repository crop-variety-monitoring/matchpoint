

min_dist <- function(x) {
	x <- as.matrix(x)
	stopifnot(all(nrow(x) == ncol(x)))
	stopifnot(length(unique(colnames(x))) < ncol(x))
	diag(x) <- NA
	i <- apply(x, 1, which.min)
	data.frame(group=colnames(x)[i], value=x[cbind(1:nrow(x), i)])
}


min_self_dist <- function(x) {
	x <- as.matrix(x)
	stopifnot(nrow(x) == ncol(x))
	stopifnot(length(unique(colnames(x))) < ncol(x))
	diag(x) <- NA
	vars <- colnames(x)
	u <- unique(vars)
	z <- rep(NA, ncol(x))
	for (i in 1:length(u)) {
		j <- which(u[i] == vars)
		if (length(j) > 1) {
			dv <- x[j, j, drop=FALSE]
			z[j] <- j[apply(dv, 1, which.min)]
		} 
	}
	v <- x[cbind(1:nrow(x), z)]
	data.frame(group=vars, value=v)
}


nngb_dist <- function(x, fun=min) {
	x <- as.matrix(x)
	stopifnot(all(nrow(x) == ncol(x)))
	#stopifnot(all(colnames(x) == rownames(x)))
	stopifnot(length(unique(colnames(x))) < ncol(x))
	vars <- colnames(x)
	u <- unique(vars)
	d <- sapply(1:length(u), \(i) {
		j <- which(u[i] == vars)
		dv <- x[j, ,drop=FALSE]
		dv[, j] <- NA
		fun(dv, na.rm=TRUE)
	})
	if (NCOL(d) == 1) {
		data.frame(group=u, value=d)
	} else { # e.g. fun=quantile
		data.frame(group=u, t(d))		
	}
}


self_dist <- function(x, fun=mean) {
	x <- as.matrix(x)
	stopifnot(nrow(x) == ncol(x))
	#stopifnot(all(colnames(x) == rownames(x)))
	stopifnot(length(unique(colnames(x))) < ncol(x))
	diag(x) <- NA
	vars <- colnames(x)
	u <- unique(vars)
	d <- sapply(1:length(u), \(i) {
		j <- which(u[i] == vars)
		dv <- x[j, j, drop=FALSE]
		fun(dv, na.rm=TRUE)
	})
	d[!is.finite(d)] <- 0
	if (NCOL(d) == 1) {
		data.frame(group=u, value=d)
	} else { # e.g. fun=quantile
		data.frame(group=u, t(d))		
	}
}


other_dist <- function(x, fun=mean) {
	x <- as.matrix(x)
	stopifnot(nrow(x) == ncol(x))
	stopifnot(length(unique(colnames(x))) < ncol(x))
	diag(x) <- NA
	vars <- colnames(x)
	x <- aggregate(x, list(vars), mean, na.rm=TRUE)[,-1]
	x <- aggregate(t(x), list(vars), mean, na.rm=TRUE)
	vars <- x[,1]
	x <- as.matrix(x[,-1])
	colnames(x) <- vars
	stopifnot(length(vars) == length(unique(vars)))

	d <- sapply(1:length(vars), \(i) {
		dv <- x[i, -i, drop=FALSE]
		fun(dv, na.rm=TRUE)
	})

	d[!is.finite(d)] <- 0
	if (NCOL(d) == 1) {
		data.frame(group=vars, value=d)
	} else { # e.g. fun=quantile
		data.frame(group=vars, t(d))		
	}
}

