
fake_it <- function(dart_file, outdir=".") {

#  dart_file <- "input/NGA/DCas23-7954_SNP.csv"
#  outdir <- "C:/github/cropvarmon/matchpoint/inst/ex"
	
	dart_file2 <- gsub(".csv$", "_2row.csv", dart_file)
	r1 <- data.frame(data.table::fread(dart_file, header=FALSE), check.names=FALSE)
	r2 <- data.frame(data.table::fread(dart_file2, header=FALSE), check.names=FALSE)

	srow <- sum(r1[1:30, 1] == "*") + 1
	scol <- sum(r1[1, 1:50] == "*") + 1

	take_anon_sample <- function(x, sid, nc=16) {
		y <- x[(srow+1):nrow(x), 1:nc]
		y <- y[y[,1] %in% sid, ]
		x <- x[1:(nrow(y)+srow), ]
		x[(srow+1):nrow(x), 1:nc] <- y
		x
	}
	
	s <- sort(sample(srow:nrow(r1), (nrow(r1)-srow-1)/2))
	sid <- r1[s,1]
	r1 <- take_anon_sample(r1, sid, scol-1)
	r2 <- take_anon_sample(r2, sid, scol)

	geno_file <- gsub("SNP.csv$", "variety-info.csv", dart_file)
	geno <- data.frame(data.table::fread(geno_file), check.names=FALSE)
	geno <- geno[sample(nrow(geno)), ]
	k <- match(geno$sample, r1[7,])
	geno <- geno[!is.na(k), ]

	i <- geno$variety == ""
	g <- geno[!i, ]
	g$variety <- paste0("var_", 1:nrow(g))
	g <- g[gtools::mixedorder(g$variety), ]

	gg <- geno[i, ]
	gg$longitude <- round(gg$longitude + stats::runif(nrow(gg), -0.5, 0.5), 2) 
	gg$latitude <- round(gg$latitude + stats::runif(nrow(gg), -0.5, 0.5), 2) 

	s <- sample(nrow(gg)/2)
	ss <- grep("undetermined", gg$inventory)
	s <- unique(c(s, ss))

	remid <- gg$sample[s]
	smpls <- unlist(r1[srow,])
	j <- smpls %in% remid

	r1s <- r1[, !j]
	r2s <- r2[, !j[-1]]
	gs <- gg[-s,]
	
	gs$inventory <- as.integer(as.factor(gs$inventory))

	gen <- rbind(g, gs)

	utils::write.csv(gen, file.path(outdir, "DCas00-0000_variety-info.csv"), row.names=FALSE, na="")

	utils::write.table(r1s, file.path(outdir, "DCas00-0000_SNP.csv"), row.names=FALSE, col.names=FALSE, na="", sep=",")
	utils::write.table(r2s, file.path(outdir, "DCas00-0000_SNP_2row.csv"), row.names=FALSE, col.names=FALSE, na="", sep=",")
}

