

fix_filename <- function(filename) {
	filename <- trimws(filename)
	if (filename != "") {
		filename <- tools::file_path_sans_ext(filename)
		outdir <- dirname(filename)
		dir.create(outdir, FALSE, TRUE)
	}
	filename
}


fix_duplicate_names <- function(x, verbose=FALSE) {
	cn <- colnames(x)
	tab <- table(cn)
	dups <- tab[tab > 1]
	if (length(dups) > 0) {
		if (verbose) message("adding _D1 _D2 to duplicates")
		for (i in 1:length(dups)) {
			n <- dups[i]
			sfx <- paste0("_D", 1:n)
			nm <- names(n)
			cn[cn == nm] <- paste0(nm, sfx)
		}
		colnames(x) <- cn
	}
	x
}

remove_unknown_samples <- function(x, plate.id, verbose=FALSE) {
	cns <- colnames(x)[-1]
	i <- match(cns, plate.id)
	j <- is.na(i)
	if (any(j)) {
		if (verbose) message(paste(sum(j), "unmatched genotypes removed"))
		x <- x[, -(which(j)+1)]
	}
	x
}


get_clean_data <- function(snp, genotypes, markers, filename, verbose=FALSE) {
	snp <- remove_unknown_samples(snp, genotypes$plate.id, verbose=verbose)
	i <- match(colnames(snp)[-1], genotypes$plate.id)
	snp <- fix_duplicate_names(snp, verbose=verbose)
	cns <- colnames(snp)[-1]
	ref.id <- cns[genotypes$reference[i]]
	field.id <- cns[!genotypes$reference[i]]
	markers <- markers[, c("MarkerName", "Chr", "Pos")]
	imark <- match(toupper(snp[,1]), toupper(markers$MarkerName))
	if (any(is.na(imark))) {
		unk <- snp[is.na(imark), 1]
		message("unknown markers in snp file:\n", paste(unk, collapse=", "))
	}
	markers <- markers[imark, ]
	markers$Chr[is.na(markers$Chr)] <- ""
	markers$Pos[is.na(markers$Pos)] <- ""
	filename <- fix_filename(filename)
	if (filename == "") {
		gds <- paste0(tempfile(), "_geno.gds")	
	} else {
		gds <- paste0(filename, "_geno.gds")
	}
	list(snp=snp, markers=markers, ref.id=ref.id, field.id=field.id, filename=filename, gds=gds)
}



match_IBS <- function(SNPs, genotypes, markers, MAF_cutoff=0.05, SNP_Missing_Rate=0.2, 
				Ref_Missing_Rate=0.2, Sample_Missing_Rate=0.2, 
				Ref_Heterozygosity_Rate=0.1, Sample_Heterozygosity_Rate=0.1,
				IBS_cutoff=0.5, Inb_method="mom.visscher", 
				threads=1, verbose=FALSE, filename="") {


	input <- get_clean_data(SNPs, genotypes, markers, filename=filename, verbose=verbose)
	
	meta1 <- data.frame(metric=c("n references", "n samples", "n markers"), 
						value=c(length(input$ref.id), length(input$field.id), nrow(input$snp)))
	
	### recode SNP to be the number of A alleles
	# There are four possible values stored in the variable genotype: 0, 1, 2 and 3.
	# For bi-allelic SNP sites, 
	# 0 = two B alleles
	# 1 = one A and one B allele
	# 2 = two A alleles
	# 3 = missing 
	# For multi-allelic sites, it is a count of the reference allele (3 meaning no call). 
	# “Bit2” indicates that each byte encodes up to four SNP genotypes 
	# since one byte consists of eight bits.

	dmat <- as.matrix(input$snp[, -1])
	dmat[] <- match(dmat, c(0,2,1,NA)) - 1
	
	## create a gds file
	SNPRelate::snpgdsCreateGeno(input$gds, genmat=dmat,
					sample.id = colnames(dmat), 
					snp.id = input$markers$MarkerName,
					snp.chromosome = input$markers$Chr,
					snp.position = input$markers$Pos,
					snp.allele = NULL, snpfirstdim=TRUE)
	
		
	# open the gds file
	genofile <- SNPRelate::snpgdsOpen(input$gds)
	
	# calculate SNP maf and missing rate
	frq <- SNPRelate::snpgdsSNPRateFreq(genofile, sample.id=NULL, snp.id=NULL, with.id=TRUE)
	snpinfo <- data.frame(snpid=frq$snp.id, maf=frq$MinorFreq, mr=frq$MissingRate)
	d0 <- snpinfo[snpinfo$maf >= MAF_cutoff & snpinfo$mr <= SNP_Missing_Rate, ]
	
	# reformatting SNP info
	snpinfo$drp_maf <- (snpinfo$maf < MAF_cutoff) | is.na(snpinfo$maf)
	snpinfo$drp_mr <- snpinfo$mr > SNP_Missing_Rate
	snpinfo$drp <- snpinfo$drp_maf | snpinfo$drp_mr
 	
	# calculate sample missing rate
	sm <- SNPRelate::snpgdsSampMissRate(genofile, sample.id=input$field.id, snp.id=d0$snpid, with.id=TRUE)
	#s0 <- subset(fld, sm <= Sample_Missing_Rate)
	
	rm <- SNPRelate::snpgdsSampMissRate(genofile, sample.id=input$ref.id, snp.id=d0$snpid, with.id=TRUE)

	#r0 <- subset(rf, rm <= Ref_Missing_Rate)
	meta2 <- data.frame(
		metric=c("Marker MAF Cutoff", "Marker Missing Cutoff", "Marker Final", "Marker Final Coverage", "Field Sample Total", "Field Sample Missing Cutoff", "Field Final", "Reference Sample Total", "Reference Sample Missing Cutoff", "Reference Final"), 
		value=c(MAF_cutoff, SNP_Missing_Rate, nrow(d0), ".", length(input$field.id), Sample_Missing_Rate, sum(sm <= Sample_Missing_Rate), length(input$ref.id), Ref_Missing_Rate, sum(sm <= Sample_Missing_Rate))
	)
	
	
	# calculate het rate
	h <- apply(dmat, 2, \(x) sum(x==1) / sum(x!=3))
	h <- data.frame(sid=colnames(dmat), h=h)

	clamp_merge_inb <- function(x, y) {
		# normalize the range between 0 and 1
		x$inb[x$inb < 0] = 0
		x$inb[x$inb > 1] = 1
		x$het <- 1 - x$inb
		ic <- merge(x, h, by="sid")
		merge(y, ic, by.x="row.names", by.y="sid")
	}
	
	drop_rates <- function(x) {
		x$drp_miss <- x$mr > Sample_Missing_Rate
		x$drp_het  <- x$het > Sample_Heterozygosity_Rate
		x$drp <- x$drp_miss | x$drp_het
		x
	}

	# calculate inbreeding coefficient and then heterozygosity for field data
	inb_fld <- SNPRelate::snpgdsIndInb(genofile, sample.id=input$field.id, snp.id=d0$snp.id, autosome.only=FALSE,remove.monosnp=TRUE, maf=NaN,missing.rate=NaN, method=Inb_method, allele.freq=NULL, out.num.iter=TRUE, verbose=verbose)

	ic0 <- data.frame(sid=inb_fld$sample.id, inb=inb_fld$inbreeding)
	fld <- clamp_merge_inb(ic0, data.frame(sm)) 
	fld <- fld[, c("Row.names", "sm", "h")]
	names(fld) <- c("fld_smp_id", "mr", "het")
	fld <- drop_rates(fld)

	# calculate inbreeding coefficient and then heterozygosity for ref data
	inb_ref <- SNPRelate::snpgdsIndInb(genofile, sample.id=input$ref.id, snp.id=d0$snp.id, autosome.only=FALSE, remove.monosnp=TRUE, maf=NaN, missing.rate=NaN, method=Inb_method, allele.freq=NULL, out.num.iter=TRUE, verbose=verbose)
	
	ic1 <- data.frame(sid=inb_ref$sample.id, inb=inb_ref$inbreeding)
	ref <- clamp_merge_inb(ic1, data.frame(rm)) 
	ref <- ref[, c("Row.names", "rm", "h")]
	names(ref) <- c("ref_smp_id", "mr", "het")
	ref <- drop_rates(ref)

	meta3 <- data.frame(
		metric=c("Sample Het Avg", "Sample Het SD", "Sample Het Max", "Sample Het Min", "Ref Heterozygosity Cutoff", "Sample Heterozygosity Cutoff"), 
		value=c(mean(fld$h), stats::sd(fld$h), max(fld$h), min(fld$h), Ref_Heterozygosity_Rate, Sample_Heterozygosity_Rate)
	)
	meta <- rbind(meta1, meta2, meta3)
	output <- list(metadata=meta)
	
	# calculate pair-wise IBS
	ibs <- SNPRelate::snpgdsIBS(genofile, num.thread=threads, autosome.only=FALSE, maf= MAF_cutoff, missing.rate=SNP_Missing_Rate, verbose=verbose)
	SNPRelate::snpgdsClose(genofile)
	
	out1 <- round(ibs[[3]], 4)
	colnames(out1) <- rownames(out1) <- names(input$snp)[-1]
	i <- colnames(out1) %in% input$ref.id
	out2 <- out1[-i, i]
	
	d <- data.frame(field_id=rownames(out2), 
					ref_id=rep(colnames(out2), each=nrow(out2)),
					variety="", IBS=as.vector(out2))

	out2 <- data.frame(FieldSample=colnames(out1)[-i], out2, check.names=FALSE)
	rownames(out2) <- NULL

	mref_id <- gsub("_D.$", "", d$ref_id)
	d$variety <- genotypes$genotype[match(mref_id, genotypes$plate.id)]
	d <- d[with(d, order(field_id, -IBS)),]

# matches
	best <- d[!duplicated(d$field_id),]
	best <- merge(best, sm, by.x="field_id", by.y="row.names", all.y=TRUE)
	names(best)[5] <- "Sample_SNP_Missing_Rate"
	output[["best_match"]] <- best

	IBS_cutoff <- IBS_cutoff[1]
	d <- d[d$IBS > IBS_cutoff, ]
	output[[paste0("IBS_", IBS_cutoff)]] <- d
		
	nr <- as.data.frame(table(field_id=d$field_id))
	rept <- stats::aggregate(d[, "IBS", drop=FALSE], d[, "field_id", drop=FALSE], 
			\(i) c(mean(i, na.rm=TRUE), stats::sd(i, na.rm=TRUE), length(stats::na.omit(i))))
	rept <- data.frame(rept[1], rept[[2]])
	colnames(rept)[-1] <- c("avg", "sd", "nr")

	meta4 <- data.frame(
		metric = paste0(c("Samples with match using IBS=", "Avg number matches per sample using IBS=", "Avg of IBS value with IBS=", "SD of IBS value with IBS="), IBS_cutoff), 
		value = c(nrow(rept), mean(rept$nr), mean(rept$avg, na.rm=TRUE), mean(rept$sd, na.rm=TRUE) ))

	output$meta <- rbind(output$meta, meta4)

	output[["similarity"]] <- out2
	out1 <- data.frame(ID=colnames(out1), 1-out1, check.names=FALSE)
	rownames(out1) <- NULL
	output[["distance"]] <- out1
	output[["snpinfo"]] <- snpinfo
	colnames(ref)[1] <- "ref_id"
	output[["ref"]] <- ref
	colnames(fld)[1] <- "field_id"
	output[["field"]] <- fld

	if (input$filename != "") {
		xlsx <- paste0(input$filename, ".xlsx")
		writexl::write_xlsx(output, xlsx)
		invisible(output)
	} else {
		output$gds <- input$gds
		output
	}
}


