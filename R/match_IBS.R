

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

remove_unknown_samples <- function(x, sample.id, verbose=FALSE) {
	cns <- colnames(x)[-1]
	i <- match(cns, sample.id)
	j <- is.na(i)
	if (any(j)) {
		if (verbose) message(paste(sum(j), "unmatched genotypes removed"))
		x <- x[, -(which(j)+1)]
	}
	x
}


get_clean_data <- function(snp, genotypes, markers, filename, verbose=FALSE) {
	snp <- remove_unknown_samples(snp, genotypes$sample, verbose=verbose)
	i <- match(colnames(snp)[-1], genotypes$sample)
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
				Ref_Heterozygosity_Rate=1, Sample_Heterozygosity_Rate=1,
				IBS_cutoff=0.5, Inb_method="mom.visscher", 
				threads=1, verbose=FALSE, filename="") {


	input <- matchpoint:::get_clean_data(SNPs, genotypes, markers, filename=filename, verbose=verbose)
	
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
	snpinfo$drop_maf <- (snpinfo$maf < MAF_cutoff) | is.na(snpinfo$maf)
	snpinfo$drop_mr <- snpinfo$mr > SNP_Missing_Rate
	snpinfo$drop <- snpinfo$drop_maf | snpinfo$drop_mr
 	
	# calculate sample missing rate
	smp_mr <- SNPRelate::snpgdsSampMissRate(genofile, sample.id=input$field.id, snp.id=d0$snpid, with.id=TRUE)
	
	ref_mr <- SNPRelate::snpgdsSampMissRate(genofile, sample.id=input$ref.id, snp.id=d0$snpid, with.id=TRUE)

	meta2 <- data.frame(
		metric=c("IBS_cutoff", "MAF_cutoff", "SNP_Missing_Raste Cutoff", "Marker Final", "Field Sample Total", "Field Sample Missing Cutoff", "Field Final", "Reference Sample Total", "Reference Sample Missing Cutoff", "Reference Final"), 
		value=c(IBS_cutoff, MAF_cutoff, SNP_Missing_Rate, nrow(d0), length(input$field.id), Sample_Missing_Rate, sum(smp_mr <= Sample_Missing_Rate), length(input$ref.id), Ref_Missing_Rate, sum(ref_mr <= Sample_Missing_Rate))
	)

	# calculate inbreeding coefficient and then heterozygosity for field data
	inb_fld <- SNPRelate::snpgdsIndInb(genofile, sample.id=input$field.id, snp.id=d0$snp.id, autosome.only=FALSE,remove.monosnp=TRUE, maf=NaN,missing.rate=NaN, method=Inb_method, allele.freq=NULL, out.num.iter=TRUE, verbose=verbose)

	# calculate inbreeding coefficient and then heterozygosity for ref data
	inb_ref <- SNPRelate::snpgdsIndInb(genofile, sample.id=input$ref.id, snp.id=d0$snp.id, autosome.only=FALSE, remove.monosnp=TRUE, maf=NaN, missing.rate=NaN, method=Inb_method, allele.freq=NULL, out.num.iter=TRUE, verbose=verbose)
	

	# calculate het rate
	h <- apply(dmat, 2, \(x) sum(x==1) / sum(x!=3))
	h <- data.frame(sid=colnames(dmat), h=h)

	clamp_merge_drop <- function(het, sid, inb, mr, mr_thold, het_thold, name) {
		x <- data.frame(sid=sid, inb=inb)
		mr <- data.frame(mr=mr)
		# normalize the range between 0 and 1
		x$inb[x$inb < 0] = 0
		x$inb[x$inb > 1] = 1
		x$het <- 1 - x$inb
		ic <- merge(x, het, by="sid")
		out <- merge(mr, ic, by.x="row.names", by.y="sid")
		out <- out[, c("Row.names", "mr", "h")]
		names(out)[c(1,3)] <- c(name, "het")
		out$drop_miss <- out$mr > mr_thold
		out$drop_het <- out$het > het_thold
		out$drop <- out$drop_miss | out$drop_het
		out
	}
	
	fld <- clamp_merge_drop(h, inb_fld$sample.id, inb_fld$inbreeding, smp_mr, 
			Sample_Missing_Rate, Sample_Heterozygosity_Rate, "fld_id") 
	ref <- clamp_merge_drop(h, inb_ref$sample.id, inb_ref$inbreeding, ref_mr, 
			Ref_Missing_Rate, Ref_Heterozygosity_Rate, "ref_id") 

	meta3 <- data.frame(
		metric=c("Sample Het Avg", "Sample Het SD", "Sample Het Max", "Sample Het Min", "Ref Heterozygosity Cutoff", "Sample Heterozygosity Cutoff"), 
		value=c(mean(fld$het), 	stats::sd(fld$het), max(fld$het), min(fld$het), Ref_Heterozygosity_Rate, Sample_Heterozygosity_Rate)
	)
	meta <- rbind(meta1, meta2, meta3)
	output <- list(metadata=meta)
	
	# calculate pair-wise IBS

### RH instead of
#	ibs <- SNPRelate::snpgdsIBS(genofile, num.thread=threads, autosome.only=FALSE,
#		maf= MAF_cutoff, missing.rate=.5, verbose=verbose)

### RH omitting genotypes / markers with too much missing data 
	use_ids <- c(fld$fld_id[!fld$drop_miss], ref$ref_id[!ref$drop_miss])
	ibs <- SNPRelate::snpgdsIBS(genofile, sample.id=use_ids, 
		num.thread=threads, autosome.only=FALSE, verbose=verbose,
		maf=MAF_cutoff, missing.rate=SNP_Missing_Rate)


	SNPRelate::snpgdsClose(genofile)

	out_all <- round(ibs[[3]], 4)
	colnames(out_all) <- rownames(out_all) <- ibs[[1]]
	i <- which(colnames(out_all) %in% input$ref.id)
	out_match <- out_all[-i, i]
	
	d <- data.frame(field_id=rownames(out_match), 
					ref_id=rep(colnames(out_match), each=nrow(out_match)),
					variety="", IBS=as.vector(out_match))

	out_match <- data.frame(FieldSample=colnames(out_all)[-i], out_match, check.names=FALSE)
	rownames(out_match) <- NULL

	mref_id <- gsub("_D.$", "", d$ref_id)
	d$variety <- genotypes$variety[match(mref_id, genotypes$sample)]
	d <- d[with(d, order(field_id, -IBS)),]

# matches
	best <- d[!duplicated(d$field_id),]
	best <- merge(best, smp_mr, by.x="field_id", by.y="row.names", all.y=TRUE)
	names(best)[5] <- "Sample_SNP_Missing_Rate"
	output[["best_match"]] <- best

	output[[paste0("IBS")]] <- d[d$IBS > IBS_cutoff[1], ]
		
	nr <- as.data.frame(table(field_id=d$field_id))
	rept <- stats::aggregate(d[, "IBS", drop=FALSE], d[, "field_id", drop=FALSE], 
			\(i) c(mean(i, na.rm=TRUE), stats::sd(i, na.rm=TRUE), length(stats::na.omit(i))))
	rept <- data.frame(rept[1], rept[[2]])
	colnames(rept)[-1] <- c("avg", "sd", "nr")

	meta4 <- data.frame(
		metric = paste0(c("Samples with match using IBS=", "Avg number matches per sample using IBS=", "Avg of IBS value with IBS=", "SD of IBS value with IBS="), IBS_cutoff), 
		value = c(nrow(rept), mean(rept$nr), mean(rept$avg, na.rm=TRUE), mean(rept$sd, na.rm=TRUE) ))

	output$metadata <- rbind(output$metadata, meta4)
	output$metadata$value <- round(output$metadata$value, 5)
	
	output[["similarity"]] <- out_match
	out_all <- data.frame(ID=colnames(out_all), 1-out_all, check.names=FALSE)
	rownames(out_all) <- NULL
	output[["distance"]] <- out_all
	output[["snpinfo"]] <- snpinfo

	ref$variety <- genotypes$variety[match(gsub("_D.$", "", ref$ref_id), genotypes$sample)]
	output[["ref"]] <- ref
	output[["field"]] <- fld

	if (input$filename != "") {
		xlsx <- paste0(input$filename, "_IBS.xlsx")
		writexl::write_xlsx(output, xlsx, format_headers=FALSE)
		saveRDS(ibs, paste0(input$filename, "_IBS.rds"))
		invisible(output)
	} else {
		output$gds <- input$gds
		output
	}
}


