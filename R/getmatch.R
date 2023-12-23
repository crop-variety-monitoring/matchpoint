


getmatch <- function(df, cutoff=0.9, outprefix="reports/result_c0.9", verbose=TRUE){
  # df: a data.frame contains IBS results
  # cutoff: IBS cutoff for a matching 
  nref <- length(unique(df$FID1))
  nfield <- length(unique(df$FID2))
  sub <- df[df$IBS > cutoff, ]
  out <- sub[, c("FID2", "FID1", "IBS")]
  names(out) <- c("field_id", "ref_id", "IBS")
  out <- out[with(out, order(field_id, -IBS)),]

  #res <- ddply(out, .(field_id), nrow)
	res <- data.frame(table(out$field_id))
  
  
  if (verbose) {
	message(sprintf("# Matching [ n=%s ] field samples with [ n=%s ] ref samples using [ IBS cutoff=%s ] ...", nfield, nref, cutoff))
	message(sprintf("# As a result, [ n=%s/%s ] field samples matched with at least one ref samples (on average matching [n= %s] refs)!", 
      nrow(res), nfield, round(mean(res$Freq), 1)))
  }
  
  if(is.null(outprefix)){
    if (verbose) message(sprintf("# return a list "))
  } else {
    long_res <- paste0(outprefix, "_all_matches.csv")
    data.table::fwrite(out, long_res, sep=",", row.names = FALSE, quote=FALSE)
  }
  
  out1 <- out[!duplicated(out$field_id),]
  
  if(is.null(outprefix)){
    return(list(out, out1))
  }else{
    best_res <- paste0(outprefix, "_best.csv")
    data.table::fwrite(out1, best_res, sep=",", row.names = FALSE, quote=FALSE)
    if (verbose) message(sprintf("# Output files: All matches [ %s ] & Best matches [ %s ]", long_res, best_res))
  }
}

# rescale x from 0 to 1
rescale01 <- function(x) {                              
  # Create user-defined function
  (x - min(x)) / (max(x) - min(x))
}

het_rate <- function(x){
  a <- as.data.frame(table(x))
  b <- merge(data.frame(x=c(0,1,2,3), val=c(0,0,0,0)), a, by="x", all.x=TRUE)
  if(nrow(a) < 4){
    b[is.na(b$Freq), ]$Freq <- 0
  }
  return(b[b$x==1,]$Freq/(b[b$x==0,]$Freq + b[b$x==1,]$Freq + b[b$x==2,]$Freq))
}

ref_field_IBS <- function(ref_file, field_file, MAF_cutoff, SNP_Missing_Rate, 
                          Ref_Missing_Rate, Sample_Missing_Rate, 
                          Ref_Heterozygosity_Rate, Sample_Heterozygosity_Rate,
                          IBS_cutoff, outdir, out_prefix, Inb_method = "mom.visscher", 
		cpus=1, verbose=TRUE){
  # Inb_method=c("mom.weir", "mom.visscher", "mle", "gcta1", "gcta2", "gcta3"),
 
	ref <- data.table::fread(ref_file)
	field <- data.table::fread(field_file)
  
	if (verbose) {
		message(sprintf("# ref data: [ %s ] rows and [ %s ] columns",  nrow(ref), ncol(ref)))
		message(sprintf("# field data: [ %s ] rows and [ %s ] columns",  nrow(field), ncol(field)))
	}
  
	stopifnot(nrow(ref) == nrow(field))
  
	
	# assign UID to ref samples
	ref.id = paste0("ref",1:(ncol(ref)-3), "_", names(ref)[-c(1:3)])
	names(ref) <- c("SNPID", "Chr", "Pos", ref.id)
	
	
	# assign UID to field samples
	field.id = paste0("field",1:(ncol(field)-3), "_", names(field)[-c(1:3)])
	names(field) <- c("SNPID", "Chr", "Pos", field.id)

	if (verbose) {
		message(sprintf("# The first three columns ref are: [ (1) %s, (2) %s, (3) %s ]; \n
           and renamed to [ (1) SNPID, (2) Chr, (3) Pos]",  names(ref)[1], names(ref)[2], names(ref)[3]))
  
		message(sprintf("# The first three columns of the field data are: [ (1) %s, (2) %s, (3) %s ];\n
            and renamed to [ (1) SNPID, (2) Chr, (3) Pos]",  names(field)[1], names(field)[2], names(field)[3]))
		
		message(sprintf("# Assigned unique IDs to both ref and field samples ..."))
		message(sprintf("# Calculating SNP info ..."))
	}
	
  
  
  df9 <- merge(ref, field[, -2:-3], by="SNPID")
  
  meta1 <- data.frame(Metric=c("Marker Total Ref", "Marker Total Sample",  "Marker Common"), 
                      Value=c(nrow(ref), nrow(field), nrow(df9)))
  
  ### recode SNP to be the number of A alleles
  # There are four possible values stored in the variable genotype: 0, 1, 2 and 3. 
  # For bi-allelic SNP sites, “0” indicates two B alleles, “1” indicates one A allele and one B allele, “2” indicates 
  # two A alleles, and “3” is a missing genotype. 
  # For multi-allelic sites, it is a count of the reference allele (3 meaning no call). 
  # “Bit2” indicates that each byte encodes up to four SNP genotypes since one byte consists of eight bits.
  tem <- df9[, -c(1:3)]
  tem[tem == "-"] <- 3 # - to missing
  tem[tem == 2] <- 9 # 2 het to 9 temperately
  tem[tem == 1] <- 2 # 1 homo alt to 2
  tem[tem == 9] <- 1 # 9 het back to 1
  df0 <- cbind(df9[, 1:3], tem)
  
  df <- apply(df0[,-c(1:3)], 2, as.numeric)
  # calculate het rate
  h <- apply(df, 2, het_rate)
  
   dir.create(outdir, FALSE, TRUE)
  
  obj_gds <- file.path(outdir, paste0(out_prefix, "_geno.gds"))
  ## create a gds file
  snpgdsCreateGeno(obj_gds, genmat = df,
                   sample.id = names(df0)[-c(1:3)], snp.id = df0$SNPID,
                   snp.chromosome = df0$Chr,
                   snp.position = df0$Pos,
                   snp.allele = NULL, snpfirstdim=TRUE)
  
  
  # open the gbs obj
  genofile <- snpgdsOpen(obj_gds)
  
  # calculate SNP maf and missing rate
  frq <- snpgdsSNPRateFreq(genofile, sample.id=NULL, snp.id=NULL, with.id=TRUE)
  d <- data.frame(snpid=frq$snp.id, maf=frq$MinorFreq, mr=frq$MissingRate)
  d0 <- subset(d, maf >= MAF_cutoff & mr <= SNP_Missing_Rate)
  
  # reformatting SNP info
  snpinfo <- d
  snpinfo$drp_maf <- 0
  snpinfo$drp_mr <- 0
  snpinfo$drp <- 0
  if(nrow(snpinfo[snpinfo$maf < MAF_cutoff, ]) >0){
    snpinfo[snpinfo$maf < MAF_cutoff, ]$drp_maf <- 1
  }
  if(nrow(snpinfo[snpinfo$mr > SNP_Missing_Rate, ]) >0){
    snpinfo[snpinfo$mr > SNP_Missing_Rate, ]$drp_mr <- 1
  }
  if(nrow(snpinfo[snpinfo$drp_maf >0 | snpinfo$drp_mr >0,])){
    snpinfo[snpinfo$drp_maf >0 | snpinfo$drp_mr >0,]$drp <- 1
  }
  
  # calculate sample missing rate
  sm <- snpgdsSampMissRate(genofile, sample.id=field.id, snp.id=d0$snpid, with.id=TRUE)
  fld <- as.data.frame(sm)
  #s0 <- subset(fld, sm <= Sample_Missing_Rate)
  
  rm <- snpgdsSampMissRate(genofile, sample.id=ref.id, snp.id=d0$snpid, with.id=TRUE)
  rf <- as.data.frame(rm)
  #r0 <- subset(rf, rm <= Ref_Missing_Rate)
  meta2 <- data.frame(Metric=c("Marker MAF Cutoff", "Marker Missing Cutoff", "Marker Final", "Marker Final Coverage", 
                               "Field Sample Total", "Field Sample Missing Cutoff", "Field Final",
                               "Reference Sample Total", "Reference Sample Missing Cutoff", "Reference Final"
                               ), 
                     Value=c(MAF_cutoff, SNP_Missing_Rate, nrow(d0), ".",
                             ncol(field)-3, Sample_Missing_Rate, nrow(subset(fld, sm <= Sample_Missing_Rate)),
                             ncol(ref)-3, Ref_Missing_Rate, nrow(subset(fld, sm <= Sample_Missing_Rate))))
  
  # calculate inbreeding coefficient and then heterozygosity for field data
  inb <- snpgdsIndInb(genofile, sample.id=field.id, snp.id=d0$snp.id, autosome.only=FALSE,remove.monosnp=TRUE,
                      maf=NaN,missing.rate=NaN, 
                      method=Inb_method, 
                      allele.freq=NULL, out.num.iter=TRUE,verbose=TRUE)
  
  # normalize the range between 0 and 1
  ic <- data.frame(sid=inb$sample.id, inb=inb$inbreeding)
  rg <- range(range(ic$inb))
  if(rg[1] < 0){
    ic[inb$inb < 0, ]$inb = 0
  }
  if(rg[2] > 1){
    ic[inb$inb >1, ]$inb = 1
  }
  ic$het <- 1 - ic$inb
  #range(ic$inb)
  h1 <- as.data.frame(h)
  h1$sid <- row.names(h1)
  ic <- merge(ic, h1, by="sid")
  fld <- merge(fld, ic, by.x="row.names", by.y="sid")
  fld <- fld[, c("Row.names", "sm", "h")]
  names(fld) <- c("fld_smp_id", "mr", "het")
  fld$drp_miss <- 0
  fld$drp_het <- 0
  if(nrow(fld[fld$mr > Sample_Missing_Rate, ]) > 0){
    fld[fld$mr > Sample_Missing_Rate, ]$drp_miss <- 1
  }
  if(nrow(fld[fld$het > Sample_Heterozygosity_Rate, ]) > 0){
    fld[fld$het > Sample_Heterozygosity_Rate, ]$drp_miss <- 1
  }
  if(nrow(fld[fld$drp_miss >0 | fld$drp_het >0,])){
    fld[fld$drp_miss >0 | fld$drp_het >0,]$drp <- 1
  }
  
  # calculate inbreeding coefficient and then heterozygosity for ref data
  inb_ref <- snpgdsIndInb(genofile, sample.id=ref.id, snp.id=d0$snp.id, autosome.only=FALSE,remove.monosnp=TRUE,
                      maf=NaN, missing.rate=NaN, 
                      method=Inb_method, 
                      allele.freq=NULL, out.num.iter=TRUE,verbose=TRUE)
  
  # normalize the range between 0 and 1
  ic1 <- data.frame(sid=inb_ref$sample.id, inb=inb$inbreeding)
  rg <- range(range(ic1$inb))
  if(rg[1] < 0){
    ic[inb$inb < 0, ]$inb = 0
  }
  if(rg[2] > 1){
    ic[inb$inb >1, ]$inb = 1
  }
  ic$het <- 1 - ic$inb
  #range(ic$inb)
  h1 <- as.data.frame(h)
  h1$sid <- row.names(h1)
  ic1 <- merge(ic1, h1, by="sid")
  rf <- merge(rf, ic1, by.x="row.names", by.y="sid")
  rf <- rf[, c("Row.names", "sm", "h")]
  names(rf) <- c("ref_smp_id", "mr", "het")
  rf$drp_miss <- 0
  rf$drp_het <-  0
  if(nrow(rf[rf$mr > Ref_Missing_Rate, ]) > 0){
    rf[rf$mr > Ref_Missing_Rate, ]$drp_miss <- 1
  }
  if(nrow(rf[rf$het > Ref_Heterozygosity_Rate, ]) > 0){
    rf[rf$het > Ref_Heterozygosity_Rate, ]$drp_miss <- 1
  }
  if(nrow(rf[rf$drp_miss >0 | rf$drp_het >0,])){
    rf[rf$drp_miss >0 | rf$drp_het >0,]$drp <- 1
  }
  meta3 <- data.frame(Metric=c("Sample Het Avg", "Sample Het SD", "Sample Het Max", "Sample Het Min",
                               "Ref Heterozygosity Cutoff", "Sample Heterozygosity Cutoff"), 
                      Value=c(mean(ic1$h), sd(ic1$h), max(ic1$h), min(ic1$h),
                              Ref_Heterozygosity_Rate, Sample_Heterozygosity_Rate))
  
  # calculate pair-wise IBS
  ibs <- snpgdsIBS(genofile, num.thread=cpus, autosome.only=FALSE, maf= MAF_cutoff, missing.rate=SNP_Missing_Rate)
  out <- ibs[[3]]
  dimnames(out) <- dimnames(out) <- list(names(df0)[-c(1:3)], names(df0)[-c(1:3)])
  xy <- t(combn(colnames(out), 2))
  
  # convert to data.frame
  d <- data.frame(xy, dist=out[xy])
  names(d) <- c("ref", "field", "IBS")
  
  # subset ref vs. field
  d$type1 <-  "field"
  d$type2 <- "ref"
  d[d$ref %in% names(ref)[-c(1:3)], ]$type1 <- "ref"
  d[d$field %in% names(field)[-c(1:3)], ]$type2 <- "field"
  
  d <- subset(d, type1 == "ref" & type2 =="field")
  
  names(d)[1:3] <- c("FID1", "FID2", "IBS")
  
  meta <- rbind(meta1, meta2, meta3)
  outlist <- list()
  outlist[["metadata"]] <- meta
  #iidx = 1
  for(i in IBS_cutoff){
    lout <- getmatch(d, cutoff=i, outprefix=NULL, verbose=verbose)
    
    ### best match
    best_match <- paste0("IBS_cutoff_", i, "_best_match")
    lout2 <- merge(lout[[2]], sm, by.x="field_id", by.y="row.names", all.y=T)
    names(lout2)[4] <- "Sample_SNP_Missing_Rate"
    outlist[[best_match]] <- lout2
    #iidx <- iidx + 1
    
    all_match <- paste0("IBS_cutoff_", i, "_all_match")
    outlist[[all_match]] <- lout[[1]]
    
    report <- plyr::ddply(lout[[1]], .(field_id), summarise,
                    avg = mean(IBS, na.rm=TRUE),
                    sd = sd(IBS, na.rm=TRUE))
    nr <- plyr::ddply(lout[[1]], .(field_id), nrow)
    meta4 <- data.frame(Metric=c(paste0("Samples with match using IBS=",i), paste0("Avg number matches per sample using IBS=",i), paste0("Avg of IBS value with IBS=", i), paste0("SD of IBS value with IBS=", i)), 
       Value=c(nrow(report), mean(nr$V1), mean(report$avg, na.rm = TRUE), mean(report$sd, na.rm = TRUE) ))
    meta <- rbind(meta, meta4)
  }
  outlist[["metadata"]] <- meta
  outlist[["ref"]] <- rf
  outlist[["fld"]] <- fld
  outlist[["snpinfo"]] <- snpinfo

  snpgdsClose(genofile)
  
  #file.remove(obj_gds)
  
  xlsx <- file.path(outdir, paste0(out_prefix, "_matches.xlsx"))
  writexl::write_xlsx(outlist, xlsx)
  if (verbose) message(sprintf("# Analysis Finished !!! results can be found in excel file [ %s ]", xlsx))
}
