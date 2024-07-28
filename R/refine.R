

refine_reference <- function(dst, lump=.01, split=.05) {
	dcn <- colnames(dst)
	splmp <- matchpoint:::split_lump(dst, lump, split) 
	out <- list(names=data.frame(old=splmp$old, new=splmp$new), pars=data.frame(splmp[c("split", "lump")]))
	nothr <- matchpoint:::min_dist(dst)
	nsame <- matchpoint:::min_self_dist(dst)
	near <- data.frame(nothr, nsame[,2])
	colnames(near)[2:3] <- c("near_other", "near_same")
	out$near <- near
	out
}