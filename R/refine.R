

refine_reference <- function(dst, lump=.01, split=.05) {
	dcn <- colnames(dst)
	nothr1 <- matchpoint:::min_dist(dst)
	nsame1 <- matchpoint:::min_self_dist(dst)
	splmp <- matchpoint:::split_lump(dst, lump, split) 
	out <- list(varieties=data.frame(old_name=splmp$old, new_name=splmp$new), parameters=data.frame(splmp[c("split", "lump")]))
	nothr2 <- matchpoint:::min_dist(dst)
	nsame2 <- matchpoint:::min_self_dist(dst)
	near <- data.frame(nsame1, nsame2[,2], nothr1[,1], nothr2[,2])
	colnames(near)[2:5] <- c("same_old", "same_new", "other_old", "other_new")
	out$nearest <- near
	out
}