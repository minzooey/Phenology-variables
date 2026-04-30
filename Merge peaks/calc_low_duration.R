calc_low_duration <- function(x, LT) {
  flag <- !is.na(x) & x <= LT
  runs <- rle(flag)
  if (!any(runs$values)) return(0L)
  max(runs$lengths[runs$values])
}