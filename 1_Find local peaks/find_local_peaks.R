#' Find raw local maxima, including boundary peaks
#' @param x     Numeric vector of (log-transformed) abundance
#' @param min_height  Minimum height threshold (absolute, in log units)
#' @return Integer vector of peak indices
find_local_peaks <- function(x, min_height) {
  n <- length(x)
  if (n < 2) return(integer(0))
  
  peak_idx <- integer(0)
  
  # left boundary
  if (!is.na(x[1]) && !is.na(x[2]) && x[1] >= x[2] && x[1] >= min_height)
    peak_idx <- c(peak_idx, 1L)
  
  # interior peaks (strict increase on left, non-strict on right to handle plateaus)
  if (n >= 3) {
    inner <- which(
      x[2:(n - 1)] > x[1:(n - 2)] &
        x[2:(n - 1)] >= x[3:n] &
        x[2:(n - 1)] >= min_height
    ) + 1L
    peak_idx <- c(peak_idx, inner)
  }
  
  # right boundary
  if (!is.na(x[n]) && !is.na(x[n - 1]) && x[n] > x[n - 1] && x[n] >= min_height)
    peak_idx <- c(peak_idx, n)
  
  sort(unique(peak_idx))
}