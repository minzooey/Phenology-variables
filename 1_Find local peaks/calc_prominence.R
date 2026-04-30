#' Compute prominence and trough values for one peak
#'
#' Prominence = peak height minus the highest of the two basin floors.
#' For the global maximum (no higher neighbour in either direction),
#' prominence equals the peak height itself.
#'
#' @param x  Full abundance vector
#' @param p  Index of the focal peak
#' @return Named tibble with prominence and trough metrics
calc_prominence <- function(x, p) {
  n <- length(x)
  h <- x[p]
  
  # Indices of neighbours that are strictly higher than peak p
  left_higher  <- which(x[seq_len(p - 1)] > h)           # already absolute indices
  right_higher <- which(x[(p + 1):n] > h) + p            # FIX B3: convert to absolute
  
  no_higher_left  <- length(left_higher)  == 0
  no_higher_right <- length(right_higher) == 0
  
  left_bound  <- if (no_higher_left)  1L else max(left_higher)
  right_bound <- if (no_higher_right) n  else min(right_higher)
  
  left_base  <- min(x[left_bound:p],  na.rm = TRUE)
  right_base <- min(x[p:right_bound], na.rm = TRUE)
  
  # FIX B1: prominence = h - base (not h <- base)
  prom <- if (no_higher_left && no_higher_right) {
    h
  } else {
    h - max(left_base, right_base)
  }
  
  # FIX B2: operator precedence with explicit parentheses
  left_trough  <- if (p > 1) min(x[(p - 1L):p],       na.rm = TRUE) else h
  right_trough <- if (p < n) min(x[p:(p + 1L)],       na.rm = TRUE) else h
  
  tibble(
    Left_Higher_Index  = if (no_higher_left)  NA_integer_ else left_bound,
    Right_Higher_Index = if (no_higher_right) NA_integer_ else right_bound,
    Trough_Before      = left_trough,
    Trough_After       = right_trough,
    Peak_Prominence    = prom
  )
}