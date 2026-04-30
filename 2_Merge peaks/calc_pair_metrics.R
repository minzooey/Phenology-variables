calc_pair_metrics <- function(x, left_idx, right_idx, par) {
  seg        <- x[left_idx:right_idx]
  trough     <- min(seg, na.rm = TRUE)
  trough_idx <- left_idx + which.min(seg) - 1L
  
  h1 <- x[left_idx]; h2 <- x[right_idx]
  s_peak <- min(h1, h2, na.rm = TRUE)
  l_peak <- max(h1, h2, na.rm = TRUE)
  
  trough_ratio    <- trough / s_peak
  LT              <- par$MTR * s_peak
  LD              <- calc_low_duration(seg, LT)
  recovery_ratio  <- (s_peak - trough) / max(l_peak - trough, .Machine$double.eps)
  peak_dist       <- right_idx - left_idx
  strong_boundary <- trough_ratio <= par$MTR && LD >= par$MDT &&
    peak_dist >= par$MPD && recovery_ratio >= par$MRR
  
  tibble(
    left_idx = left_idx, right_idx = right_idx,
    trough = trough, trough_idx = trough_idx,
    trough_ratio = trough_ratio, low_dur = LD,
    recovery_ratio = recovery_ratio, peak_dist = peak_dist,
    strong_boundary = strong_boundary,
    merge_flag = !strong_boundary,
    merge_score =
      pmin(trough_ratio, 1) +
      (1 - pmin(LD / par$MDT, 1)) +
      (1 - pmin(peak_dist / par$MPD, 1)) +
      (1 - pmin(recovery_ratio / par$MRR, 1))
  )
}