#' Summarise peak candidates for one species
summarize_peak_candidates <- function(df_sp, peak_idx) {
  if (length(peak_idx) == 0) return(tibble())
  
  x <- df_sp$Abundance
  n <- length(x)
  
  prom_df <- map_dfr(peak_idx, ~ calc_prominence(x, .x))
  
  tibble(
    Peak_Index  = peak_idx,
    Peak_DOY    = df_sp$DOY[peak_idx],
    Peak_Date   = df_sp$Date[peak_idx],
    Peak_Height = x[peak_idx]
  ) %>%
    bind_cols(prom_df) %>%
    mutate(
      Prev_Peak_Index = lag(Peak_Index),
      Next_Peak_Index = lead(Peak_Index),
      Dist_from_prev  = Peak_Index - Prev_Peak_Index,
      Dist_to_next    = Next_Peak_Index - Peak_Index,
      Boundary_Left   = Peak_Index == 1L,
      Boundary_Right  = Peak_Index == n
    )
}