#' Detect peak candidates for one species
detect_peak_candidates <- function(data_log, data_raw, sp_name, par) {
  df_sp <- data_log %>%
    transmute(Date, DOY, Abundance = as.numeric(.data[[sp_name]]))
  
  x      <- df_sp$Abundance
  x_max  <- max(x, na.rm = TRUE)
  x_range <- diff(range(x, na.rm = TRUE))
  
  if (!is.finite(x_max) || x_max <= 0 || x_range <= 0) return(tibble())
  
  min_height <- par$MHR * x_max
  min_prom   <- par$MPR * x_range
  
  peak_idx <- find_local_peaks(x, min_height)
  
  summarize_peak_candidates(df_sp, peak_idx) %>%
    mutate(
      Species              = sp_name,
      Peak_ID              = row_number(),
      Peak_Height_Raw_SM5  = data_raw[[sp_name]][Peak_Index],
      Max_Log_Species      = x_max,
      Range_Log_Species    = x_range,
      Min_Height_Log       = min_height,
      Min_Prom_Log         = min_prom
    ) %>%
    dplyr::select(
      Species, Peak_ID,
      Peak_Index, Peak_DOY, Peak_Date,
      Peak_Height, Peak_Height_Raw_SM5,
      Peak_Prominence,
      Dist_from_prev, Dist_to_next,
      Trough_Before, Trough_After,
      Boundary_Left, Boundary_Right,
      Max_Log_Species, Range_Log_Species, Min_Height_Log, Min_Prom_Log
    )
}