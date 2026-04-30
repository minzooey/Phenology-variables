run_mdt_sensitivity <- function(peak_candidates, data_log, base_par, MDT_cand) {
  map_dfr(MDT_cand, function(mdt) {
    par <- base_par; par$MDT <- mdt
    events <- run_peak_merging(peak_candidates, data_log, par)
    events %>%
      summarise(
        MDT = mdt,
        Total_Events              = n(),
        Mean_Events_per_Species   = mean(table(Species)),
        Median_Events_per_Species = median(table(Species)),
        Total_Merged_Peaks        = sum(Group_N_Peaks - 1),
        Mean_Group_N_Peaks        = mean(Group_N_Peaks)
      )
  }) %>%
    arrange(MDT) %>%
    mutate(Delta_Events = Total_Events - lag(Total_Events),
           Delta_Ratio  = abs(Delta_Events) / lag(Total_Events))
}