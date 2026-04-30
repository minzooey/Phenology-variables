make_duration_sensitivity <- function(event_info, thresholds = c(14, 21, 28),
                                      min_phase_length = 5) {
  map_dfr(thresholds, function(th) {
    event_info %>%
      mutate(
        Threshold_Days = th,
        Kept = Event_Length >= th &
          Increase_Window_Length >= min_phase_length &
          Decrease_Window_Length >= min_phase_length
      ) %>%
      summarise(
        Threshold_Days      = th,
        Total_Events        = n(),
        Kept_Events         = sum(Kept, na.rm = TRUE),
        Excluded_Events     = sum(!Kept, na.rm = TRUE),
        Kept_Ratio          = Kept_Events / Total_Events,
        Median_Event_Length = median(Event_Length, na.rm = TRUE),
        Min_Event_Length    = min(Event_Length, na.rm = TRUE),
        Max_Event_Length    = max(Event_Length, na.rm = TRUE)
      )
  })
}
