add_event_fit_flags <- function(event_info, par) {
  event_info %>%
    mutate(
      Too_Short_Event           = Event_Length < par$min_fit_length,
      Too_Short_Increase_Window = Increase_Window_Length < par$min_phase_length,
      Too_Short_Decrease_Window = Decrease_Window_Length < par$min_phase_length,
      Truncated_Event           = Truncated_Left | Truncated_Right,
      Fit_Status = case_when(
        Too_Short_Event           ~ "excluded_short_event",
        Too_Short_Increase_Window ~ "excluded_short_increase_phase",
        Too_Short_Decrease_Window ~ "excluded_short_decrease_phase",
        TRUE                      ~ "kept_for_spline_LOOCV"
      ),
      Use_For_Spline_LOOCV = Fit_Status == "kept_for_spline_LOOCV"
    )
}
