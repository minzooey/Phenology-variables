make_event_data <- function(event_info, data_log, data_raw, keep_only_fit_ready = TRUE) {
  event_tbl <- if (keep_only_fit_ready && "Use_For_Spline_LOOCV" %in% names(event_info)) {
    filter(event_info, Use_For_Spline_LOOCV)
  } else event_info
  
  map_dfr(seq_len(nrow(event_tbl)), function(i) {
    sp  <- event_tbl$Species[i]
    idx <- event_tbl$Start_Index[i]:event_tbl$End_Index[i]
    pk  <- event_tbl$Peak_Index[i]
    
    tibble(
      Species          = sp,
      Event_ID         = event_tbl$Event_ID[i],
      Global_Event_ID  = event_tbl$Global_Event_ID[i],
      Fit_Status       = event_tbl$Fit_Status[i],
      Date             = data_log$Date[idx],
      DOY              = data_log$DOY[idx],
      Index            = idx,
      Rel_Day          = idx - pk,
      Abundance_log    = data_log[[sp]][idx],
      Abundance_raw_sm5 = data_raw[[sp]][idx],
      Peak_Index        = pk,
      Peak_DOY          = event_tbl$Peak_DOY[i],
      Peak_Date         = event_tbl$Peak_Date[i],
      Left_Trough_Index  = event_tbl$Left_Trough_Index[i],
      Right_Trough_Index = event_tbl$Right_Trough_Index[i],
      Event_Start_Index  = event_tbl$Start_Index[i],
      Event_End_Index    = event_tbl$End_Index[i],
      Event_Start_Date   = event_tbl$Event_Start_Date[i],
      Event_End_Date     = event_tbl$Event_End_Date[i]
    )
  })
}