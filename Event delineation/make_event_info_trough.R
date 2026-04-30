make_event_info_trough <- function(peak_tbl, data_log, par) {
  # FIX B4: use correct column names (Group_First / Group_Last, not *_Peak_Index)
  peak_tbl %>%
    arrange(Species, Group_First, Group_Last) %>%
    group_by(Species) %>%
    group_modify(~ {
      sp <- .y$Species[[1]]
      x  <- data_log[[sp]]
      n  <- length(x)
      df <- arrange(.x, Group_First, Group_Last)
      
      first_pk <- df$Group_First
      last_pk  <- df$Group_Last
      
      map_dfr(seq_len(nrow(df)), function(i) {
        left_limit  <- if (i == 1) 1L else floor((last_pk[i-1] + first_pk[i]) / 2) + 1L
        right_limit <- if (i == nrow(df)) n else floor((last_pk[i] + first_pk[i+1]) / 2)
        
        find_trough_idx <- function(seg_x, offset) {
          if (all(is.na(seg_x))) return(offset + length(seg_x) - 1L)
          offset + which.min(seg_x) - 1L
        }
        
        left_trough <- if (left_limit < first_pk[i]) {
          seg <- x[left_limit:first_pk[i]]
          find_trough_idx(seg, left_limit)
        } else first_pk[i]
        
        right_trough <- if (last_pk[i] < right_limit) {
          seg <- x[last_pk[i]:right_limit]
          find_trough_idx(seg, last_pk[i])
        } else last_pk[i]
        
        event_width <- right_trough - left_trough + 1L
        base_keep   <- pmin(pmax(round(event_width * par$base_prop), par$base_min), par$base_max)
        start_idx   <- max(left_limit,  left_trough  - base_keep)
        end_idx     <- min(right_limit, right_trough + base_keep)
        
        df[i, ] %>%
          mutate(
            Event_ID               = i,
            Global_Event_ID        = paste(sp, i, sep = "_"),
            Left_Search_Limit      = left_limit,
            Right_Search_Limit     = right_limit,
            Left_Trough_Index      = left_trough,
            Right_Trough_Index     = right_trough,
            Left_Trough_Log        = x[left_trough],
            Right_Trough_Log       = x[right_trough],
            Trough_to_Trough_Width = event_width,
            Baseline_Keep          = base_keep,
            Start_Index            = start_idx,
            End_Index              = end_idx,
            Event_Length           = end_idx - start_idx + 1L,
            Increase_Window_Length = Peak_Index - start_idx + 1L,
            Decrease_Window_Length = end_idx - Peak_Index + 1L,
            Event_Start_Date       = data_log$Date[start_idx],
            Event_End_Date         = data_log$Date[end_idx],
            Event_Start_DOY        = data_log$DOY[start_idx],
            Event_End_DOY          = data_log$DOY[end_idx],
            Truncated_Left         = start_idx == 1L,
            Truncated_Right        = end_idx == n,
            Boundary_Based_On_Merged_Group = Group_N_Peaks > 1
          )
      })
    }) %>%
    ungroup()
}