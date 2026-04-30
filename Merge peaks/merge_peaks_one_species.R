merge_peaks_one_species <- function(peak_tbl, x, par) {
  peak_tbl <- arrange(peak_tbl, Peak_Index)
  if (nrow(peak_tbl) == 0) return(tibble())
  
  groups <- peak_tbl %>%
    mutate(
      Event_ID     = row_number(),
      Group_First  = Peak_Index,
      Group_Last   = Peak_Index,
      Merged_Peaks = as.character(Peak_ID)
    )
  
  repeat {
    if (nrow(groups) <= 1) break
    
    pair_tbl <- map_dfr(seq_len(nrow(groups) - 1), function(i) {
      calc_pair_metrics(x, groups$Group_Last[i], groups$Group_First[i + 1], par) %>%
        mutate(left_row = i, right_row = i + 1)
    })
    
    merge_pair <- pair_tbl %>%
      filter(merge_flag) %>%
      arrange(desc(merge_score)) %>%
      slice(1)
    
    if (nrow(merge_pair) == 0) break
    
    i <- merge_pair$left_row; j <- merge_pair$right_row
    merged <- bind_rows(groups[i, ], groups[j, ])
    
    new_group <- merged %>%
      arrange(desc(Peak_Height), Peak_Index) %>%
      slice(1) %>%
      mutate(
        Group_First  = min(merged$Group_First),
        Group_Last   = max(merged$Group_Last),
        Merged_Peaks = paste(merged$Merged_Peaks, collapse = ",")
      )
    
    groups <- groups %>%
      slice(-c(i, j)) %>%
      bind_rows(new_group) %>%
      arrange(Group_First) %>%
      mutate(Event_ID = row_number())
  }
  
  groups %>% mutate(
    Event_Peak_ID = row_number(),
    Group_N_Peaks = lengths(strsplit(Merged_Peaks, ","))
  )
}