run_peak_merging <- function(peak_candidates, data_log, par) {
  peak_candidates %>%
    group_split(Species) %>%
    map_dfr(function(df) {
      sp <- unique(df$Species)
      if (length(sp) != 1) stop("Each split must contain exactly one species.")
      if (!sp %in% names(data_log)) {
        warning("Species not found in data_log: ", sp); return(tibble())
      }
      merge_peaks_one_species(df, data_log[[sp]], par) %>% mutate(Species = sp)
    }) %>%
    arrange(Species, Group_First)
}