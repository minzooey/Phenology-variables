## =============================================================================
## Phenology of Dinoflagellates — Jangcheon Harbor (411-day fixed-station)
## =============================================================================
## Description : Computes cardinal phenology variables (DBS, DMF, DMA, DMM, DBE
##               and derived metrics) for 32 dinoflagellate taxa from daily
##               5-day smoothed abundance time-series.
##
## Pipeline    : 1) Data import & log10 transformation
##               2) Local peak detection (find_local_peaks)
##               3) Peak merging via MDT-based score (merge_peaks_one_species)
##               4) Event delineation with trough-to-trough windows
##               5) Cubic smoothing spline + residual bootstrap CI (LOOCV / fixed-spar)
##               6) Cardinal phenology variable extraction
##
## Inputs      : JC_abundance.xlsx   — daily raw counts per species
##               JC_envs_daily.xlsx  — daily environmental covariates
##               dino_5day.xlsx      — 5-day smoothed abundances (pre-processed)
##
## Outputs     : OUT_DIR/
##                 1_local_peak_candidates.xlsx
##                 2_peak_merging_results_MDT_sensitivity.xlsx
##                 3_event_specific_dataset_filtered.xlsx
##                 4_event_spline_LOOCV_fitting.xlsx
##                 5_phenology_variables.xlsx
##
## Author      : Juhee, Min
## Affiliation : Department of Oceanography, Chonnam National University, Gwangju, Korea, 61186
## Contact     : minzooey@gmail.com
## Last updated: 2026-04-20
## License     : MIT
## =============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

# ── 0) Packages ───────────────────────────────────────────────────────────────
pkgs <- c(
  "openxlsx", "dplyr", "tibble", "purrr", "tidyr",   # data wrangling
  "ggplot2", "gridExtra",                              # visualisation
  "lubridate", "pracma"                                # date / math helpers
)

new_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(new_pkgs)) install.packages(new_pkgs)
invisible(lapply(pkgs, library, character.only = TRUE))

# ── 1) Paths ──────────────────────────────────────────────────────────────────
INPUT_DIR <- "/Users/minjuhee/Desktop/HPLC/7_Jangcheon/3_Phenology"          # change to your working directory
ENVS_FILE <- "JC_envs_daily.xlsx"
ABUN_FILE <- "JC_abundance.xlsx"
SMOO_FILE <- "dino_5day.xlsx"   # 5-day Gaussian-smoothed abundance

OUT_DIR <- file.path(INPUT_DIR, "ONE", paste0("Output_", format(Sys.Date(), "%y%m%d")))
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# ── 2) Load & align data ──────────────────────────────────────────────────────
dino_daily <- read.xlsx(file.path(INPUT_DIR, ABUN_FILE), sheet = "Dino_daily") %>%
  mutate(Date = make_date(Year, Month, Day)) %>%
  relocate(Date, .before = "Alex") %>%
  dplyr::select(-any_of("Undefined"))

# Species columns: all numeric columns except Year/Month/Day/Date
SP_COLS <- names(dino_daily)[sapply(dino_daily, is.numeric) &
                               !names(dino_daily) %in% c("Year", "Month", "Day")]

dino_sm5 <- read.xlsx(file.path(INPUT_DIR, SMOO_FILE)) %>%
  mutate(
    Date = make_date(Year, Month, Day),
    # DOY: days since first observation (1-indexed), robust to date gaps
    DOY  = as.integer(Date - min(Date, na.rm = TRUE)) + 1L
  ) %>%
  dplyr::select(-any_of("Undefined")) %>%
  na.omit()

# Verify both data frames share the same species columns
sp_cols_sm5 <- intersect(SP_COLS, names(dino_sm5))
if (!all(SP_COLS %in% names(dino_sm5))) {
  warning("Some species in dino_daily are absent from dino_sm5. Using intersection.")
}
SP_COLS <- sp_cols_sm5

## — 2a) Long-format for plotting (raw smoothed counts)
dino_long <- dino_sm5 %>%
  pivot_longer(all_of(SP_COLS), names_to = "Species", values_to = "Abundance") %>%
  mutate(Abundance = replace_na(Abundance, 0))

p_raw <- ggplot(dino_long, aes(x = DOY, y = Abundance)) +
  annotate("rect", xmin = 1,   xmax = 46,  ymin = -Inf, ymax = Inf,
           fill = "#d7e8bc", alpha = 0.6) +
  annotate("rect", xmin = 366, xmax = 411, ymin = -Inf, ymax = Inf,
           fill = "#d7e8bc", alpha = 0.6) +
  geom_line(color = "blue", linewidth = 0.6) +
  geom_vline(xintercept = 274, color = "red", linetype = "dashed", linewidth = 0.6) +
  facet_wrap(~Species, scales = "free_y", ncol = 6) +
  labs(x = NULL, y = expression("Abundance (cells L"^{-1}*")"),
       title = "Daily abundance of dinoflagellate species (5-day smoothing)") +
  theme_bw(base_size = 12) +
  theme(strip.background = element_rect(fill = "grey95", color = "black"),
        strip.text   = element_text(face = "bold"),
        axis.text.x  = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank())
print(p_raw)
ggsave(file.path(OUT_DIR,"daily abun.pdf"), plot = p_raw, width = 12, height = 7)

## — 2b) Log10-transform species columns only
dino_sm5_log10 <- dino_sm5 %>%
  mutate(across(all_of(SP_COLS), ~ log10(. + 1)))

# ── 3) Helper functions: local peak detection ──────────────────────────────────
functs1 <- list.files(path = file.path(INPUT_DIR, "ONE/1_Find local peaks"), pattern = "\\.R$", full.names = TRUE, ignore.case = TRUE)

lapply(functs1, source)

# ── 4) Detect peaks for all species ───────────────────────────────────────────
# Include species with ≥ 5 positive observations
sp_names <- names(which(
  colSums(dino_daily[, SP_COLS] > 0, na.rm = TRUE) >= 5
))

peak_pars <- list(
  MHR = 0.15,   # minimum peak height: 15 % of species-specific max
  MPR = 0.15,   # minimum peak prominence: 15 % of species-specific range
  MPD = 15      # minimum peak distance: 15 days (used in sensitivity only)
)

peak_candidates <- map_dfr(sp_names,
  ~ detect_peak_candidates(dino_sm5_log10, dino_sm5, .x, peak_pars)
)

write.xlsx(peak_candidates,
  file.path(OUT_DIR, "1_local_peak_candidates.xlsx"), overwrite = TRUE)

## — visual check
dino_log_long <- dino_sm5_log10 %>%
  pivot_longer(all_of(sp_names), names_to = "Species", values_to = "LogAbundance")

p_peaks <- ggplot(dino_log_long, aes(x = Date, y = LogAbundance)) +
  annotate("rect", xmin = as.Date("2020-04-02"), xmax = as.Date("2020-05-17"),
           ymin = -Inf, ymax = Inf, fill = "#d7e8bc", alpha = 0.6) +
  annotate("rect", xmin = as.Date("2021-03-03"), xmax = as.Date("2021-05-17"),
           ymin = -Inf, ymax = Inf, fill = "#d7e8bc", alpha = 0.6) +
  geom_line(color = "blue", linewidth = 0.5) +
  geom_point(data = peak_candidates,
             aes(x = Peak_Date, y = Peak_Height), color = "red", size = 1.8) +
  facet_wrap(~Species, scales = "free_y", ncol = 6) +
  labs(x = "Date", y = expression(log[10](x + 1)),
       title = "Step 1. Local peak candidates") +
  theme_bw(base_size = 12) +
  theme(strip.background = element_rect(fill = "grey95", color = "black"),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank())
print(p_peaks)
ggsave(file.path(OUT_DIR,"local peaks.pdf"), plot = p_peaks, width = 12, height = 7)

# ── 5) Peak merging (MDT-based score) ─────────────────────────────────────────
peak_merg <- list(
  MTR       = 0.50,
  MPD       = 15,
  MRR       = 0.50,
  MIN_EVENT = 7
)
MDT_cand <- c(3, 5, 7, 10, 14, 21)

## — load useful functions
functs2 <- list.files(path = file.path(INPUT_DIR, "ONE/2_Merge peaks"), pattern = "\\.R$", full.names = TRUE, ignore.case = TRUE)

lapply(functs2, source)

## — sensitivity analysis
mdt_sens <- run_mdt_sensitivity(peak_candidates, dino_sm5_log10, peak_merg, MDT_cand)
print(mdt_sens)

p_mdt <- ggplot(mdt_sens, aes(x = MDT, y = Total_Events)) +
  geom_line() + geom_point(size = 2) +
  theme_bw(base_size = 14) +
  labs(x = "MDT (days)", y = "N detected bloom events",
       title = "Sensitivity of event detection to MDT")
print(p_mdt)
ggsave(file.path(OUT_DIR,"MDT_sensitivity.pdf"), plot = p_mdt, width = 6, height = 4)

## — apply chosen MDT
peak_merg_f       <- peak_merg
peak_merg_f$MDT   <- 10

final_event_peaks <- run_peak_merging(peak_candidates, dino_sm5_log10, peak_merg_f)

## — unpack merged peak IDs
event_peak_map <- final_event_peaks %>%
  dplyr::select(Species, Event_Peak_ID, Merged_Peaks) %>%
  # FIX W3: replace deprecated separate_rows()
  tidyr::separate_longer_delim(Merged_Peaks, delim = ",") %>%
  mutate(Peak_ID = as.integer(Merged_Peaks)) %>%
  dplyr::select(Species, Peak_ID, Event_Peak_ID)

representative_peak_map <- final_event_peaks %>%
  dplyr::select(Species, Event_Peak_ID, Representative_Peak_ID = Peak_ID)

all_peak_results <- peak_candidates %>%
  left_join(event_peak_map, by = c("Species", "Peak_ID")) %>%
  left_join(representative_peak_map, by = c("Species", "Event_Peak_ID")) %>%
  mutate(
    Representative = Peak_ID == Representative_Peak_ID,
    Peak_Status    = ifelse(Representative, "Final", "Merged")
  ) %>%
  arrange(Species, Peak_Index)

final_event_peaks_out <- final_event_peaks %>%
  mutate(Representative = TRUE, Peak_Status = "Final") %>%
  arrange(Species, Group_First)

removed_peaks <- filter(all_peak_results, Peak_Status == "Merged")

openxlsx::write.xlsx(
  list(MDT_sensitivity  = mdt_sens,
       Peak_candidates  = peak_candidates,
       All_peaks        = all_peak_results,
       Final_event_peaks = final_event_peaks_out,
       Removed_peaks    = removed_peaks),
  file = file.path(OUT_DIR, "2_peak_merging_results_MDT_sensitivity.xlsx"),
  overwrite = TRUE
)

## — visualise final peaks
p_merged <- ggplot(dino_log_long, aes(x = Date, y = LogAbundance)) +
  annotate("rect", xmin = as.Date("2020-04-02"), xmax = as.Date("2020-05-17"),
           ymin = -Inf, ymax = Inf, fill = "#d7e8bc", alpha = 0.6) +
  annotate("rect", xmin = as.Date("2021-03-03"), xmax = as.Date("2021-05-17"),
           ymin = -Inf, ymax = Inf, fill = "#d7e8bc", alpha = 0.6) +
  geom_line(color = "blue", linewidth = 0.5) +
  geom_point(data = removed_peaks,
             aes(x = Peak_Date, y = Peak_Height), color = "grey60", size = 1.5) +
  geom_point(data = final_event_peaks_out,
             aes(x = Peak_Date, y = Peak_Height), color = "red", size = 2) +
  facet_wrap(~Species, scales = "free_y", ncol = 6) +
  labs(x = "Date", y = expression(log[10](x + 1)),
       title = "Step 2. Final event peaks after MDT-based score merging") +
  theme_bw(base_size = 12) +
  theme(strip.background = element_rect(fill = "grey95", color = "black"),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank())
print(p_merged)
ggsave(file.path(OUT_DIR,"merged peaks.pdf"), plot = p_merged, width = 12, height = 7)

# ── 6) Event delineation (trough-to-trough windows) ──────────────────────────
event_par <- list(
  base_prop          = 0.25,
  base_min           = 15,
  base_max           = 30,
  min_fit_length     = 21,
  min_phase_length   = 5,
  length_sensitivity = c(14, 21, 28)
)

## — load useful functions
functs3 <- list.files(path = file.path(INPUT_DIR, "ONE/3_Event delineation"), pattern = "\\.R$", full.names = TRUE, ignore.case = TRUE)

lapply(functs3, source)

## — FIX B5: use final_event_peaks_out (not undefined final_event_peaks_re)
event_info_all <- make_event_info_trough(final_event_peaks_out, dino_sm5_log10, event_par) %>%
  add_event_fit_flags(event_par)

event_data_for_fit         <- make_event_data(event_info_all, dino_sm5_log10, dino_sm5, TRUE)
event_duration_sensitivity <- make_duration_sensitivity(
  event_info_all, event_par$length_sensitivity, event_par$min_phase_length)

write.xlsx(
  list(event_info_all             = event_info_all,
       event_info_for_fit         = filter(event_info_all, Use_For_Spline_LOOCV),
       event_info_excluded        = filter(event_info_all, !Use_For_Spline_LOOCV),
       event_data_for_fit         = event_data_for_fit,
       event_duration_sensitivity = event_duration_sensitivity),
  file.path(OUT_DIR, "3_event_specific_dataset_filtered.xlsx"), overwrite = TRUE
)

## — visual check
p_event_check <- ggplot(event_data_for_fit, aes(x = Date, y = Abundance_log)) +
  geom_line(color = "blue", linewidth = 0.5) +
  geom_point(
    data = filter(event_info_all, Use_For_Spline_LOOCV),
    aes(x = Peak_Date, y = Peak_Height), color = "red", size = 2, inherit.aes = FALSE
  ) +
  geom_vline(
    data = filter(event_info_all, Use_For_Spline_LOOCV),
    aes(xintercept = Event_Start_Date), linetype = "dotted", color = "grey40", inherit.aes = FALSE
  ) +
  geom_vline(
    data = filter(event_info_all, Use_For_Spline_LOOCV),
    aes(xintercept = Event_End_Date), linetype = "dotted", color = "grey40", inherit.aes = FALSE
  ) +
  facet_wrap(~Global_Event_ID, scales = "free_y", ncol = 8) +
  labs(x = NULL, y = expression(log[10](x + 1)),
       title = "Step 3. Event-specific datasets kept for Spline + LOOCV") +
  theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey95", color = "black"),
        strip.text = element_text(face = "bold", size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank())
print(p_event_check)
ggsave(file.path(OUT_DIR,"event-specific peaks.pdf"), plot = p_event_check, width = 12, height = 7)

# ── 7) Cubic smoothing spline + residual bootstrap CI ─────────────────────────
spline_par <- list(
  min_n          = 8,
  n_boot         = 200,
  fit_mode       = "fixed_spar",   # "loocv" or "fixed_spar"
  spar_fixed     = 0.9,
  spar_min       = 0.55,
  spar_max       = 1.20,
  peak_diff_cut  = 0.5,
  baseline_prop  = 0.15,
  baseline_min_n = 3
)

## — load useful functions
functs4 <- list.files(path = file.path(INPUT_DIR, "ONE/4_Fit & phenology variables"), pattern = "\\.R$", full.names = TRUE, ignore.case = TRUE)

lapply(functs4, source)

#' Fit cubic smoothing spline with residual bootstrap CI for one event
#'
#' Efficiency: bootstrap matrix filled with vapply for type-safety;
#'             CI computed in one apply pass.
event_data_input <- if (exists("event_data_for_fit")) event_data_for_fit else event_data
if (!"Global_Event_ID" %in% names(event_data_input))
  event_data_input <- mutate(event_data_input, Global_Event_ID = paste(Species, Event_ID, sep = "_"))

event_spline <- event_data_input %>%
  group_by(Species, Event_ID) %>%
  group_split() %>%
  map_dfr(function(dat) {
    res <- fit_event_spline(dat, spline_par)
    if (is.null(res)) {
      dat %>% arrange(Index) %>%
        mutate(Fitted_log = NA_real_, CI_lwr = NA_real_, CI_upr = NA_real_,
               Spline_spar = NA_real_, Spline_df = NA_real_, Spline_cv = NA_real_,
               Fit_Mode = spline_par$fit_mode, Fit_OK = FALSE)
    } else res
  }) %>%
  ungroup()

## — quality flagging
event_quality <- event_spline %>%
  group_by(Species, Event_ID) %>%
  summarise(
    Global_Event_ID  = first(Global_Event_ID),
    n                = n(),
    Fit_OK           = first(Fit_OK),
    Event_Length     = n(),
    Peak_Index       = first(Peak_Index),
    Event_Start_Index = first(Event_Start_Index),
    Event_End_Index   = first(Event_End_Index),
    Days_Before_Peak  = sum(Index < Peak_Index, na.rm = TRUE),
    Days_After_Peak   = sum(Index > Peak_Index, na.rm = TRUE),
    weak_left_baseline = sum(
      Index < Peak_Index &
      Abundance_log <= spline_par$baseline_prop * max(Abundance_log, na.rm = TRUE),
      na.rm = TRUE) < spline_par$baseline_min_n,
    weak_right_baseline = sum(
      Index > Peak_Index &
      Abundance_log <= spline_par$baseline_prop * max(Abundance_log, na.rm = TRUE),
      na.rm = TRUE) < spline_par$baseline_min_n,
    short_event       = Event_Length < spline_par$min_n,
    short_left_phase  = Days_Before_Peak < 3,
    short_right_phase = Days_After_Peak  < 3,
    max_obs_log   = suppressWarnings(max(Abundance_log, na.rm = TRUE)),
    max_fit_log   = suppressWarnings(max(Fitted_log,    na.rm = TRUE)),
    peak_fit_diff = abs(max_obs_log - max_fit_log),
    poor_peak_fit = ifelse(is.na(peak_fit_diff), TRUE, peak_fit_diff > spline_par$peak_diff_cut),
    Spline_spar = first(Spline_spar),
    Spline_df   = first(Spline_df),
    Spline_cv   = first(Spline_cv),
    .groups = "drop"
  ) %>%
  mutate(
    quality_flag = case_when(
      !Fit_OK                                ~ "fit_failed",
      short_event                            ~ "short_event",
      short_left_phase & short_right_phase   ~ "short_both_phases",
      short_left_phase                       ~ "short_left_phase",
      short_right_phase                      ~ "short_right_phase",
      weak_left_baseline & weak_right_baseline ~ "weak_both_baselines",
      weak_left_baseline                     ~ "weak_left_baseline",
      weak_right_baseline                    ~ "weak_right_baseline",
      poor_peak_fit                          ~ "poor_peak_fit",
      TRUE                                   ~ "good"
    )
  )

event_spline <- event_spline %>%
  left_join(dplyr::select(event_quality, Species, Event_ID, quality_flag),
            by = c("Species", "Event_ID"))

write.xlsx(
  list(spline_fitted  = event_spline,
       spline_summary = arrange(event_quality, Species, Event_ID),
       event_quality  = event_quality),
  file.path(OUT_DIR, "4_event_spline_LOOCV_fitting.xlsx"), overwrite = TRUE
)

p_spline_check <- ggplot(event_spline, aes(x = Date)) +
  geom_ribbon(aes(ymin = CI_lwr, ymax = CI_upr), fill = "grey80", alpha = 0.5) +
  geom_point(aes(y = Abundance_log), color = "blue", pch = 1, size = 1) +
  geom_line(aes(y = Fitted_log), color = "red", linewidth = 0.7) +
  geom_point(
    data = event_spline %>%
      group_by(Species, Event_ID) %>%
      slice_min(abs(Index - Peak_Index), n = 1, with_ties = FALSE) %>%
      ungroup(),
    aes(x = Peak_Date, y = Abundance_log),
    color = "black", size = 1.6, inherit.aes = FALSE
  ) +
  facet_wrap(~Global_Event_ID, scales = "free_y", ncol = 8) +
  labs(x = NULL, y = expression(log[10](x + 1)),
       title = "Step 4. Event-wise cubic smoothing spline (LOOCV / fixed-spar)") +
  theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey95", color = "black"),
        strip.text = element_text(face = "bold", size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank())
print(p_spline_check)
ggsave(file.path(OUT_DIR,"cubic smoothin spline.pdf"), plot = p_spline_check, width = 12, height = 7)

# ── 8) Cardinal phenology variables ──────────────────────────────────────────
phenology_vars <- event_spline %>%
  group_by(Species, Event_ID) %>%
  group_split() %>%
  map_dfr(calc_pheno_one_event)

write.xlsx(phenology_vars,
  file.path(OUT_DIR, "5_phenology_variables.xlsx"), overwrite = TRUE)

## — visualise phenology landmarks
pheno_points <- phenology_vars %>%
  dplyr::select(Species, Event_ID,
                DBS_Date, DMF_Date, DMA_Date, DMM_Date, DBE_Date,
                XO, XF, MA, XM, XE) %>%
  pivot_longer(ends_with("_Date"), names_to = "Variable", values_to = "Date") %>%
  mutate(
    Value = case_when(
      Variable == "DBS_Date" ~ XO,
      Variable == "DMF_Date" ~ XF,
      Variable == "DMA_Date" ~ MA,
      Variable == "DMM_Date" ~ XM,
      Variable == "DBE_Date" ~ XE
    ),
    Variable = gsub("_Date", "", Variable),
    Global_Event_ID = paste0(Species,"_",Event_ID)
  )

p_pheno <- ggplot(event_spline, aes(x = Date, y = Fitted_log)) +
  geom_line(color = "red", linewidth = 0.7) +
  geom_point(aes(y = Abundance_log), color = "blue", size = 0.8) +
  geom_point(data = pheno_points,
             aes(x = Date, y = Value, shape = Variable),
             size = 2, inherit.aes = FALSE) +
  facet_wrap(~Global_Event_ID, scales = "free_y", ncol = 8) +
  labs(x = NULL, y = expression(log[10](x + 1)),
       title = "Step 5. Cardinal phenology variables") +
  theme_bw(base_size = 11)
print(p_pheno)
ggsave(file.path(OUT_DIR,"phenology variables.pdf"), plot = p_pheno, width = 12, height = 7)

message("\n✓ All outputs written to: ", OUT_DIR)
