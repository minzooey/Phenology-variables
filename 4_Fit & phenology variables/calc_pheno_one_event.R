#' Extract cardinal phenology variables for one bloom event
#'
#' Variables
#'   DBS  Day of Bloom Start      (max positive curvature, ascending limb)
#'   DMF  Day of Maximum Fitness  (max first derivative, ascending)
#'   DMA  Day of Maximum Abundance
#'   DMM  Day of Maximum Mortality (min first derivative, descending)
#'   DBE  Day of Bloom End        (min negative curvature, descending limb)
#'   MF   Maximum Fitness rate    (max Δlog-abundance day⁻¹)
#'   MM   Maximum Mortality rate  (min Δlog-abundance day⁻¹)
#'   ONS  Onset duration          DBS → DMF
#'   CLI  Climax duration         DMF → DMA
#'   DEC  Decline duration        DMA → DMM
#'   END  Termination duration    DMM → DBE
#'   IL   Increase length         DBS → DMA
#'   DL   Decrease length         DMA → DBE
#'   BL   Bloom length            DBS → DBE
#'   HAB  Habitat duration        DMF → DMM
#'   SI   Steepness Index (increase)
#'   SD   Steepness Index (decrease)
calc_pheno_one_event <- function(df) {
  df <- df %>% arrange(Index) %>% filter(Fit_OK, !is.na(Fitted_log))
  if (nrow(df) < 5) return(tibble())
  
  y  <- df$Fitted_log
  x  <- df$Index
  dy <- diff(y)
  d2y <- c(NA_real_, diff(dy), NA_real_)
  
  i_dma <- which.max(y)
  
  i_dmf <- if (i_dma > 1) {
    r <- seq_len(i_dma - 1)
    r[which.max(dy[r])]
  } else i_dma
  
  i_dmm <- if (i_dma < length(y)) {
    r <- i_dma:(length(y) - 1L)
    r[which.min(dy[r])] + 1L
  } else i_dma
  
  i_dbs <- if (i_dma > 2) {
    r <- seq_len(i_dma); r[which.max(d2y[r])]
  } else 1L
  
  i_dbe <- if (i_dma < length(y) - 1L) {
    r <- i_dma:length(y); r[which.min(d2y[r])]
  } else length(y)
  
  # enforce temporal order
  i_dbs <- min(i_dbs, i_dmf, i_dma)
  i_dmf <- max(i_dbs, min(i_dmf, i_dma))
  i_dmm <- max(i_dma, min(i_dmm, i_dbe))
  i_dbe <- max(i_dmm, i_dbe)
  
  MA <- y[i_dma]; XO <- y[i_dbs]; XF <- y[i_dmf]
  XM <- y[i_dmm]; XE <- y[i_dbe]
  
  MF <- if (i_dmf <= length(dy))           dy[i_dmf]     else NA_real_
  # FIX B6: MM should be the first-difference at i_dmm (maximum decline rate)
  MM <- if (i_dmm > 1 && (i_dmm - 1) <= length(dy)) dy[i_dmm - 1L] else NA_real_
  
  ONS <- x[i_dmf] - x[i_dbs]; CLI <- x[i_dma] - x[i_dmf]
  DEC <- x[i_dmm] - x[i_dma]; END <- x[i_dbe] - x[i_dmm]
  IL  <- x[i_dma] - x[i_dbs]; DL  <- x[i_dbe] - x[i_dma]
  BL  <- x[i_dbe] - x[i_dbs]; HAB <- x[i_dmm] - x[i_dmf]
  
  SI <- if (IL > 0 && XO > 0) log(MA / XO) / IL else NA_real_
  SD <- if (DL > 0 && XE > 0) log(MA / XE) / DL else NA_real_
  
  tibble(
    Species     = first(df$Species),
    Event_ID    = first(df$Event_ID),
    quality_flag = first(df$quality_flag),
    DBS = x[i_dbs], DMF = x[i_dmf], DMA = x[i_dma], DMM = x[i_dmm], DBE = x[i_dbe],
    DBS_Date = df$Date[i_dbs], DMF_Date = df$Date[i_dmf], DMA_Date = df$Date[i_dma],
    DMM_Date = df$Date[i_dmm], DBE_Date = df$Date[i_dbe],
    XO = XO, XF = XF, MA = MA, XM = XM, XE = XE,
    MF = MF, MM = MM,
    ONS = ONS, CLI = CLI, DEC = DEC, END = END,
    HAB = HAB, IL = IL, DL = DL, BL = BL,
    SI = SI, SD = SD
  )
}