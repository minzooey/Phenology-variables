fit_event_spline <- function(df_event, par = spline_par) {
  df_event <- df_event %>% arrange(Index) %>% filter(!is.na(Abundance_log))
  
  if (nrow(df_event) < par$min_n || length(unique(df_event$Index)) < par$min_n)
    return(NULL)
  
  x <- df_event$Index
  y <- df_event$Abundance_log
  
  fit <- tryCatch({
    if (par$fit_mode == "loocv") {
      smooth.spline(x, y, cv = TRUE,
                    control.spar = list(low = par$spar_min, high = par$spar_max))
    } else if (par$fit_mode == "fixed_spar") {
      smooth.spline(x, y, spar = par$spar_fixed)
    } else stop("par$fit_mode must be 'loocv' or 'fixed_spar'.")
  }, error = function(e) NULL)
  
  if (is.null(fit)) return(NULL)
  
  y_fit <- as.numeric(predict(fit, x)$y)
  resid <- y - y_fit
  spar  <- fit$spar
  
  # Residual bootstrap — vectorised over iterations
  boot_mat <- vapply(seq_len(par$n_boot), function(b) {
    y_b    <- y_fit + sample(resid, replace = TRUE)
    b_fit  <- tryCatch(smooth.spline(x, y_b, spar = spar), error = function(e) NULL)
    if (is.null(b_fit)) rep(NA_real_, length(x))
    else as.numeric(predict(b_fit, x)$y)
  }, numeric(length(x)))
  
  ci_lwr <- apply(boot_mat, 1, quantile, probs = 0.025, na.rm = TRUE)
  ci_upr <- apply(boot_mat, 1, quantile, probs = 0.975, na.rm = TRUE)
  
  df_event %>%
    mutate(
      Fitted_log  = y_fit,
      CI_lwr      = ci_lwr,
      CI_upr      = ci_upr,
      Spline_spar = spar,
      Spline_df   = fit$df,
      Spline_cv   = fit$cv.crit,
      Fit_Mode    = par$fit_mode,
      Fit_OK      = TRUE
    )
}