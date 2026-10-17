#================================================================================#
#
# camsq_mipv_analysis(), version 1.0.0 (2026-10-17)
#
# Model-implied predicted values (MIPVs) for CAMSQ profiles.
#
# Copyright (c) 2026 Leon P. Wendt
# License: CC BY-NC-ND 4.0 (https://creativecommons.org/licenses/by-nc-nd/4.0/)
# Maintained: https://github.com/DrLeonPWendt/camsq-mipv
#
#================================================================================#
camsq_mipv_version <- "1.0.0"

camsq_mipv_analysis <- function(data, 
                                criterion,  # Criterion Variable
                                sc,    # CAMSQ Self-Certainty
                                oc,    # CAMSQ Other-Certainty
                                hdr = 99,
                                min_local_n = 20,
                                # min_local_n: recommended rule: 20 for analyses with at
                                # least 300 complete cases (default), 10 for
                                # smaller analyses. NULL disables the adaptive
                                # HDR selection.
                                local_n_profiles = c("Median", "Hypo", "Hyper", "Effective", "High-Certainty"),
                                hdr_search = seq(99, 1, by = -1),
                                hdr_min = NULL,
                                # hdr_min: NULL (default) applies no lower bound
                                # to the adaptive HDR search (search from 99% down to 1%).
                                local_n_radius = 0.50,
                                # local_n_radius: Mahalanobis distance threshold
                                # of the elliptical local-N neighborhood
                                # (recommended: 0.50).
                                local_n_region_profiles = c("Median", "Hypo", "Hyper", "Effective", "High-Certainty"),
                                local_n_region_alpha = 0.12,
                                local_n_region_linewidth = 0.7,
                                local_n_region_linetype = "solid",
                                plot_profiles = "all",
                                plot_scatter = TRUE,
                                show_local_n_region = TRUE,
                                show_min_max = FALSE,
                                profile_unstd = FALSE,
                                profile_digits = 2,
                                show_r2 = TRUE,
                                r2_position = "topleft",
                                # r2_position: corner of the box with CV-R2 and
                                # adaptive HDR: "topleft" (default, original
                                # placement), "topright", "bottomleft",
                                # "bottomright". show_r2 = FALSE removes the box.
                                show_n = TRUE,
                                # show_n: adds a separate box with the number of
                                # complete cases used in the analysis; independent
                                # of show_r2. show_n = FALSE removes it.
                                n_position = "bottomright",
                                # n_position: corner of the N box, same options as
                                # r2_position (default "bottomright").
                                r2_text_size = 2.7,
                                # r2_text_size: text size of both boxes (ggplot2 units).
                                show_outside = TRUE,
                                model = c("polynomial", "loess"),
                                poly_model = "auto",
                                # poly_model: "auto" (default) selects the model via
                                # 10-fold cross-validated R2 
                                # "linear"/"quadratic"/"cubic"
                                # forces that model without model selection.
                                cv_r2 = TRUE,
                                cv_folds = 10,
                                cv_seed = 1234,
                                loess_span = 0.75,
                                loess_degree = 1,
                                loess_family = c("gaussian", "symmetric"),
                                greyscale = FALSE,
                                grey_power = 1,
                                sc_range = c(1, 7),
                                oc_range = c(1, 7),
                                grid_n = 300,
                                # grid_n: resolution of the prediction grid
                                # over the CAMSQ score space (300 x 300).
                                label_offset = 0.040,
                                label_dir = c(-1, 0),
                                # label_dir: direction of the profile name labels
                                # relative to the circle, as c(x, y) multipliers of
                                # label_offset. The default c(-1, 0) places labels to
                                # the left (original placement). c(1, 1) places them
                                # to the upper right. A named list sets directions
                                # per profile, e.g., list(Effective = c(1, 1)),
                                # unnamed profiles keep the default.
                                profile_label_size = 2.0,
                                title = NULL,
                                scatter_alpha = 0.1,
                                scatter_size = 1.5,
                                auto_install = FALSE,
                                estimator = c("ols", "bayes"),
                                family = gaussian(),
                                prior = NULL,
                                chains = 4,
                                cores = min(4, parallel::detectCores()),
                                iter = 2000,
                                warmup = floor(iter / 2),
                                seed = 1234,
                                backend = "rstan",
                                control = list(adapt_delta = 0.95, max_treedepth = 12),
                                save_pars = brms::save_pars(all = TRUE),
                                stan_threads = FALSE,
                                ...) {
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  # 2. Validate estimator, model selection, and cross-validation settings      ####
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  model <- match.arg(model)
  estimator <- match.arg(estimator)
  loess_family <- match.arg(loess_family)
  
  if (!is.character(poly_model) || length(poly_model) != 1 || is.na(poly_model) ||
      !(poly_model %in% c("linear", "quadratic", "cubic", "auto"))) {
    stop("poly_model must be one of: 'linear', 'quadratic', 'cubic', 'auto'.")
  }
  
  if (!is.logical(cv_r2) || length(cv_r2) != 1 || is.na(cv_r2)) {
    stop("cv_r2 must be either TRUE or FALSE.")
  }
  
  if (!is.logical(show_min_max) || length(show_min_max) != 1 || is.na(show_min_max)) {
    stop("show_min_max must be either TRUE or FALSE.")
  }
  
  if (!is.numeric(loess_span) || length(loess_span) != 1 ||
      is.na(loess_span) || loess_span <= 0) {
    stop("loess_span must be a single positive numeric value.")
  }
  
  if (!is.numeric(loess_degree) || length(loess_degree) != 1 ||
      is.na(loess_degree) || loess_degree != round(loess_degree) ||
      !(loess_degree %in% c(0, 1, 2))) {
    stop("loess_degree must be one of: 0, 1, or 2.")
  }
  loess_degree <- as.integer(loess_degree)
  
  if (identical(model, "polynomial") && identical(poly_model, "auto") && identical(estimator, "bayes")) {
    stop("poly_model = 'auto' selects the model based on cross-validated R², which is currently implemented only for estimator = 'ols'. For Bayesian models, set poly_model explicitly to 'linear', 'quadratic', or 'cubic'.")
  }
  
  if (identical(model, "polynomial") && isTRUE(cv_r2) && identical(estimator, "bayes")) {
    stop("Cross-validated R² is currently implemented only for estimator = 'ols'.")
  }
  
  # Cross-validation is needed only when explicitly requested (cv_r2 = TRUE)
  # or when it is required for automatic model selection (poly_model = 'auto').
  # With an explicit poly_model (linear/quadratic/cubic) and cv_r2 = FALSE,
  # no cross-validation is computed and only the forced model is fitted.
  cv_r2_for_selection <- identical(model, "polynomial") && identical(estimator, "ols") && identical(poly_model, "auto")
  cv_r2_needed <- identical(model, "polynomial") && (isTRUE(cv_r2) || isTRUE(cv_r2_for_selection))
  
  if (isTRUE(cv_r2_needed)) {
    if (!is.numeric(cv_folds) || length(cv_folds) != 1 || is.na(cv_folds) ||
        cv_folds < 2 || cv_folds != round(cv_folds)) {
      stop("cv_folds must be a single integer >= 2.")
    }
    cv_folds <- as.integer(cv_folds)
    
    if (!is.numeric(cv_seed) || length(cv_seed) != 1 || is.na(cv_seed)) {
      stop("cv_seed must be a single numeric value.")
    }
  }
  
  ensure_pkg <- function(p) {
    if (!requireNamespace(p, quietly = TRUE)) {
      if (isTRUE(auto_install)) install.packages(p)
      else stop(sprintf("Package '%s' is required. Set auto_install=TRUE to install automatically.", p))
    }
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  }
  
  pkgs_common <- c("ks", "ggplot2", "dplyr", "tibble", "purrr", "rlang",
                   "scales", "grid", "ggnewscale")
  if (identical(model, "polynomial")) pkgs_common <- c(pkgs_common, "emmeans")
  invisible(lapply(pkgs_common, ensure_pkg))
  if (identical(model, "polynomial") && identical(estimator, "bayes")) ensure_pkg("brms")
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  # 3. Define helper functions                                                ####
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  to_name <- function(arg) if (is.character(arg)) arg else rlang::as_name(rlang::enquo(arg))
  
  coerce_numeric <- function(v) {
    if (is.numeric(v)) return(v)
    if (is.factor(v)) v <- as.character(v)
    suppressWarnings(as.numeric(v))
  }
  
  # Coerce a column to numeric and fail loudly if the coercion produces new
  # missing values, so that non-numeric input cannot silently reduce the
  # complete-case sample.
  coerce_or_stop <- function(v, nm) {
    was_na <- is.na(v)
    out <- coerce_numeric(v)
    if (any(is.na(out) & !was_na)) {
      stop(sprintf("Column '%s' contains non-numeric values that could not be converted to numeric.", nm))
    }
    out
  }
  
  calc_aicc <- function(model) {
    k <- length(stats::coef(model))
    n <- stats::nobs(model)
    aic_val <- stats::AIC(model)
    aic_val + (2 * k * (k + 1)) / (n - k - 1)
  }
  
  # Fold assignment for the 10-fold cross-validation. The user's random seed is saved and
  # restored so the function does not perturb the RNG state of the session.
  make_cv_folds <- function(n, k, seed) {
    if (k > n) stop("cv_folds cannot exceed the number of complete cases.")
    
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      get(".Random.seed", envir = .GlobalEnv)
    } else {
      NULL
    }
    
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    
    set.seed(seed)
    sample(rep(seq_len(k), length.out = n))
  }
  
  # Pooled 10-fold cross-validated R2: within each fold, the model is fitted on
  # the training data and evaluated on the held-out data, with the outcome
  # standardized by the training-fold mean and SD.
  calc_cv_r2_lm <- function(formula, data, outcome, fold_id) {
    k <- length(unique(fold_id))
    sse <- 0
    sst <- 0
    
    for (fold in seq_len(k)) {
      train <- data[fold_id != fold, , drop = FALSE]
      test  <- data[fold_id == fold, , drop = FALSE]
      
      y_train_mean <- mean(train[[outcome]], na.rm = TRUE)
      y_train_sd   <- stats::sd(train[[outcome]], na.rm = TRUE)
      
      if (!is.finite(y_train_mean) || !is.finite(y_train_sd) || y_train_sd <= 0) {
        stop("Cannot compute cross-validated R2 because the outcome has zero or invalid variance in at least one training fold.")
      }
      
      fit_cv <- stats::lm(formula, data = train)
      pred_cv <- as.numeric(stats::predict(fit_cv, newdata = test))
      obs_cv <- (test[[outcome]] - y_train_mean) / y_train_sd
      
      ok <- is.finite(obs_cv) & is.finite(pred_cv)
      sse <- sse + sum((obs_cv[ok] - pred_cv[ok])^2)
      sst <- sst + sum(obs_cv[ok]^2)
    }
    
    if (!is.finite(sst) || sst <= 0) return(NA_real_)
    1 - (sse / sst)
  }
  
  calc_r2_from_pred <- function(obs, pred) {
    ok <- is.finite(obs) & is.finite(pred)
    if (sum(ok) < 2) return(NA_real_)
    
    sse <- sum((obs[ok] - pred[ok])^2)
    sst <- sum((obs[ok] - mean(obs[ok]))^2)
    
    if (!is.finite(sst) || sst <= 0) return(NA_real_)
    1 - (sse / sst)
  }
  
  get_waic <- function(fit) {
    w <- tryCatch(brms::waic(fit), error = function(e) NULL)
    if (is.null(w)) return(NA_real_)
    as.numeric(w$estimates["waic", "Estimate"])
  }
  
  get_bayes_r2 <- function(fit) {
    r2 <- tryCatch(brms::bayes_R2(fit), error = function(e) NULL)
    if (is.null(r2)) return(NA_real_)
    mean(as.numeric(r2))
  }
  
  fmt_p <- function(p) {
    if (is.na(p)) return("NA")
    if (p < .001) return("< .001")
    sprintf("= %.3f", p)
  }
  
  sig_star <- function(p) {
    if (is.na(p)) return("")
    if (p < .001) return("***")
    if (p < .01)  return("**")
    if (p < .05)  return("*")
    ""
  }
  
  cv_txt <- function(value) {
    if (!isTRUE(cv_r2_needed)) return("")
    sprintf(" | CV-R²=%.3f", value)
  }
  
  normalize_profile_names <- function(z) {
    z <- as.character(z)
    z[z == "Certainty"] <- "High-Certainty"
    z
  }
  
  # Squared Mahalanobis distances between query coordinates and one target,
  # based on the empirical covariance of the observed CAMSQ scores (used for
  # the empirical adjustment of profile coordinates).
  calc_mahalanobis_d2 <- function(query_xy, target_xy, S_inv) {
    query_xy <- as.matrix(query_xy)
    target_xy <- as.numeric(target_xy)
    S_inv <- as.matrix(S_inv)
    
    if (ncol(query_xy) != 2 || length(target_xy) != 2) {
      stop("Mahalanobis distance can only be computed for two-dimensional coordinates.")
    }
    
    if (!all(dim(S_inv) == c(2, 2)) || any(!is.finite(S_inv))) {
      stop("S_inv must be a finite 2 x 2 inverse covariance matrix.")
    }
    
    dx <- sweep(query_xy, 2, target_xy, "-")
    rowSums((dx %*% S_inv) * dx)
  }
  
  # Local sample size: the number of observed
  # CAMSQ profiles falling within an elliptical neighborhood around each
  # (adjusted) profile coordinate, defined by a Mahalanobis distance threshold
  # (radius, in covariance-standardized units; 0.50).
  calc_local_n <- function(query_xy, obs_xy, S, radius = 1,
                           chunk_size = 1000) {
    query_xy <- as.matrix(query_xy)
    obs_xy <- as.matrix(obs_xy)
    S <- as.matrix(S)
    
    if (ncol(query_xy) != 2 || ncol(obs_xy) != 2) {
      stop("Local sample size can only be computed for two-dimensional coordinates.")
    }
    
    if (!all(dim(S) == c(2, 2)) || any(!is.finite(S))) {
      stop("S must be a finite 2 x 2 covariance matrix.")
    }
    
    if (!is.numeric(radius) || length(radius) != 1 ||
        is.na(radius) || radius <= 0) {
      stop("radius must be a single positive numeric value.")
    }
    
    S_inv <- tryCatch(solve(S), error = function(e) NULL)
    if (is.null(S_inv) || any(!is.finite(S_inv))) {
      stop("Could not invert the empirical covariance matrix for Mahalanobis local sample size estimation.")
    }
    
    radius2 <- radius^2
    out <- numeric(nrow(query_xy))
    starts <- seq(1, nrow(query_xy), by = chunk_size)
    
    for (s in starts) {
      e <- min(s + chunk_size - 1, nrow(query_xy))
      block <- query_xy[s:e, , drop = FALSE]
      
      dx <- outer(block[, 1], obs_xy[, 1], "-")
      dy <- outer(block[, 2], obs_xy[, 2], "-")
      
      d2 <- S_inv[1, 1] * dx^2 +
        (S_inv[1, 2] + S_inv[2, 1]) * dx * dy +
        S_inv[2, 2] * dy^2
      
      out[s:e] <- rowSums(d2 <= radius2)
    }
    
    out
  }
  
  # Outline of the elliptical local-N neighborhood, for optional display in
  # the surface plot.
  calc_local_n_region <- function(points_df, S, radius = 1, n = 160) {
    if (is.null(points_df) || nrow(points_df) < 1) {
      return(tibble::tibble())
    }
    
    S <- as.matrix(S)
    
    if (!all(dim(S) == c(2, 2)) || any(!is.finite(S))) {
      stop("S must be a finite 2 x 2 covariance matrix.")
    }
    
    if (!is.numeric(radius) || length(radius) != 1 ||
        is.na(radius) || radius <= 0) {
      stop("radius must be a single positive numeric value.")
    }
    
    R <- tryCatch(chol(S), error = function(e) NULL)
    if (is.null(R) || any(!is.finite(R))) {
      stop("Could not decompose the empirical covariance matrix for Mahalanobis local sample size region plotting.")
    }
    
    theta <- seq(0, 2 * pi, length.out = n)
    unit_circle <- cbind(cos(theta), sin(theta))
    ellipse_offsets <- radius * unit_circle %*% R
    
    purrr::pmap_dfr(
      list(points_df$x, points_df$y, points_df$label),
      function(cx, cy, lab) {
        tibble::tibble(
          x = cx + ellipse_offsets[, 1],
          y = cy + ellipse_offsets[, 2],
          label = lab,
          local_n_radius = radius,
          local_n_metric = "mahalanobis"
        )
      }
    )
  }
  
  get_hdr_cutoff <- function(kde_fit, cont_value) {
    cutoff_value <- tryCatch(
      ks::contourLevels(kde_fit, cont = cont_value),
      error = function(e) NA_real_
    )
    if (length(cutoff_value) != 1 || is.na(cutoff_value)) {
      return(NA_real_)
    }
    as.numeric(cutoff_value)
  }
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  # 4. Validate plotting and contour settings                                 ####
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  
  if (is.logical(plot_profiles)) plot_profiles <- if (isTRUE(plot_profiles)) "all" else FALSE
  if (!(identical(plot_profiles, FALSE) || (is.character(plot_profiles) && length(plot_profiles) == 1 &&
                                            plot_profiles %in% c("constructs", "all")))) {
    stop("plot_profiles must be FALSE, 'constructs', or 'all'.")
  }
  
  if (!is.character(r2_position) || length(r2_position) != 1 || is.na(r2_position) ||
      !(r2_position %in% c("topleft", "topright", "bottomleft", "bottomright"))) {
    stop("r2_position must be one of: 'topleft', 'topright', 'bottomleft', 'bottomright'.")
  }
  
  if (!is.logical(show_n) || length(show_n) != 1 || is.na(show_n)) {
    stop("show_n must be either TRUE or FALSE.")
  }
  
  if (!is.character(n_position) || length(n_position) != 1 || is.na(n_position) ||
      !(n_position %in% c("topleft", "topright", "bottomleft", "bottomright"))) {
    stop("n_position must be one of: 'topleft', 'topright', 'bottomleft', 'bottomright'.")
  }
  
  if (!is.numeric(r2_text_size) || length(r2_text_size) != 1 || is.na(r2_text_size) || r2_text_size <= 0) {
    stop("r2_text_size must be a single positive numeric value.")
  }
  
  hdr_is_auto <- is.character(hdr) && length(hdr) == 1 && identical(tolower(hdr), "auto")
  adaptive_hdr <- hdr_is_auto || !is.null(min_local_n)
  
  if (!hdr_is_auto) {
    if (!is.numeric(hdr) || length(hdr) != 1 || is.na(hdr) || hdr <= 0 || hdr >= 100) {
      stop("hdr must be a single numeric value between 0 and 100, or 'auto'.")
    }
    hdr <- as.numeric(hdr)
  }
  
  if (!is.null(min_local_n)) {
    if (!is.numeric(min_local_n) || length(min_local_n) != 1 || is.na(min_local_n) || min_local_n < 0) {
      stop("min_local_n must be NULL or a single non-negative numeric value.")
    }
    min_local_n <- as.numeric(min_local_n)
  }
  
  if (hdr_is_auto && is.null(min_local_n)) {
    stop("When hdr = 'auto', min_local_n must be set, for example min_local_n = 20.")
  }
  
  if (!is.numeric(hdr_search) || length(hdr_search) < 1 || anyNA(hdr_search) ||
      any(hdr_search <= 0) || any(hdr_search >= 100)) {
    stop("hdr_search must be a numeric vector with values between 0 and 100.")
  }
  hdr_search <- sort(unique(as.numeric(hdr_search)), decreasing = TRUE)
  
  # Optional lower bound for the adaptive HDR search. With hdr_min = NULL
  # (default), no lower bound is applied and the search may proceed down to
  # the smallest value in hdr_search.
  if (!is.null(hdr_min)) {
    if (!is.numeric(hdr_min) || length(hdr_min) != 1 || is.na(hdr_min) ||
        hdr_min < 1 || hdr_min > 99) {
      stop("hdr_min must be NULL or a single numeric value between 1 and 99.")
    }
    hdr_min <- as.numeric(hdr_min)
    
    if (!hdr_is_auto && hdr < hdr_min) {
      stop("hdr must not be smaller than hdr_min.")
    }
    
    hdr_search <- hdr_search[hdr_search >= hdr_min]
    
    if (length(hdr_search) < 1) {
      stop("After applying hdr_min, no candidate HDR values remain in hdr_search. Lower hdr_min or adjust hdr_search.")
    }
  }
  
  if (!is.numeric(local_n_radius) || length(local_n_radius) != 1 ||
      is.na(local_n_radius) || local_n_radius <= 0) {
    stop("local_n_radius must be a single positive numeric value.")
  }
  local_n_radius <- as.numeric(local_n_radius)
  
  valid_profile_names <- c("Median", "Hypo", "Hyper", "Effective", "High-Certainty")
  if (isTRUE(show_min_max)) {
    valid_profile_names <- c(valid_profile_names, "Max", "Min")
  }
  valid_profile_names_txt <- paste(valid_profile_names, collapse = ", ")
  
  local_n_profiles <- normalize_profile_names(local_n_profiles)
  
  if (!is.character(local_n_profiles) || length(local_n_profiles) < 1 ||
      anyNA(local_n_profiles) || !all(local_n_profiles %in% valid_profile_names)) {
    stop(sprintf(
      "local_n_profiles must contain one or more of: %s.",
      valid_profile_names_txt
    ))
  }
  local_n_profiles <- unique(local_n_profiles)
  
  if (is.null(local_n_region_profiles)) {
    local_n_region_profiles <- local_n_profiles
  }
  
  local_n_region_profiles <- normalize_profile_names(local_n_region_profiles)
  
  if (!is.logical(show_local_n_region) || length(show_local_n_region) != 1 || is.na(show_local_n_region)) {
    stop("show_local_n_region must be either TRUE or FALSE.")
  }
  
  if (!is.character(local_n_region_profiles) || length(local_n_region_profiles) < 1 ||
      anyNA(local_n_region_profiles) || !all(local_n_region_profiles %in% valid_profile_names)) {
    stop(sprintf(
      "local_n_region_profiles must contain one or more of: %s.",
      valid_profile_names_txt
    ))
  }
  local_n_region_profiles <- unique(local_n_region_profiles)
  
  if (!is.numeric(local_n_region_alpha) || length(local_n_region_alpha) != 1 ||
      is.na(local_n_region_alpha) || local_n_region_alpha < 0 || local_n_region_alpha > 1) {
    stop("local_n_region_alpha must be a single numeric value between 0 and 1.")
  }
  
  if (!is.numeric(local_n_region_linewidth) || length(local_n_region_linewidth) != 1 ||
      is.na(local_n_region_linewidth) || local_n_region_linewidth < 0) {
    stop("local_n_region_linewidth must be a single non-negative numeric value.")
  }
  
  if (!is.character(local_n_region_linetype) || length(local_n_region_linetype) != 1 ||
      is.na(local_n_region_linetype)) {
    stop("local_n_region_linetype must be a single character value.")
  }
  
  # Label direction: a numeric c(x, y) vector applies to all profiles; a named
  # list of such vectors sets directions per profile, others keep the default
  # c(-1, 0) (label to the left of the circle, original placement).
  label_dir_default <- c(-1, 0)
  if (is.numeric(label_dir) && !is.list(label_dir)) {
    if (length(label_dir) != 2 || anyNA(label_dir)) {
      stop("label_dir must be a numeric vector of length 2, e.g., c(-1, 0) or c(1, 1), or a named list of such vectors.")
    }
    label_dir_default <- as.numeric(label_dir)
    label_dir <- list()
  } else if (is.list(label_dir)) {
    if (length(label_dir) > 0) {
      if (is.null(names(label_dir)) || any(!nzchar(names(label_dir)))) {
        stop("When label_dir is a list, its elements must be named by profile.")
      }
      names(label_dir) <- normalize_profile_names(names(label_dir))
      if (!all(names(label_dir) %in% valid_profile_names)) {
        stop(sprintf("label_dir names must be among: %s.", valid_profile_names_txt))
      }
      if (!all(vapply(label_dir, function(v) is.numeric(v) && length(v) == 2 && !anyNA(v), logical(1)))) {
        stop("Each element of label_dir must be a numeric vector of length 2.")
      }
    }
  } else {
    stop("label_dir must be a numeric vector of length 2 or a named list of such vectors.")
  }
  
  criterion_nm <- to_name(criterion)
  sc_nm   <- to_name(sc)
  oc_nm   <- to_name(oc)
  
  vars <- unique(c(criterion_nm, sc_nm, oc_nm))
  df_all_vars <- data[, vars, drop = FALSE]
  
  # All analysis columns are coerced to numeric; non-numeric input fails
  # loudly instead of silently shrinking the complete-case sample.
  df_all_vars[[criterion_nm]] <- coerce_or_stop(df_all_vars[[criterion_nm]], criterion_nm)
  df_all_vars[[sc_nm]]   <- coerce_or_stop(df_all_vars[[sc_nm]], sc_nm)
  df_all_vars[[oc_nm]]   <- coerce_or_stop(df_all_vars[[oc_nm]], oc_nm)
  
  if (!is.numeric(profile_digits) || length(profile_digits) != 1 || is.na(profile_digits) ||
      profile_digits < 0 || profile_digits != round(profile_digits)) {
    stop("profile_digits must be a single non-negative integer.")
  }
  profile_digits <- as.integer(profile_digits)
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  # 5. Prepare complete-case data and estimate kernel density                 ####
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  
  # Each analysis is based on
  # participants with complete data on the CAMSQ scales and the respective
  # criterion variable. This single case base enters the kernel density
  # estimation, the RSA model, the Mahalanobis metric, and the local sample
  # size counts.
  df <- df_all_vars[stats::complete.cases(df_all_vars), , drop = FALSE]
  n_used <- nrow(df)
  
  if (n_used < 3) {
    stop("At least three complete cases are required for kernel density estimation, RSA, and local sample size estimation.")
  }
  
  observed_points <- tibble::tibble(
    x = df[[sc_nm]],
    y = df[[oc_nm]]
  )
  observed_matrix_kde <- as.matrix(observed_points)
  
  # Empirical covariance of the observed CAMSQ scores among complete cases;
  # defines the Mahalanobis metric for the profile adjustment and the
  # elliptical local-N neighborhoods. Identical to the KDE case base (single
  # case base per criterion).
  observed_matrix_local <- observed_matrix_kde
  
  S_local <- stats::cov(observed_matrix_local, use = "complete.obs")
  S_local_eigen <- tryCatch(
    eigen(S_local, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) NA_real_
  )
  
  if (!all(dim(S_local) == c(2, 2)) ||
      any(!is.finite(S_local)) ||
      any(!is.finite(S_local_eigen)) ||
      any(S_local_eigen <= 0)) {
    stop("Could not compute a positive definite empirical covariance matrix for Mahalanobis local sample size estimation based on complete cases.")
  }
  
  # Kernel density estimation of the observed CAMSQ score distribution, using a plug-in
  # bandwidth (ks::Hpi) and evaluated on the CAMSQ score range (1 to 7 on
  # both scales).
  H <- ks::Hpi(observed_matrix_kde)
  
  kde_fit <- ks::kde(
    x = observed_matrix_kde,
    H = H,
    xmin = c(1, 1),
    xmax = c(7, 7),
    compute.cont = TRUE
  )
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  # 6. Prepare model formulas                                                 ####
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  
  if (isTRUE(cv_r2_needed) && cv_folds > n_used) {
    stop("cv_folds cannot exceed the number of complete cases used for RSA.")
  }
  
  cv_fold_id <- if (isTRUE(cv_r2_needed)) make_cv_folds(n_used, cv_folds, cv_seed) else NULL
  
  # Mean and SD of the criterion, used to convert standardized predictions
  # back to the original metric (e.g., age in years via mipv_original).
  criterion_mean <- mean(df[[criterion_nm]], na.rm = TRUE)
  criterion_sd   <- stats::sd(df[[criterion_nm]], na.rm = TRUE)
  if (!is.finite(criterion_mean) || !is.finite(criterion_sd) || criterion_sd <= 0) {
    criterion_mean <- NA_real_
    criterion_sd   <- NA_real_
  }
  
  # Criterion variables are
  # z-standardized (scale() on the left-hand side), and three candidate RSA
  # models of increasing complexity are considered:
  # linear (SC + OC), quadratic (+ SC^2, SC*OC, OC^2), and cubic (+ SC^3,
  # SC^2*OC, SC*OC^2, OC^3).
  if (identical(model, "loess")) {
    df$.y_z <- as.numeric(scale(df[[criterion_nm]]))
    fml_loess <- stats::as.formula(paste0(".y_z ~ ", sc_nm, " + ", oc_nm))
  } else if (identical(estimator, "ols")) {
    fml_lin <- stats::as.formula(paste0("scale(", criterion_nm, ") ~ ", sc_nm, " + ", oc_nm))
    fml_quad_explicit <- stats::as.formula(
      paste0("scale(", criterion_nm, ") ~ ", sc_nm, " + ", oc_nm, " + ",
             "I(", sc_nm, "^2) + I(", sc_nm, "*", oc_nm, ") + I(", oc_nm, "^2)")
    )
    fml_cub_explicit <- stats::as.formula(
      paste0("scale(", criterion_nm, ") ~ ", sc_nm, " + ", oc_nm, " + ",
             "I(", sc_nm, "^2) + I(", sc_nm, "*", oc_nm, ") + I(", oc_nm, "^2) + ",
             "I(", sc_nm, "^3) + I(", sc_nm, "^2*", oc_nm, ") + I(", sc_nm, "*", oc_nm, "^2) + I(", oc_nm, "^3)")
    )
  } else {
    df$.y_z <- as.numeric(scale(df[[criterion_nm]]))
    fml_lin <- stats::as.formula(paste0(".y_z ~ ", sc_nm, " + ", oc_nm))
    fml_quad_explicit <- stats::as.formula(
      paste0(".y_z ~ ", sc_nm, " + ", oc_nm, " + ",
             "I(", sc_nm, "^2) + I(", sc_nm, "*", oc_nm, ") + I(", oc_nm, "^2)")
    )
    fml_cub_explicit <- stats::as.formula(
      paste0(".y_z ~ ", sc_nm, " + ", oc_nm, " + ",
             "I(", sc_nm, "^2) + I(", sc_nm, "*", oc_nm, ") + I(", oc_nm, "^2) + ",
             "I(", sc_nm, "^3) + I(", sc_nm, "^2*", oc_nm, ") + I(", sc_nm, "*", oc_nm, "^2) + I(", oc_nm, "^3)")
    )
  }
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  # 7. Define model-fitting functions                                         ####
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  fit_lm <- function(formula, data) stats::lm(formula, data = data)
  
  fit_loess <- function(formula, data) {
    stats::loess(
      formula = formula,
      data = data,
      span = loess_span,
      degree = loess_degree,
      family = loess_family,
      control = stats::loess.control(surface = "direct")
    )
  }
  
  fit_brms <- function(formula, data) {
    brms::brm(
      formula = formula,
      data = data,
      family = family,
      prior = prior,
      chains = chains,
      cores = cores,
      iter = iter,
      warmup = warmup,
      seed = seed,
      backend = backend,
      control = control,
      save_pars = save_pars,
      threads = if (isTRUE(stan_threads)) brms::threading(threads = cores) else NULL,
      refresh = 0,
      ...
    )
  }
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  # 8. Fit models and select model                                            ####
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  selected_model_name <- "User-Defined"
  
  r2_tbl <- NULL
  p_lin_quad <- NA_real_
  p_quad_cub <- NA_real_
  comp_tbl <- NULL
  
  if (identical(model, "loess")) {
    
    # Local regression surface.
    fml <- fml_loess
    fit <- fit_loess(fml, data = df)
    
    pred_train <- as.numeric(stats::predict(fit, newdata = df))
    loess_r2 <- calc_r2_from_pred(df$.y_z, pred_train)
    
    r2_tbl <- tibble::tibble(
      model = "LOESS",
      r2 = loess_r2,
      cv_r2 = NA_real_,
      adj_r2 = NA_real_,
      aic = NA_real_,
      bic = NA_real_,
      aicc = NA_real_,
      span = loess_span,
      degree = loess_degree,
      family = loess_family
    )
    
    selected_model_name <- "Local Regression: LOESS"
    
  } else if (identical(poly_model, "auto")) {
    
    # The three candidate models
    # are fitted, and the model with the highest 10-fold cross-validated R2
    # is selected for computing the MIPVs. 
    fit_lin  <- fit_lm(fml_lin, data = df)
    fit_quad <- fit_lm(fml_quad_explicit, data = df)
    fit_cub  <- fit_lm(fml_cub_explicit, data = df)
    
    cv_vals <- c(
      calc_cv_r2_lm(fml_lin, data = df, outcome = criterion_nm, fold_id = cv_fold_id),
      calc_cv_r2_lm(fml_quad_explicit, data = df, outcome = criterion_nm, fold_id = cv_fold_id),
      calc_cv_r2_lm(fml_cub_explicit, data = df, outcome = criterion_nm, fold_id = cv_fold_id)
    )
    
    r2_tbl <- tibble::tibble(
      model = c("Linear", "Quadratic", "Cubic"),
      r2 = c(summary(fit_lin)$r.squared,
             summary(fit_quad)$r.squared,
             summary(fit_cub)$r.squared),
      cv_r2 = cv_vals,
      adj_r2 = c(summary(fit_lin)$adj.r.squared,
                 summary(fit_quad)$adj.r.squared,
                 summary(fit_cub)$adj.r.squared),
      aic  = c(stats::AIC(fit_lin), stats::AIC(fit_quad), stats::AIC(fit_cub)),
      bic  = c(stats::BIC(fit_lin), stats::BIC(fit_quad), stats::BIC(fit_cub)),
      aicc = c(calc_aicc(fit_lin), calc_aicc(fit_quad), calc_aicc(fit_cub))
    )
    
    comp_lin_quad <- stats::anova(fit_lin, fit_quad)
    comp_quad_cub <- stats::anova(fit_quad, fit_cub)
    
    get_p <- function(anv) {
      if (is.null(anv) || nrow(anv) < 2) return(NA_real_)
      pcol <- which(names(anv) %in% c("Pr(>F)", "Pr(>Chisq)"))
      if (length(pcol) == 0) return(NA_real_)
      as.numeric(anv[nrow(anv), pcol[1]])
    }
    
    p_lin_quad <- get_p(comp_lin_quad)
    p_quad_cub <- get_p(comp_quad_cub)
    
    if (all(!is.finite(r2_tbl$cv_r2))) {
      stop("Automatic model selection failed because all cross-validated R² values are invalid.")
    }
    
    # which.max() returns the first maximum, so exact ties in CV-R2 favor the
    # simpler model (Linear before Quadratic before Cubic).
    best_idx <- which.max(ifelse(is.finite(r2_tbl$cv_r2), r2_tbl$cv_r2, -Inf))
    best_model <- r2_tbl$model[best_idx]
    
    if (best_model == "Cubic") {
      fml <- fml_cub_explicit
      fit <- fit_cub
    } else if (best_model == "Quadratic") {
      fml <- fml_quad_explicit
      fit <- fit_quad
    } else {
      fml <- fml_lin
      fit <- fit_lin
    }
    
    selected_model_name <- paste0("CV-R²-Selected: ", best_model)
    
  } else {
    
    # Forced model: poly_model = 'linear', 'quadratic', or 'cubic'.
    # Only the requested model is fitted. No model comparison is performed,
    # and cross-validated R² is computed only if cv_r2 = TRUE. 
    if (identical(poly_model, "cubic")) {
      fml <- fml_cub_explicit
      selected_model_name <- "User-Defined: Cubic"
    } else if (identical(poly_model, "quadratic")) {
      fml <- fml_quad_explicit
      selected_model_name <- "User-Defined: Quadratic"
    } else {
      fml <- fml_lin
      selected_model_name <- "User-Defined: Linear"
    }
    
    if (identical(estimator, "ols")) {
      fit <- fit_lm(fml, data = df)
      
      cv_val <- if (isTRUE(cv_r2_needed)) {
        calc_cv_r2_lm(fml, data = df, outcome = criterion_nm, fold_id = cv_fold_id)
      } else {
        NA_real_
      }
      
      r2_tbl <- tibble::tibble(
        model = sub("^User-Defined: ", "", selected_model_name),
        r2 = summary(fit)$r.squared,
        cv_r2 = cv_val,
        adj_r2 = summary(fit)$adj.r.squared,
        aic  = stats::AIC(fit),
        bic  = stats::BIC(fit),
        aicc = calc_aicc(fit)
      )
    } else {
      # Bayesian estimation (brms).
      fit <- fit_brms(fml, data = df)
      comp_tbl <- tibble::tibble(
        model = sub("^User-Defined: ", "", selected_model_name),
        bayes_r2 = get_bayes_r2(fit),
        waic = get_waic(fit)
      )
    }
  }
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  # 9. Create prediction grid and profile targets                             ####
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  
  if (!is.numeric(grid_n) || length(grid_n) != 1 || is.na(grid_n) || grid_n < 50) {
    stop("grid_n must be a single numeric value >= 50.")
  }
  grid_n <- as.integer(round(grid_n))
  
  # Prediction grid over the full CAMSQ score space (default 300 x 300).
  # Empirically adjusted profile coordinates are selected from this grid.
  grid <- base::expand.grid(
    x = seq(1, 7, length.out = grid_n),
    y = seq(1, 7, length.out = grid_n)
  )
  names(grid) <- c(sc_nm, oc_nm)
  
  if (identical(model, "loess")) {
    grid$predicted <- as.numeric(stats::predict(fit, newdata = grid))
  } else if (identical(estimator, "ols")) {
    grid$predicted <- stats::predict(fit, newdata = grid)
  } else {
    ep <- brms::posterior_epred(fit, newdata = grid)
    grid$predicted <- as.numeric(colMeans(ep))
  }
  
  grid_can <- grid %>% dplyr::rename(x = !!sc_nm, y = !!oc_nm)
  grid_can$density <- as.numeric(stats::predict(kde_fit, x = as.matrix(grid_can[, c("x", "y")])))
  
  S_profile_inv <- tryCatch(solve(S_local), error = function(e) NULL)
  
  if (is.null(S_profile_inv) || any(!is.finite(S_profile_inv))) {
    stop("Could not invert the empirical covariance matrix for Mahalanobis profile adjustment.")
  }
  
  find_closest_in_region <- function(tx, ty, region_df, S_inv) {
    if (nrow(region_df) < 1) {
      stop("The selected HDR does not contain any grid points. Try a larger hdr value or a finer grid_n.")
    }
    
    region_xy <- as.matrix(region_df[, c("x", "y")])
    target_xy <- c(tx, ty)
    
    d2 <- calc_mahalanobis_d2(
      query_xy = region_xy,
      target_xy = target_xy,
      S_inv = S_inv
    )
    
    idx <- which.min(d2)
    region_df[idx, c("x", "y"), drop = FALSE]
  }
  
  # Theoretical profile coordinates in the SC-OC space.
  theoretical_profiles <- tibble::tribble(
    ~x_theoretical, ~y_theoretical, ~label,
    5.47,           4.70,           "Median",
    1.00,           1.00,           "Hypo",
    1.00,           7.00,           "Hyper",
    7.00,           4.70,           "Effective",
    7.00,           7.00,           "High-Certainty"
  )
  
  make_profile_points <- function(cutoff_value) {
    region_df <- grid_can %>% dplyr::filter(density >= cutoff_value)
    
    if (nrow(region_df) < 1) {
      return(NULL)
    }
    
    region_pred <- region_df %>%
      dplyr::filter(is.finite(predicted))
    
    if (nrow(region_pred) < 1) {
      stop("No finite predicted values were available within the selected HDR. Check the fitted model and prediction grid.")
    }
    
    # Empirically adjusted profile coordinates:
    # for each theoretical profile, the coordinate within the HDR closest to
    # the theoretical location, using Mahalanobis distance.
    empirical_targets <- theoretical_profiles %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        closest = list(find_closest_in_region(
          tx = x_theoretical,
          ty = y_theoretical,
          region_df = region_df,
          S_inv = S_profile_inv
        )),
        x = closest$x,
        y = closest$y
      ) %>%
      dplyr::ungroup() %>%
      dplyr::select(x, y, label, x_theoretical, y_theoretical)
    
    # Optional empirical extrema of the predicted surface within the HDR.
    empirical_extrema <- tibble::tibble()
    
    if (isTRUE(show_min_max)) {
      max_target <- region_pred %>%
        dplyr::slice_max(order_by = predicted, n = 1, with_ties = FALSE) %>%
        dplyr::transmute(
          x = x,
          y = y,
          label = "Max",
          x_theoretical = NA_real_,
          y_theoretical = NA_real_
        )
      
      min_target <- region_pred %>%
        dplyr::slice_min(order_by = predicted, n = 1, with_ties = FALSE) %>%
        dplyr::transmute(
          x = x,
          y = y,
          label = "Min",
          x_theoretical = NA_real_,
          y_theoretical = NA_real_
        )
      
      empirical_extrema <- dplyr::bind_rows(max_target, min_target)
    }
    
    # All profile coordinates are selected from region_df and therefore lie
    # inside the HDR by construction.
    profile_pts_tmp <- dplyr::bind_rows(empirical_targets, empirical_extrema)
    
    # Local sample size at the adjusted coordinates.
    profile_pts_tmp$local_n <- calc_local_n(
      query_xy = as.matrix(profile_pts_tmp[, c("x", "y")]),
      obs_xy = observed_matrix_local,
      S = S_local,
      radius = local_n_radius
    )
    
    list(points = profile_pts_tmp, region = region_df)
  }
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  # 10. Select highest-density region and local-N regions                     ####
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  
  # beginning with the 99%-HDR, increasingly
  # restrictive HDRs are examined in steps of 1% (down to the 1%-HDR, no
  # lower bound by default), and the largest HDR is selected in which all
  # tested profiles satisfy the minimum local sample size criterion. If no
  # HDR satisfies the criterion, the analysis is not conducted (stop below).
  hdr_diagnostics <- NULL
  
  if (isTRUE(adaptive_hdr)) {
    if (hdr_is_auto) {
      candidate_hdr <- hdr_search
    } else {
      candidate_hdr <- hdr_search[hdr_search <= hdr]
      candidate_hdr <- sort(unique(c(hdr, candidate_hdr)), decreasing = TRUE)
    }
    
    if (length(candidate_hdr) < 1) {
      stop("No candidate HDR values are available. Check hdr and hdr_search.")
    }
    
    hdr_diagnostics_list <- vector("list", length(candidate_hdr))
    selected <- NULL
    
    for (i in seq_along(candidate_hdr)) {
      hdr_i <- candidate_hdr[i]
      cutoff_i <- get_hdr_cutoff(kde_fit, hdr_i)
      
      if (!is.finite(cutoff_i)) {
        hdr_diagnostics_list[[i]] <- tibble::tibble(
          level = hdr_i,
          hdr_cutoff = NA_real_,
          min_local_n = NA_real_,
          passed = FALSE
        )
        next
      }
      
      tmp <- make_profile_points(cutoff_i)
      
      if (is.null(tmp)) {
        hdr_diagnostics_list[[i]] <- tibble::tibble(
          level = hdr_i,
          hdr_cutoff = cutoff_i,
          min_local_n = NA_real_,
          passed = FALSE
        )
        next
      }
      
      local_support_i <- tmp$points %>%
        dplyr::filter(label %in% local_n_profiles)
      
      min_local_i <- min(local_support_i$local_n, na.rm = TRUE)
      passed_i <- is.finite(min_local_i) && min_local_i >= min_local_n
      
      hdr_diagnostics_list[[i]] <- tibble::tibble(
        level = hdr_i,
        hdr_cutoff = cutoff_i,
        min_local_n = min_local_i,
        passed = passed_i
      )
      
      if (isTRUE(passed_i)) {
        selected <- list(
          level = hdr_i,
          hdr_cutoff = cutoff_i,
          points = tmp$points,
          region = tmp$region
        )
        break
      }
    }
    
    hdr_diagnostics <- dplyr::bind_rows(hdr_diagnostics_list)
    
    # Infeasibility outcome.
    if (is.null(selected)) {
      best_observed <- hdr_diagnostics %>%
        dplyr::filter(is.finite(min_local_n)) %>%
        dplyr::arrange(dplyr::desc(min_local_n)) %>%
        dplyr::slice(1)
      
      if (nrow(best_observed) == 1) {
        stop(sprintf(
          "No candidate HDR satisfied min_local_n = %.0f for profiles: %s. Best observed minimum local N was %.0f at the %.1f%%-HDR. Reduce min_local_n, expand hdr_search, lower or remove hdr_min, or increase local_n_radius.",
          min_local_n,
          paste(local_n_profiles, collapse = ", "),
          best_observed$min_local_n,
          best_observed$level
        ))
      } else {
        stop(sprintf(
          "No candidate HDR satisfied min_local_n = %.0f. Reduce min_local_n, expand hdr_search, lower or remove hdr_min, or increase local_n_radius.",
          min_local_n
        ))
      }
    }
    
    selected_hdr <- selected$level
    hdr_cutoff <- selected$hdr_cutoff
    profile_pts <- selected$points
    
  } else {
    # Fixed (non-adaptive) HDR (min_local_n = NULL).
    selected_hdr <- hdr
    hdr_cutoff <- get_hdr_cutoff(kde_fit, selected_hdr)
    if (length(hdr_cutoff) != 1 || is.na(hdr_cutoff)) stop("Could not compute a single HDR cutoff.")
    
    tmp <- make_profile_points(hdr_cutoff)
    if (is.null(tmp)) {
      stop("The selected HDR does not contain any grid points. Try a larger hdr value or a finer grid_n.")
    }
    
    profile_pts <- tmp$points
  }
  
  grid_can$in_cont <- grid_can$density >= hdr_cutoff
  
  local_support <- profile_pts %>%
    dplyr::mutate(local_n = as.integer(local_n)) %>%
    dplyr::select(label, x, y, local_n)
  
  local_n_region_df <- tibble::tibble()
  
  if (isTRUE(show_local_n_region)) {
    local_n_region_df <- profile_pts %>%
      dplyr::filter(label %in% local_n_region_profiles) %>%
      calc_local_n_region(
        S = S_local,
        radius = local_n_radius
      )
  }
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  # 11. Compute profile predictions and plotting values                       ####
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  
  # Model-implied predicted values (MIPVs) at the empirically adjusted
  # profile coordinates, on the standardized criterion metric (z) and
  # back-transformed to the original metric (pred_orig; used, e.g., to
  # report age in years).
  newdata_pts <- data.frame(profile_pts$x, profile_pts$y)
  names(newdata_pts) <- c(sc_nm, oc_nm)
  
  if (identical(model, "loess")) {
    profile_pts$pred <- as.numeric(stats::predict(fit, newdata = newdata_pts))
  } else if (identical(estimator, "ols")) {
    profile_pts$pred <- stats::predict(fit, newdata = newdata_pts)
  } else {
    ep_pts <- brms::posterior_epred(fit, newdata = newdata_pts)
    profile_pts$pred <- as.numeric(colMeans(ep_pts))
  }
  
  profile_pts$pred_orig <- if (is.finite(criterion_mean) && is.finite(criterion_sd)) profile_pts$pred * criterion_sd + criterion_mean else NA_real_
  
  fmt_unstd <- paste0("%.", profile_digits, "f")
  fmt_std   <- paste0("%+.", profile_digits, "f")
  profile_pts$pred_label <- if (isTRUE(profile_unstd) && all(is.finite(profile_pts$pred_orig))) {
    sprintf(fmt_unstd, profile_pts$pred_orig)
  } else {
    sprintf(fmt_std, profile_pts$pred)
  }
  
  z_lim    <- c(-2.1, 2.1)
  z_breaks <- seq(-2, 2, by = 1)
  
  if (isTRUE(greyscale)) {
    transform_fun_fill <- function(z) sign(z) * abs(z)^grey_power
    grid_can$fill_val <- transform_fun_fill(grid_can$predicted)
    profile_pts$fill_val      <- transform_fun_fill(profile_pts$pred)
  } else {
    grid_can$fill_val <- grid_can$predicted
    profile_pts$fill_val      <- profile_pts$pred
  }
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  # 12. Prepare plotting data and color scales                                ####
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  
  x_limits <- if (!is.null(sc_range)) sc_range else c(1, 7)
  y_limits <- if (!is.null(oc_range)) oc_range else c(1, 7)
  
  grid_plot <- grid_can %>%
    dplyr::filter(x >= x_limits[1], x <= x_limits[2],
                  y >= y_limits[1], y <= y_limits[2])
  
  df_in  <- grid_plot %>% dplyr::filter(in_cont)
  df_out <- grid_plot %>% dplyr::filter(!in_cont)
  
  # Local-N regions for display: an ellipse that would extend beyond the
  # response scale is clipped to the plotting limits, so its outline runs
  # along the axis where it would leave the scale. The unclipped outline is
  # kept in local_n_region_df (returned in the output). No observed profile
  # lies outside the scale, so the clipped region contains the same cases.
  local_n_region_plot <- local_n_region_df
  if (nrow(local_n_region_plot) > 0) {
    local_n_region_plot <- local_n_region_plot %>%
      dplyr::mutate(
        x = pmin(pmax(x, x_limits[1]), x_limits[2]),
        y = pmin(pmax(y, y_limits[1]), y_limits[2])
      )
  }
  
  if (!is.numeric(label_offset) || length(label_offset) != 1 || is.na(label_offset) || label_offset < 0) {
    stop("label_offset must be a single numeric value >= 0.")
  }
  label_offset <- as.numeric(label_offset)
  
  if (!is.numeric(profile_label_size) || length(profile_label_size) != 1 ||
      is.na(profile_label_size) || profile_label_size < 0) {
    stop("profile_label_size must be a single numeric value >= 0.")
  }
  profile_label_size <- as.numeric(profile_label_size)
  
  # Label placement: offset of label_offset (as a fraction of the axis
  # ranges) in the direction given by label_dir for each profile. The default
  # direction c(-1, 0) with right-aligned text reproduces the original
  # placement to the left of the circle. Horizontal alignment follows the
  # direction (right-aligned when placed left, left-aligned when placed right,
  # centered when placed straight above or below).
  x_off <- label_offset * diff(x_limits)
  y_off <- label_offset * diff(y_limits)
  
  dir_x <- vapply(profile_pts$label, function(l) {
    if (!is.null(label_dir[[l]])) label_dir[[l]][1] else label_dir_default[1]
  }, numeric(1))
  dir_y <- vapply(profile_pts$label, function(l) {
    if (!is.null(label_dir[[l]])) label_dir[[l]][2] else label_dir_default[2]
  }, numeric(1))
  
  label_df <- profile_pts %>%
    dplyr::mutate(
      label_x = x + dir_x * x_off,
      label_y = y + dir_y * y_off,
      label_hjust = ifelse(dir_x > 0, 0, ifelse(dir_x < 0, 1, 0.5)),
      label_txt = label
    )
  
  hex_to_rgb01 <- function(hex) as.numeric(grDevices::col2rgb(hex) / 255)
  rgb01_to_hex <- function(r, g, b) grDevices::rgb(r, g, b)
  mix_hex <- function(from_hex, to_hex, t) {
    t <- pmin(pmax(t, 0), 1)
    from <- hex_to_rgb01(from_hex)
    to   <- hex_to_rgb01(to_hex)
    rgb01_to_hex(
      from[1] + (to[1] - from[1]) * t,
      from[2] + (to[2] - from[2]) * t,
      from[3] + (to[3] - from[3]) * t
    )
  }
  
  neg_col <- "#556359"
  pos_col <- "#755E69"
  white0  <- "#f3f3f3"
  zmax_out <- max(abs(z_lim))
  fade_pow <- 0.60
  
  df_out_col <- df_out %>%
    dplyr::mutate(
      pred_for_col = ifelse(is.finite(predicted), predicted, 0),
      t = pmin(1, (abs(pred_for_col) / zmax_out))^fade_pow,
      outside_col = ifelse(
        pred_for_col < 0,
        vapply(t, function(tt) mix_hex(white0, neg_col, tt), character(1)),
        vapply(t, function(tt) mix_hex(white0, pos_col, tt), character(1))
      )
    )
  
  custom_colors <- c("#0d4a1b", "#0d4a1b", "#258e3d", "#71da89", "#ffffff",
                     "#ee65a3", "#ab0e63", "#7a0745", "#7a0745")
  color_breaks <- c(-2.1, -1.5, -1.0, -0.5, 0, 0.5, 1, 1.5, 2.1)
  keep <- color_breaks >= -2 & color_breaks <= 2
  colors_legend <- custom_colors[keep]
  breaks_legend <- color_breaks[keep]
  
  # Response labels of the CAMSQ (from never [1] to always [7]).
  tick_labels <- c(
    "(1)\nNever",
    "(2)\nAlmost\nnever",
    "(3)\nSome-\ntimes",
    "(4)\nHalf of\nthe time",
    "(5)\nOften",
    "(6)\nAlmost\nalways",
    "(7)\nAlways"
  )
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  # 13. Build plot                                                            ####
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  
  # Response surface plot:
  # predicted values shown fully within the selected HDR (dashed contour),
  # attenuated outside; circles mark the empirically adjusted profile
  # coordinates, colored and labeled by their predicted value.
  p <- ggplot2::ggplot()
  
  if (isTRUE(show_outside)) {
    if (!isTRUE(greyscale)) {
      p <- p +
        ggplot2::geom_tile(
          data = df_out_col,
          ggplot2::aes(x = x, y = y, fill = outside_col),
          alpha = 1,
          inherit.aes = FALSE,
          show.legend = FALSE
        ) +
        ggplot2::scale_fill_identity() +
        ggnewscale::new_scale_fill()
    } else {
      p <- p +
        ggplot2::geom_tile(
          data = df_out,
          ggplot2::aes(x = x, y = y),
          fill = "grey88",
          alpha = 1,
          inherit.aes = FALSE
        )
    }
  }
  
  p <- p +
    ggplot2::geom_tile(
      data = df_in,
      ggplot2::aes(x = x, y = y, fill = fill_val),
      alpha = 0.95,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_contour(
      data = grid_plot,
      ggplot2::aes(x = x, y = y, z = density),
      breaks = hdr_cutoff,
      color = "gray25",
      linewidth = 0.8,
      linetype = "dashed",
      inherit.aes = FALSE
    )
  
  if (isTRUE(show_local_n_region) && nrow(local_n_region_plot) > 0) {
    p <- p +
      ggplot2::geom_polygon(
        data = local_n_region_plot,
        ggplot2::aes(x = x, y = y, group = label),
        fill = "white",
        color = "black",
        alpha = local_n_region_alpha,
        linewidth = local_n_region_linewidth,
        linetype = local_n_region_linetype,
        inherit.aes = FALSE,
        show.legend = FALSE
      )
  }
  
  if (isTRUE(plot_scatter)) {
    p <- p +
      ggplot2::geom_point(
        data = observed_points,
        ggplot2::aes(x = x, y = y),
        color = "black",
        alpha = scatter_alpha,
        size = scatter_size,
        inherit.aes = FALSE
      )
  }
  
  if (!identical(plot_profiles, FALSE)) {
    profile_pts_plot <- profile_pts
    label_df_plot <- label_df
    
    if (is.character(plot_profiles) && identical(plot_profiles, "constructs")) {
      keep_labs <- c("Hypo", "Hyper", "Effective")
      profile_pts_plot <- profile_pts_plot %>% dplyr::filter(label %in% keep_labs)
      label_df_plot <- label_df_plot %>% dplyr::filter(label_txt %in% keep_labs)
    }
    
    p <- p +
      ggplot2::geom_point(
        data = profile_pts_plot,
        ggplot2::aes(x = x, y = y, fill = fill_val),
        shape = 21, size = 7, color = "black", stroke = 0.7,
        inherit.aes = FALSE
      ) +
      ggplot2::geom_text(
        data = profile_pts_plot,
        ggplot2::aes(x = x, y = y, label = pred_label),
        size = 2, color = "black", inherit.aes = FALSE
      )
    
    construct_labs <- c("Hypo", "Hyper", "Effective")
    label_df_construct <- label_df_plot %>% dplyr::filter(label_txt %in% construct_labs)
    label_df_nonconstruct <- label_df_plot %>% dplyr::filter(!label_txt %in% construct_labs)
    
    if (nrow(label_df_nonconstruct) > 0) {
      p <- p +
        ggplot2::geom_label(
          data = label_df_nonconstruct,
          ggplot2::aes(x = label_x, y = label_y, label = label_txt, hjust = label_hjust),
          size = profile_label_size,
          color = "black",
          fill = "white",
          label.r = grid::unit(0.15, "lines"),
          vjust = 0.5,
          inherit.aes = FALSE
        )
    }
    
    if (nrow(label_df_construct) > 0) {
      p <- p +
        ggplot2::geom_label(
          data = label_df_construct,
          ggplot2::aes(x = label_x, y = label_y, label = label_txt, hjust = label_hjust),
          size = profile_label_size,
          color = "black",
          fill = "white",
          label.r = grid::unit(0.15, "lines"),
          vjust = 0.5,
          inherit.aes = FALSE
        )
    }
  }
  
  if (isTRUE(show_r2)) {
    clean_model_name <- sub("^.*: ", "", selected_model_name)
    
    # Corner of the R2 box (3% inset from the plot limits); the default
    # "topleft" is the original placement.
    r2_left   <- r2_position %in% c("topleft", "bottomleft")
    r2_top    <- r2_position %in% c("topleft", "topright")
    r2_x      <- if (r2_left) x_limits[1] + 0.03 * diff(x_limits) else x_limits[2] - 0.03 * diff(x_limits)
    r2_y      <- if (r2_top)  y_limits[2] - 0.03 * diff(y_limits) else y_limits[1] + 0.03 * diff(y_limits)
    r2_hjust  <- if (r2_left) 0 else 1
    r2_vjust  <- if (r2_top)  1 else 0
    
    if (identical(model, "loess")) {
      p <- p + ggplot2::annotate(
        "label",
        x = r2_x, y = r2_y,
        label = sprintf("LOESS R² = %.3f\nAdaptive HDR = %.1f%%",
                        r2_tbl$r2[1],
                        selected_hdr),
        size = r2_text_size, color = "black", fill = "white",
        label.r = grid::unit(0.15, "lines"),
        hjust = r2_hjust, vjust = r2_vjust
      )
    } else if (identical(estimator, "ols")) {
      
      if (isTRUE(cv_r2_needed)) {
        # Cross-validation was computed: annotate with CV-R².
        selected_cv_r2 <- r2_tbl$cv_r2[match(clean_model_name, r2_tbl$model)]
        
        if (length(selected_cv_r2) != 1 || !is.finite(selected_cv_r2)) {
          selected_cv_r2 <- r2_tbl$cv_r2[1]
        }
        
        annot_label <- sprintf("CV-R² = %.3f (%s)\nAdaptive HDR = %.1f%%",
                               selected_cv_r2,
                               clean_model_name,
                               selected_hdr)
      } else {
        # Forced model with cv_r2 = FALSE: annotate with in-sample R² instead,
        # since no cross-validated R² was computed.
        selected_r2 <- r2_tbl$r2[match(clean_model_name, r2_tbl$model)]
        
        if (length(selected_r2) != 1 || !is.finite(selected_r2)) {
          selected_r2 <- r2_tbl$r2[1]
        }
        
        annot_label <- sprintf("R² = %.3f (%s)\nAdaptive HDR = %.1f%%",
                               selected_r2,
                               clean_model_name,
                               selected_hdr)
      }
      
      p <- p + ggplot2::annotate(
        "label",
        x = r2_x, y = r2_y,
        label = annot_label,
        size = r2_text_size, color = "black", fill = "white",
        label.r = grid::unit(0.15, "lines"),
        hjust = r2_hjust, vjust = r2_vjust
      )
    } else {
      bayes_r2_mean <- get_bayes_r2(fit)
      p <- p + ggplot2::annotate(
        "label",
        x = r2_x, y = r2_y,
        label = sprintf("Bayesian R² = %.3f (%s)\nAdaptive HDR = %.1f%%",
                        bayes_r2_mean,
                        clean_model_name,
                        selected_hdr),
        size = r2_text_size, color = "black", fill = "white",
        label.r = grid::unit(0.15, "lines"),
        hjust = r2_hjust, vjust = r2_vjust
      )
    }
  }
  
  if (isTRUE(show_n)) {
    # Separate box with the number of complete cases (default "bottomright"),
    # independent of the R2 box.
    n_left <- n_position %in% c("topleft", "bottomleft")
    n_top  <- n_position %in% c("topleft", "topright")
    n_x    <- if (n_left) x_limits[1] + 0.03 * diff(x_limits) else x_limits[2] - 0.03 * diff(x_limits)
    n_y    <- if (n_top)  y_limits[2] - 0.03 * diff(y_limits) else y_limits[1] + 0.03 * diff(y_limits)
    p <- p + ggplot2::annotate(
      "label",
      x = n_x, y = n_y,
      label = sprintf("N = %d", n_used),
      size = r2_text_size, color = "black", fill = "white",
      label.r = grid::unit(0.15, "lines"),
      hjust = if (n_left) 0 else 1, vjust = if (n_top) 1 else 0
    )
  }
  
  if (!isTRUE(greyscale)) {
    p <- p +
      ggplot2::scale_fill_gradientn(
        colors = colors_legend,
        values = scales::rescale(breaks_legend),
        limits = z_lim,
        breaks  = z_breaks,
        labels  = z_breaks,
        oob     = scales::squish,
        name    = "Predicted\nValue (z)"
      )
  } else {
    transform_fun_scale <- function(z) sign(z) * abs(z)^grey_power
    lim_trans    <- transform_fun_scale(z_lim)
    breaks_trans <- transform_fun_scale(z_breaks)
    p <- p +
      ggplot2::scale_fill_gradient2(
        low = "grey10", mid = "white", high = "grey10",
        midpoint = 0,
        limits = range(lim_trans),
        breaks = breaks_trans,
        labels = z_breaks,
        oob = scales::squish,
        name = "Predicted\nValue (z)"
      )
  }
  
  p <- p +
    ggplot2::scale_x_continuous(limits = x_limits, breaks = 1:7, labels = tick_labels, expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = y_limits, breaks = 1:7, labels = tick_labels, expand = c(0, 0)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_minimal()
  
  if (!is.null(title) && nzchar(title)) {
    p <- p + ggplot2::labs(title = title, x = "CAMSQ Self-Certainty", y = "CAMSQ Other-Certainty")
  } else {
    p <- p + ggplot2::labs(x = "CAMSQ Self-Certainty", y = "CAMSQ Other-Certainty")
  }
  
  # Rendering warnings from ggplot2 (e.g., about tiles at the plot boundary)
  # are suppressed so that the function emits no warnings; this is scoped to
  # printing the plot and does not affect any statistical output.
  suppressWarnings(print(p))
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  # 14. Print coordinates and model output                                    ####
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  empirical_coords <- profile_pts %>%
    dplyr::mutate(
      local_n = as.integer(local_n)
    ) %>%
    dplyr::select(x, y, label, local_n)
  
  empirical_coords_print <- empirical_coords %>%
    dplyr::mutate(
      x = round(x, 3),
      y = round(y, 3)
    )
  
  message("\nEmpirical CAMSQ profile coordinates (rounded to 3 d.p.):")
  for (i in seq_len(nrow(empirical_coords_print))) {
    message(sprintf("  %-15s  x = %.3f,  y = %.3f,  local N = %.0f",
                    empirical_coords_print$label[i],
                    empirical_coords_print$x[i],
                    empirical_coords_print$y[i],
                    empirical_coords_print$local_n[i]))
  }
  
  if (isTRUE(adaptive_hdr)) {
    message(sprintf(
      "\nAdaptive highest-density region selected: %.1f%%.",
      selected_hdr
    ))
    message(sprintf(
      "Minimum required local N: %.0f for profiles: %s.",
      min_local_n,
      paste(local_n_profiles, collapse = ", ")
    ))
    if (!is.null(hdr_min)) {
      message(sprintf(
        "Lower bound applied to the HDR search: %.1f%%.",
        hdr_min
      ))
    } else {
      message("No lower bound applied to the HDR search (searched down to the smallest value in hdr_search).")
    }
    message(sprintf(
      "Local N Mahalanobis distance threshold: %.2f covariance-standardized units.",
      local_n_radius
    ))
  } else {
    message(sprintf(
      "\nHighest-density region used: %.1f%%.",
      selected_hdr
    ))
    message(sprintf(
      "Local N Mahalanobis distance threshold: %.2f covariance-standardized units.",
      local_n_radius
    ))
  }
  
  if (isTRUE(show_local_n_region)) {
    message(sprintf(
      "Displayed hard-count local-N regions: profiles = %s, Mahalanobis distance threshold = %.2f.",
      paste(local_n_region_profiles, collapse = ", "),
      local_n_radius
    ))
  }
  
  message(sprintf("\nCases used for all analyses (KDE, RSA, local N): %d complete cases.\n", n_used))
  
  if (isTRUE(cv_r2_needed)) {
    cv_msg_suffix <- if (identical(model, "polynomial") && identical(estimator, "ols") && identical(poly_model, "auto")) {
      " for automatic model selection"
    } else {
      ""
    }
    message(sprintf(
      "Cross-validated R² estimated using %d-fold cross-validation%s.\n",
      cv_folds,
      cv_msg_suffix
    ))
  } else if (identical(model, "polynomial") && identical(estimator, "ols")) {
    message(sprintf(
      "Cross-validated R² was not computed (poly_model = '%s', cv_r2 = FALSE).\n",
      poly_model
    ))
  }
  
  if (!is.null(title) && nzchar(title)) message(sprintf("Criterion: %s\n", title))
  
  if (identical(model, "loess")) {
    
    message("\n--- LOESS surface metrics ---")
    message(sprintf("LOESS: R²=%.3f | span=%.3f | degree=%d | family=%s",
                    r2_tbl$r2[1],
                    loess_span,
                    loess_degree,
                    loess_family))
    
    message(sprintf("\nModel used for Plot/Profile Predictions: **%s**",
                    sub("^.*: ", "", selected_model_name)))
    message("Note: Cross-validation, polynomial model comparison, and emmeans are not performed when model = 'loess'.")
    message("The poly_model and estimator arguments are ignored when model = 'loess'.\n")
    
    message("\n--- Selected LOESS model summary ---")
    print(summary(fit))
    
  } else if (identical(estimator, "ols")) {
    
    if (identical(poly_model, "auto")) {
      message("\n--- Model comparison: Linear vs Quadratic vs Cubic ---")
      message(sprintf("Linear:    R²=%.3f%s | Adj.R²=%.3f | AIC=%.1f | BIC=%.1f | AICc=%.1f",
                      r2_tbl$r2[r2_tbl$model == "Linear"],
                      cv_txt(r2_tbl$cv_r2[r2_tbl$model == "Linear"]),
                      r2_tbl$adj_r2[r2_tbl$model == "Linear"],
                      r2_tbl$aic[r2_tbl$model == "Linear"],
                      r2_tbl$bic[r2_tbl$model == "Linear"],
                      r2_tbl$aicc[r2_tbl$model == "Linear"]))
      message(sprintf("Quadratic: R²=%.3f%s | Adj.R²=%.3f | AIC=%.1f | BIC=%.1f | AICc=%.1f",
                      r2_tbl$r2[r2_tbl$model == "Quadratic"],
                      cv_txt(r2_tbl$cv_r2[r2_tbl$model == "Quadratic"]),
                      r2_tbl$adj_r2[r2_tbl$model == "Quadratic"],
                      r2_tbl$aic[r2_tbl$model == "Quadratic"],
                      r2_tbl$bic[r2_tbl$model == "Quadratic"],
                      r2_tbl$aicc[r2_tbl$model == "Quadratic"]))
      message(sprintf("Cubic:     R²=%.3f%s | Adj.R²=%.3f | AIC=%.1f | BIC=%.1f | AICc=%.1f",
                      r2_tbl$r2[r2_tbl$model == "Cubic"],
                      cv_txt(r2_tbl$cv_r2[r2_tbl$model == "Cubic"]),
                      r2_tbl$adj_r2[r2_tbl$model == "Cubic"],
                      r2_tbl$aic[r2_tbl$model == "Cubic"],
                      r2_tbl$bic[r2_tbl$model == "Cubic"],
                      r2_tbl$aicc[r2_tbl$model == "Cubic"]))
      
      message(sprintf("\nNested F-test (Linear vs Quadratic): p %s %s",
                      fmt_p(p_lin_quad), sig_star(p_lin_quad)))
      message(sprintf("Nested F-test (Quadratic vs Cubic):  p %s %s",
                      fmt_p(p_quad_cub), sig_star(p_quad_cub)))
    } else {
      message("\n--- Model metrics (requested model only) ---")
      message(sprintf("%s: R²=%.3f%s | Adj.R²=%.3f | AIC=%.1f | BIC=%.1f | AICc=%.1f",
                      r2_tbl$model[1],
                      r2_tbl$r2[1],
                      cv_txt(r2_tbl$cv_r2[1]),
                      r2_tbl$adj_r2[1],
                      r2_tbl$aic[1],
                      r2_tbl$bic[1],
                      r2_tbl$aicc[1]))
    }
    
    message(sprintf("\nModel used for Plot/EMMs: **%s**", sub("^.*: ", "", selected_model_name)))
    if (identical(poly_model, "auto")) {
      message("(Selected automatically based on highest cross-validated R²)")
      message("Note: Higher cross-validated R² indicates better out-of-sample predictive performance.\n")
    } else {
      message("(Forced by the user via the poly_model argument; no automatic model selection performed)")
    }
    
    message("\n--- Selected regression model summary ---")
    print(summary(fit))
    
  } else {
    
    message("\n--- Model metrics (requested model only, Bayesian) ---")
    message(sprintf("%s: Bayesian R²=%.3f | WAIC=%.1f",
                    comp_tbl$model[1], comp_tbl$bayes_r2[1], comp_tbl$waic[1]))
    
    message(sprintf("\nModel used for Plot/EMMs: **%s**", sub("^.*: ", "", selected_model_name)))
    message("Note: Lower WAIC indicates better expected out-of-sample predictive accuracy.\n")
    message("\n--- Selected Bayesian regression model summary ---")
    print(summary(fit))
  }
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  # 15. Profile predictions and contrast definitions                          ####
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  
  profile_levels_core <- c("Median", "Hypo", "Hyper", "Effective", "High-Certainty")
  
  profile_predictions <- profile_pts %>%
    dplyr::mutate(
      profile = factor(label, levels = c(profile_levels_core, "Max", "Min")),
      mipv_z = pred,
      mipv_original = pred_orig,
      local_n = as.integer(local_n)
    ) %>%
    dplyr::arrange(profile) %>%
    dplyr::select(profile, x, y, mipv_z, mipv_original, local_n)
  
  found_labels <- as.character(profile_predictions$profile)
  
  model_comparison_out <- if (identical(model, "loess")) {
    list(
      metrics = r2_tbl,
      p_lin_quad = NA_real_,
      p_quad_cub = NA_real_
    )
  } else if (identical(estimator, "ols")) {
    list(
      metrics = r2_tbl,
      p_lin_quad = p_lin_quad,
      p_quad_cub = p_quad_cub
    )
  } else {
    list(
      metrics = comp_tbl
    )
  }
  
  if (!all(profile_levels_core %in% found_labels)) {
    message("Not all core profiles could be placed within the selected HDR; contrasts are not computed.")
    return(invisible(list(
      ref_grid = NULL,
      profile_predictions = profile_predictions,
      results = NULL,
      fit = fit,
      empirical_coords = empirical_coords,
      local_support = local_support,
      local_n_region = local_n_region_df,
      plot = p,
      model_comparison = model_comparison_out,
      selected_model = sub("^.*: ", "", selected_model_name),
      model = model,
      estimator = estimator,
      show_min_max = show_min_max,
      hdr = list(
        level = selected_hdr,
        cutoff = hdr_cutoff,
        adaptive = adaptive_hdr,
        min_local_n = min_local_n,
        hdr_min = hdr_min,
        local_n_profiles = local_n_profiles,
        local_n_radius = local_n_radius,
        local_n_metric = "mahalanobis",
        local_n_covariance = S_local,
        diagnostics = hdr_diagnostics
      ),
      local_n_region_settings = list(
        shown = show_local_n_region,
        profiles = local_n_region_profiles,
        radius = local_n_radius,
        metric = "mahalanobis",
        alpha = local_n_region_alpha,
        linewidth = local_n_region_linewidth,
        linetype = local_n_region_linetype
      )
    )))
  }
  
  if (isTRUE(show_min_max)) {
    if ("Max" %in% found_labels) {
      max_row <- profile_predictions %>% dplyr::filter(profile == "Max")
      message(sprintf(
        "Max: empirical maximum of the predicted surface within the HDR at x = %.3f, y = %.3f (predicted z = %.3f).",
        max_row$x[1], max_row$y[1], max_row$mipv_z[1]
      ))
    }
    if ("Min" %in% found_labels) {
      min_row <- profile_predictions %>% dplyr::filter(profile == "Min")
      message(sprintf(
        "Min: empirical minimum of the predicted surface within the HDR at x = %.3f, y = %.3f (predicted z = %.3f).",
        min_row$x[1], min_row$y[1], min_row$mipv_z[1]
      ))
    }
  }
  
  print(profile_predictions)
  
  # The nine pairwise comparisons,
  # expressed as contrast weights over the profile order
  # (Median, Hypo, Hyper, Effective, High-Certainty).
  # The comparisons are grouped by focal construct (effective mentalizing,
  # hypermentalizing, hypomentalizing); comparisons appearing in two groups
  # (e.g., Hyper_vs_Effective and Effective_vs_Hyper) are sign-reversed
  # duplicates of the same contrast, so the three groups jointly cover the
  # nine unique comparisons.
  effective_test <- list(
    "Effective_vs_High-Certainty" = c(0, 0, 0, 1, -1),
    "Effective_vs_Median"         = c(-1, 0, 0, 1, 0),
    "Effective_vs_Hypo"           = c(0, -1, 0, 1, 0),
    "Effective_vs_Hyper"          = c(0, 0, -1, 1, 0)
  )
  
  hyper_test <- list(
    "Hyper_vs_High-Certainty" = c(0, 0, 1, 0, -1),
    "Hyper_vs_Median"         = c(-1, 0, 1, 0, 0),
    "Hyper_vs_Hypo"           = c(0, -1, 1, 0, 0),
    "Hyper_vs_Effective"      = c(0, 0, 1, -1, 0)
  )
  
  hypo_test <- list(
    "Hypo_vs_High-Certainty" = c(0, 1, 0, 0, -1),
    "Hypo_vs_Median"         = c(-1, 1, 0, 0, 0),
    "Hypo_vs_Hyper"          = c(0, 1, -1, 0, 0),
    "Hypo_vs_Effective"      = c(0, 1, 0, -1, 0)
  )
  
  if (identical(model, "loess")) {
    
    # LOESS contrasts: differences in
    # predicted values without standard errors or p values.
    make_loess_contrasts <- function(contrast_list, label) {
      est <- purrr::map_dbl(contrast_list, function(w) {
        sum(w * profile_predictions$mipv_z[match(profile_levels_core, profile_predictions$profile)])
      })
      est_orig <- if (all(is.finite(profile_predictions$mipv_original[match(profile_levels_core, profile_predictions$profile)]))) {
        purrr::map_dbl(contrast_list, function(w) {
          sum(w * profile_predictions$mipv_original[match(profile_levels_core, profile_predictions$profile)])
        })
      } else {
        rep(NA_real_, length(contrast_list))
      }
      
      out <- tibble::tibble(
        contrast = names(contrast_list),
        estimate = est,
        estimate_original = est_orig
      )
      
      message(sprintf("\n--- %s (LOESS: descriptive differences of predicted values; no SEs/p-values) ---", label))
      print(as.data.frame(out), row.names = FALSE)
      out
    }
    
    eff_res   <- make_loess_contrasts(effective_test, "Effective Mentalizing Contrasts")
    hyper_res <- make_loess_contrasts(hyper_test, "Hypermentalizing Contrasts")
    hypo_res  <- make_loess_contrasts(hypo_test, "Hypomentalizing Contrasts")
    
    return(invisible(list(
      ref_grid = NULL,
      profile_predictions = profile_predictions,
      results = list(
        effective_mentalizing = eff_res,
        hypermentalizing = hyper_res,
        hypomentalizing = hypo_res
      ),
      fit = fit,
      empirical_coords = empirical_coords,
      local_support = local_support,
      local_n_region = local_n_region_df,
      plot = p,
      model_comparison = model_comparison_out,
      selected_model = sub("^.*: ", "", selected_model_name),
      model = model,
      estimator = estimator,
      show_min_max = show_min_max,
      hdr = list(
        level = selected_hdr,
        cutoff = hdr_cutoff,
        adaptive = adaptive_hdr,
        min_local_n = min_local_n,
        hdr_min = hdr_min,
        local_n_profiles = local_n_profiles,
        local_n_radius = local_n_radius,
        local_n_metric = "mahalanobis",
        local_n_covariance = S_local,
        diagnostics = hdr_diagnostics
      ),
      local_n_region_settings = list(
        shown = show_local_n_region,
        profiles = local_n_region_profiles,
        radius = local_n_radius,
        metric = "mahalanobis",
        alpha = local_n_region_alpha,
        linewidth = local_n_region_linewidth,
        linetype = local_n_region_linetype
      )
    )))
  }
  
  # Reference grid at the empirically adjusted profile coordinates: the
  # contrasts below compare model-implied predicted values (MIPVs) at these
  # coordinates on the standardized criterion metric.
  empirical_coords_core <- empirical_coords %>%
    dplyr::filter(label %in% profile_levels_core) %>%
    dplyr::mutate(label = factor(label, levels = profile_levels_core)) %>%
    dplyr::arrange(label)
  
  self_coords  <- unique(empirical_coords_core$x)
  other_coords <- unique(empirical_coords_core$y)
  
  grid_at_list <- list(self_coords, other_coords)
  names(grid_at_list) <- c(sc_nm, oc_nm)
  
  rg_full <- emmeans::ref_grid(fit, at = grid_at_list)
  
  # Match the profile coordinates in the reference grid using a numerical
  # tolerance. The coordinates are passed to emmeans unchanged, so the
  # tolerance only guards against floating-point representation issues.
  profile_row_idx <- vapply(
    seq_len(nrow(empirical_coords_core)),
    function(i) {
      idx <- which(
        abs(rg_full@grid[[sc_nm]] - empirical_coords_core$x[i]) < 1e-8 &
          abs(rg_full@grid[[oc_nm]] - empirical_coords_core$y[i]) < 1e-8
      )
      if (length(idx) != 1) {
        stop("Internal error: empirical profile coordinates could not be uniquely matched in the emmeans reference grid.")
      }
      idx
    },
    integer(1)
  )
  
  rg_profiles <- rg_full[profile_row_idx]
  
  rg_profiles@grid$profile <- factor(
    empirical_coords_core$label,
    levels = profile_levels_core
  )
  
  if (!identical(as.character(rg_profiles@grid$profile), as.character(empirical_coords_core$label))) {
    stop("Internal error: profile ordering mismatch between empirical coordinates and the emmeans reference grid.")
  }
  
  message("\n--- Reference grid at empirical profile coordinates ---")
  print(rg_profiles)
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  # 16. Estimate contrasts                                                    ####
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  
  if (identical(estimator, "ols")) {
    
    # Exact unadjusted p values
    # with a two-sided alpha of .05. 
    eff_res <- emmeans::contrast(rg_profiles, method = effective_test) %>%
      summary(infer = c(TRUE, TRUE), adjust = "none")
    message("\n--- Effective Mentalizing Contrasts (OLS) ---")
    print(eff_res)
    
    hyper_res <- emmeans::contrast(rg_profiles, method = hyper_test) %>%
      summary(infer = c(TRUE, TRUE), adjust = "none")
    message("\n--- Hypermentalizing Contrasts (OLS) ---")
    print(hyper_res)
    
    hypo_res <- emmeans::contrast(rg_profiles, method = hypo_test) %>%
      summary(infer = c(TRUE, TRUE), adjust = "none")
    message("\n--- Hypomentalizing Contrasts (OLS) ---")
    print(hypo_res)
    
  } else {
    
    # Bayesian contrasts: posterior
    # medians with 95% HPD intervals instead of frequentist p values.
    summarize_contrast <- function(rg, contrast_list, label) {
      ctr <- emmeans::contrast(rg, method = contrast_list)
      out <- summary(ctr, point.est = stats::median, level = 0.95)
      
      nm <- names(out)
      nm[nm == "lower.HPD"] <- "lower_hpd"
      nm[nm == "upper.HPD"] <- "upper_hpd"
      nm[nm == "lower.CL"] <- "lower_hpd"
      nm[nm == "upper.CL"] <- "upper_hpd"
      names(out) <- nm
      
      message(sprintf("\n--- %s (Bayesian: posterior median and 95%% HPD interval) ---", label))
      print(out)
      out
    }
    
    eff_res   <- summarize_contrast(rg_profiles, effective_test, "Effective Mentalizing Contrasts")
    hyper_res <- summarize_contrast(rg_profiles, hyper_test, "Hypermentalizing Contrasts")
    hypo_res  <- summarize_contrast(rg_profiles, hypo_test, "Hypomentalizing Contrasts")
  }
  
  invisible(list(
    ref_grid = rg_profiles,
    profile_predictions = profile_predictions,
    results = list(
      effective_mentalizing = eff_res,
      hypermentalizing = hyper_res,
      hypomentalizing = hypo_res
    ),
    fit = fit,
    empirical_coords = empirical_coords,
    local_support = local_support,
    local_n_region = local_n_region_df,
    plot = p,
    model_comparison = model_comparison_out,
    selected_model = sub("^.*: ", "", selected_model_name),
    model = model,
    estimator = estimator,
    show_min_max = show_min_max,
    hdr = list(
      level = selected_hdr,
      cutoff = hdr_cutoff,
      adaptive = adaptive_hdr,
      min_local_n = min_local_n,
      hdr_min = hdr_min,
      local_n_profiles = local_n_profiles,
      local_n_radius = local_n_radius,
      local_n_metric = "mahalanobis",
      local_n_covariance = S_local,
      diagnostics = hdr_diagnostics
    ),
    local_n_region_settings = list(
      shown = show_local_n_region,
      profiles = local_n_region_profiles,
      radius = local_n_radius,
      metric = "mahalanobis",
      alpha = local_n_region_alpha,
      linewidth = local_n_region_linewidth,
      linetype = local_n_region_linetype
    )
  ))
}