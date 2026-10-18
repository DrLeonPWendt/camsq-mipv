# camsq-mipv

`camsq_mipv_analysis()` is an R function that models a criterion variable as a function of the two scales of the Certainty About Mental States Questionnaire (CAMSQ), Self-Certainty (SC) and Other-Certainty (OC). It calculates model-implied predicted values (MIPVs) of the criterion across the entire response surface spanned by SC and OC, and it tests differences in MIPVs between five theoretically defined CAMSQ profiles, which represent positions within this two-dimensional space:

- **Effective**: high SC, moderate OC
- **Hypo**: low SC, low OC
- **Hyper**: low SC, high OC
- **High-Certainty**: high SC, high OC
- **Median**: moderate SC and OC

MIPVs are based on polynomial regression models estimated within a response surface analysis (RSA) framework. Three candidate models are considered: a linear model (additive linear effects of SC and OC), a quadratic model (adding SC², OC², SC×OC), and a cubic model (adding SC³, SC²×OC, SC×OC², OC³). The model with the highest 10-fold cross-validated R² (CV-R²) is selected, and the criterion is z-standardized so that differences in MIPVs are on a standardized metric. Because the theoretical profile coordinates partly lie at scale boundaries that are rarely observed, empirical adjustments are applied: nonparametric kernel density estimation determines the highest density region (HDR) of the observed CAMSQ scores, and for each profile the coordinate within the HDR closest to its theoretical location (Mahalanobis distance) is selected. The HDR is chosen adaptively, starting at the 99%-HDR and examining increasingly restrictive HDRs in steps of 1%, selecting the largest HDR in which all profiles reach a minimum local sample size, defined as the number of observations within an elliptical neighborhood (Mahalanobis distance 0.50) around each adjusted coordinate. The function returns the MIPVs at the adjusted profile coordinates with confidence intervals, significance tests for all pairwise contrasts between profiles, the model comparison, the fitted model, and a plot of the response surface with the adjusted profiles located in the observed data.

The function was developed for and is documented in

> Wendt, L. P., Kiel, L., Fonagy, P., Lind, M., Lukowitsky, M., Nolte, T., Wester, R., Rubel, J., Montague, P. R., & Zimmermann, J. (2026). *Revisiting the construct validity of self-reported mentalizing.* Manuscript in preparation.

Please cite this paper when using the function.

## Installation

No package installation is needed. Download `camsq_mipv_analysis.R` from this repository (or from a specific release) and source it:

```r
source("camsq_mipv_analysis.R")
```

The function loads the following packages and stops with a message if one is missing: `ks`, `ggplot2`, `dplyr`, `tibble`, `purrr`, `rlang`, `scales`, `grid`, `ggnewscale`, `emmeans`. Install them once with `install.packages()`, or call the function with `auto_install = TRUE`.

## Minimal example

A minimal call needs a data frame, the name of the criterion variable, and the names of the two CAMSQ scale scores (item means, range 1 to 7):

```r
res <- camsq_mipv_analysis(
  data        = dat,
  criterion   = "PDS_tot",
  sc          = "CAMSQ_sc",
  oc          = "CAMSQ_oc",
  min_local_n = 20,
  title       = "Personality functioning (PDS-ICD-11)"
)

res$results               # profile contrasts (estimate in criterion SD units, SE, 95% CI, p), three groups
res$profile_predictions   # MIPVs of the five profiles with 95% CIs
res$model_comparison      # CV-R², R², AIC, BIC, and nested model tests of the three models
res$hdr$level             # selected HDR level
res$fit                   # lm object of the selected model
res$plot                  # ggplot object, e.g., ggsave("plot.png", res$plot)
```

`min_local_n` is the minimum local sample size. Based on a power analysis, we recommend 20 for analyses with N ≥ 300 and 10 for smaller samples.

## Full example

The same analysis with some useful arguments written out at their default values, each with a short explanation:

```r
res <- camsq_mipv_analysis(

  # data and variables
  data                = dat,                 # data frame
  criterion           = "PDS_tot",           # criterion variable
  sc                  = "CAMSQ_sc",          # CAMSQ Self-Certainty (mean of the ten items, range 1 to 7)
  oc                  = "CAMSQ_oc",          # CAMSQ Other-Certainty (mean of the ten items, range 1 to 7)

  # profiles and highest density region (HDR)
  min_local_n         = 20,                  # minimum N within the local neighborhood of each profile (NULL fixes the HDR at hdr)
  hdr                 = 99,                  # HDR level in percent at which the search starts (or the fixed level if min_local_n = NULL)
  hdr_min             = NULL,                # lower bound of the search, NULL = none
  local_n_radius      = 0.50,                # Mahalanobis distance that defines the elliptical local neighborhood of a profile
  local_n_profiles    = c("Median",
                          "Hypo",
                          "Hyper",
                          "Effective",
                          "High-Certainty"), # profiles that must reach min_local_n

  # model estimation
  poly_model          = "auto",              # "auto" selects the model by cross-validated R²,
                                             # "linear", "quadratic", or "cubic" forces a model
  cv_folds            = 10,                  # folds of the cross-validation
  cv_seed             = 1234,                # seed of the fold assignment
  estimator           = "ols",               # ordinary least squares

  # plot: general
  title               = "Personality functioning (PDS-ICD-11)",   # plot title
  sc_range            = c(1, 7),             # score range of SC over which the surface is displayed
  oc_range            = c(1, 7),             # score range of OC over which the surface is displayed
  grid_n              = 300,                 # resolution of the prediction grid (300 x 300)
  greyscale           = FALSE,               # greyscale instead of color

  # plot: scatterplot
  plot_scatter        = TRUE,                # show the observations
  scatter_alpha       = 0.1,                 # transparency of the observations
  scatter_size        = 1.5,                 # size of the observations
  show_outside        = TRUE,                # draw the attenuated surface outside the HDR

  # plot: white box 1 (CV-R² and HDR)
  show_r2             = TRUE,                # box with CV-R² and selected HDR
  r2_position         = "topleft",           # corner of that box ("topleft", "topright", "bottomleft", "bottomright")
  r2_text_size        = 2.7,                 # text size of both white boxes

  # plot: white box 2 (N)
  show_n              = TRUE,                # box with the number of complete cases
  n_position          = "bottomright",       # corner of that box

  # plot: profiles
  plot_profiles       = "all",               # "all" = five profiles, "constructs" = Effective, Hyper, Hypo only, FALSE = none
  show_local_n_region = TRUE,                # draw the local neighborhoods of the profiles
  label_dir           = c(-1, 0),            # direction of the profile labels as c(x, y); a named list sets single profiles,
                                             # e.g., list(Effective = c(1, 0)) places the Effective label to the right
  label_offset        = 0.040,               # distance of the labels from the profile circles (fraction of the axis range)
  profile_label_size  = 2.0,                 # text size of the profile labels
  profile_digits      = 2,                   # decimals of the MIPV labels
  profile_unstd       = FALSE,               # label the profiles with raw-scale instead of standardized MIPVs

  # plot: other
  show_min_max        = FALSE                # mark the minimum and maximum of the surface within the HDR
)
```

All remaining arguments are documented in the header of `camsq_mipv_analysis.R`.

## Versions

Releases follow semantic versioning. `v1.0.0` is the version used for all analyses in the manuscript.

## License

Copyright (c) 2026 Leon P. Wendt. Licensed under [CC BY-NC-ND 4.0](https://creativecommons.org/licenses/by-nc-nd/4.0/). You may use and redistribute the unmodified function for non-commercial purposes with attribution. Distribution of modified versions and commercial use require written permission of the author.

## Contact

Leon P. Wendt, University of Kassel, [ORCID 0000-0003-2229-2860](https://orcid.org/0000-0003-2229-2860). Please report problems through the Issues tab of this repository.
