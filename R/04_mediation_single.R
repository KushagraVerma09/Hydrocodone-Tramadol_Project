################################################################################
## 15. MEDIATION ANALYSIS
##     X = log ROE   ->   M = 5-HT4 / 5-HT6 related activity   ->   Y = NRS
##
## MODEL
##   a  : X -> M            (exposure changes receptor-related activity)
##   b  : M -> Y | X        (receptor-related activity changes pain,
##                           holding exposure constant)
##   c' : X -> Y | M        (direct effect of exposure not routed through M)
##   c  : X -> Y            (total effect = c' + ab)
##   ab : indirect / mediated effect  <- the quantity of interest
##
## WHY EVERYTHING IS ESTIMATED SEPARATELY BY DRUG AND BY RECEPTOR
##   The c path reverses sign between subgroups (hydrocodone +, tramadol -,
##   interaction p = 0.008). A product ab that is positive in one group and
##   negative in the other averages toward zero when the groups are pooled, so
##   a pooled model would report "no mediation" even if both groups mediated
##   strongly in opposite directions. Formally: drug moderates the mediation.
##   Fitting one model per group is the stratified estimate of that moderated
##   mediation; section 15d tests the moderation itself.
##
## SCALING
##   ROE differs ~40x in scale between the two drugs, so raw-unit a, b and ab
##   are not comparable across groups. All primary estimates are standardized
##   WITHIN GROUP (z-scored on the estimation sample), which puts every path in
##   SD units and makes the hydrocodone/tramadol contrast meaningful.
##   Unstandardized point estimates are reported alongside for reference.
##
## INFERENCE
##   ab is a product of two normal-ish estimates and is not itself normal, so
##   CIs come from a nonparametric bootstrap (percentile method), not from a
##   Sobel test. Each replicate resamples subjects with replacement and refits
##   all three regressions. With n = 17-19 per cell, BCa is unreliable, so
##   percentile intervals are used deliberately.
##
## CAUSAL CAVEAT
##   These are cross-sectional data. Mediation notation implies a causal
##   ordering that the design cannot establish; treat ab as "variance in pain
##   statistically routed through receptor activity", not as a demonstrated
##   mechanism. Sequential-ignorability (no unmeasured M-Y confounding) is
##   untestable here.
##
## Self-contained: needs `master`, `FOCUS_GROUPS`, `OUT_DIR`, `fmt_p` and
## `MED_LABELS` (all from 00_config.R/01_load_data.R) and the packages loaded
## there (dplyr, purrr, tibble, ggplot2, writexl). No new packages required.
##
## Run standalone: Rscript R/04_mediation_single.R
## Sourced by: Run_All.R.
################################################################################

.HTK_PROJECT_DIR <- "/Users/kushagraverma/Work/Projects/Hydrocodone+Tramadol_Project"
if (!exists(".HTK_CONFIG_LOADED")) source(file.path(.HTK_PROJECT_DIR, "R", "00_config.R"))
if (!exists("master"))              source(file.path(.HTK_PROJECT_DIR, "R", "01_load_data.R"))

## ---- 15a. CONFIG -------------------------------------------------------------

MED_X       <- "log_ROE"                      # exposure
MED_M       <- c("R_5HT4", "R_5HT6", "R_MOR")          # mediators, fitted one at a time
MED_Y       <- c("nrs")                       # outcome(s); add "PC2" if wanted
MED_GROUPS  <- FOCUS_GROUPS                   # Hydrocodone group / Tramadol group

## Covariates entered into all three regressions. Leave empty for the primary
## model: with n = 17-19 per cell, each covariate costs a df you cannot spare.
## Sex is the obvious candidate (87% vs 52% female) but is nearly collinear
## with group -- run it as a sensitivity check, not as the headline model.
MED_COVARS  <- character(0)                   # e.g. c("age", "female")

MED_BOOT    <- 5000                           # bootstrap replicates
MED_SEED    <- 42
MED_MIN_N   <- 12                             # skip cells smaller than this

## MED_LABELS / med_lab() and fmt_p() come from 00_config.R -- shared with
## 06_parallel_mediation_sem.R and 07_nested_regression.R so a new/renamed
## receptor label only needs to be added there.

MED_DIR <- file.path(OUT_DIR, "mediation")
dir.create(MED_DIR, showWarnings = FALSE)

set.seed(MED_SEED)

## ---- 15b. CORE MACHINERY -----------------------------------------------------

## z-score the listed columns; returns NULL if any column is constant, which
## can happen in a degenerate bootstrap replicate.
med_scale <- function(df, vars) {
  for (v in vars) {
    s <- stats::sd(df[[v]])
    if (!is.finite(s) || s == 0) return(NULL)
    df[[v]] <- (df[[v]] - mean(df[[v]])) / s
  }
  df
}

## Look up a coefficient by variable name. R strips backticks from coefficient
## names when the variable name is syntactically valid, but keeps them when it
## is not, so a formula written with backticks can yield either form. This
## tries both instead of assuming.
med_pick <- function(nms, v) {
  hit <- match(v, nms)
  if (is.na(hit)) hit <- match(paste0("`", v, "`"), nms)
  if (is.na(hit)) stop("No coefficient found for variable: ", v,
                       "  (available: ", paste(nms, collapse = ", "), ")")
  hit
}

## Fit the three regressions and return every path coefficient.
med_fit_paths <- function(df, xv, mv, yv, covs = character(0)) {
  cov_rhs <- if (length(covs)) paste("+", paste(covs, collapse = " + ")) else ""

  m_a <- lm(as.formula(sprintf("`%s` ~ `%s` %s", mv, xv, cov_rhs)), data = df)
  m_b <- lm(as.formula(sprintf("`%s` ~ `%s` + `%s` %s", yv, mv, xv, cov_rhs)), data = df)
  m_c <- lm(as.formula(sprintf("`%s` ~ `%s` %s", yv, xv, cov_rhs)), data = df)
  
  ca <- summary(m_a)$coefficients
  cb <- summary(m_b)$coefficients
  cc <- summary(m_c)$coefficients
  
  ia <- med_pick(rownames(ca), xv)
  ib <- med_pick(rownames(cb), mv)
  ic <- med_pick(rownames(cb), xv)
  id <- med_pick(rownames(cc), xv)
  
  a  <- ca[ia, 1];  a_se  <- ca[ia, 2];  a_p  <- ca[ia, 4]
  b  <- cb[ib, 1];  b_se  <- cb[ib, 2];  b_p  <- cb[ib, 4]
  cp <- cb[ic, 1];  cp_se <- cb[ic, 2];  cp_p <- cb[ic, 4]
  ct <- cc[id, 1];  c_se  <- cc[id, 2];  c_p  <- cc[id, 4]
  
  ## SE of ab by the multivariate delta method (Sobel). Reported for reference
  ## only -- it assumes ab is normal, which it is not; the bootstrap CI is the
  ## interval to trust and to report.
  ab_se <- sqrt(b^2 * a_se^2 + a^2 * b_se^2)
  
  list(a = a, b = b, cp = cp, c = ct, ab = a * b,
       a_se = a_se, b_se = b_se, cp_se = cp_se, c_se = c_se, ab_se = ab_se,
       a_p = a_p, b_p = b_p, cp_p = cp_p, c_p = c_p,
       r2_m = summary(m_a)$r.squared, r2_y = summary(m_b)$r.squared)
}

## Nonparametric bootstrap of all paths. Re-standardizes inside every replicate
## so the interval reflects the standardized quantity actually being reported.
med_boot <- function(df, xv, mv, yv, covs, R, standardize = TRUE) {
  n   <- nrow(df)
  out <- matrix(NA_real_, nrow = R, ncol = 5,
                dimnames = list(NULL, c("a", "b", "cp", "c", "ab")))
  for (i in seq_len(R)) {
    dd <- df[sample.int(n, n, replace = TRUE), , drop = FALSE]
    if (standardize) {
      dd <- med_scale(dd, c(xv, mv, yv))
      if (is.null(dd)) next
    }
    res <- tryCatch(med_fit_paths(dd, xv, mv, yv, covs), error = function(e) NULL)
    if (!is.null(res)) out[i, ] <- c(res$a, res$b, res$cp, res$c, res$ab)
  }
  out[stats::complete.cases(out), , drop = FALSE]
}

## Percentile CI plus a two-sided bootstrap p (proportion of replicates on the
## far side of zero, doubled).
med_ci <- function(v, conf = 0.95) {
  v <- v[is.finite(v)]
  if (!length(v)) return(c(lo = NA, hi = NA, p = NA))
  q <- stats::quantile(v, c((1 - conf) / 2, 1 - (1 - conf) / 2), names = FALSE)
  p <- 2 * min(mean(v <= 0), mean(v >= 0))
  c(lo = q[1], hi = q[2], p = min(p, 1))
}


## ---- 15c. SINGLE-MEDIATOR MODELS, ONE PER GROUP x RECEPTOR -------------------

med_rows   <- list()
med_bootd  <- list()
med_data   <- list()

for (yv in MED_Y) {
  for (g in MED_GROUPS) {
    for (mv in MED_M) {
      
      cell <- sprintf("%s | %s -> %s -> %s", g, MED_X, mv, yv)
      keep <- c(MED_X, mv, yv, MED_COVARS)
      
      d <- master[master$Plot_Group == g, , drop = FALSE]
      d <- d[stats::complete.cases(d[, keep, drop = FALSE]), keep, drop = FALSE]
      d <- d[apply(d, 1, function(r) all(is.finite(r))), , drop = FALSE]
      
      if (nrow(d) < MED_MIN_N) {
        message(sprintf("Skipping %s -- n = %d < MED_MIN_N", cell, nrow(d)))
        next
      }
      
      ## point estimates, unstandardized and standardized
      raw <- med_fit_paths(d, MED_X, mv, yv, MED_COVARS)
      ds  <- med_scale(d, c(MED_X, mv, yv))
      std <- med_fit_paths(ds, MED_X, mv, yv, MED_COVARS)
      
      ## bootstrap the standardized solution
      bd <- med_boot(d, MED_X, mv, yv, MED_COVARS, MED_BOOT, standardize = TRUE)
      ci_ab <- med_ci(bd[, "ab"]);  ci_a  <- med_ci(bd[, "a"])
      ci_b  <- med_ci(bd[, "b"]);   ci_cp <- med_ci(bd[, "cp"])
      ci_c  <- med_ci(bd[, "c"])
      
      med_rows[[cell]] <- tibble(
        outcome   = yv,
        group     = g,
        mediator  = mv,
        n         = nrow(d),
        ## standardized (primary). SE = model standard error;
        ## lo/hi = bootstrap percentile 95% CI.
        a_std     = std$a,  a_se  = std$a_se,  a_p   = std$a_p,
        a_lo      = ci_a[["lo"]],  a_hi  = ci_a[["hi"]],
        b_std     = std$b,  b_se  = std$b_se,  b_p   = std$b_p,
        b_lo      = ci_b[["lo"]],  b_hi  = ci_b[["hi"]],
        cprime_std= std$cp, cp_se = std$cp_se, cp_p  = std$cp_p,
        cp_lo     = ci_cp[["lo"]], cp_hi = ci_cp[["hi"]],
        c_total_std = std$c, c_se = std$c_se,  c_p   = std$c_p,
        c_lo      = ci_c[["lo"]],  c_hi  = ci_c[["hi"]],
        ab_std    = std$ab,
        ab_se_delta = std$ab_se,
        ab_lo     = ci_ab[["lo"]],
        ab_hi     = ci_ab[["hi"]],
        ab_boot_p = ci_ab[["p"]],
        ab_sig    = !is.na(ci_ab[["lo"]]) & (ci_ab[["lo"]] * ci_ab[["hi"]] > 0),
        prop_med  = ifelse(abs(std$c) > 0.05, std$ab / std$c, NA_real_),
        ## unstandardized (reference)
        a_raw = raw$a, b_raw = raw$b, cprime_raw = raw$cp,
        c_total_raw = raw$c, ab_raw = raw$ab,
        R2_mediator = std$r2_m, R2_outcome = std$r2_y,
        n_boot_ok = nrow(bd)
      )
      
      med_bootd[[cell]] <- bd
      med_data[[cell]]  <- cbind(d, group = g, mediator = mv, outcome = yv)
      
      cat(sprintf("\n===== MEDIATION: %s =====\n", cell))
      cat(sprintf("n = %d\n", nrow(d)))
      cat(sprintf("  a  (X->M)      = %7.3f  SE %.3f  p = %s\n",
                  std$a,  std$a_se,  fmt_p(std$a_p)))
      cat(sprintf("  b  (M->Y|X)    = %7.3f  SE %.3f  p = %s\n",
                  std$b,  std$b_se,  fmt_p(std$b_p)))
      cat(sprintf("  c' (X->Y|M)    = %7.3f  SE %.3f  p = %s   <- DIRECT effect\n",
                  std$cp, std$cp_se, fmt_p(std$cp_p)))
      cat(sprintf("  c  (X->Y)      = %7.3f  SE %.3f  p = %s   <- TOTAL effect\n",
                  std$c,  std$c_se,  fmt_p(std$c_p)))
      cat(sprintf("  ab (indirect)  = %7.3f  95%% CI [%.3f, %.3f]  %s\n",
                  std$ab, ci_ab[["lo"]], ci_ab[["hi"]],
                  ifelse(ci_ab[["lo"]] * ci_ab[["hi"]] > 0,
                         "*** CI excludes zero", "CI includes zero")))
      cat(sprintf("  check: c' + ab = %.3f  (should equal c = %.3f)\n",
                  std$cp + std$ab, std$c))
    }
  }
}

med_tbl <- bind_rows(med_rows)


## ---- 15d. MODERATED MEDIATION: does the mediation differ by drug? ------------
##
## The stratified models in 15c show two separate pictures but never test
## whether they differ. This does. The index of moderated mediation is
##      ab(hydrocodone) - ab(tramadol)
## bootstrapped with resampling STRATIFIED WITHIN GROUP so both cell sizes are
## preserved in every replicate. Path differences (a, b, c') are bootstrapped
## the same way. Standardization is done within group inside each replicate,
## for the scale reason given in the header.

modmed_rows <- list()

for (yv in MED_Y) {
  for (mv in MED_M) {
    
    keep <- c(MED_X, mv, yv, MED_COVARS)
    dg <- lapply(MED_GROUPS, function(g) {
      d <- master[master$Plot_Group == g, , drop = FALSE]
      d <- d[stats::complete.cases(d[, keep, drop = FALSE]), keep, drop = FALSE]
      d[apply(d, 1, function(r) all(is.finite(r))), , drop = FALSE]
    })
    names(dg) <- MED_GROUPS
    if (any(vapply(dg, nrow, 1L) < MED_MIN_N)) next
    
    one_rep <- function(resample) {
      est <- lapply(dg, function(d) {
        dd <- if (resample) d[sample.int(nrow(d), nrow(d), replace = TRUE), , drop = FALSE] else d
        dd <- med_scale(dd, c(MED_X, mv, yv))
        if (is.null(dd)) return(NULL)
        tryCatch(med_fit_paths(dd, MED_X, mv, yv, MED_COVARS), error = function(e) NULL)
      })
      if (any(vapply(est, is.null, TRUE))) return(rep(NA_real_, 4))
      c(ab = est[[1]]$ab - est[[2]]$ab,
        a  = est[[1]]$a  - est[[2]]$a,
        b  = est[[1]]$b  - est[[2]]$b,
        cp = est[[1]]$cp - est[[2]]$cp)
    }
    
    pt <- one_rep(FALSE)
    names(pt) <- c("ab", "a", "b", "cp")
    bs <- t(vapply(seq_len(MED_BOOT), function(i) one_rep(TRUE), numeric(4)))
    colnames(bs) <- c("ab", "a", "b", "cp")
    bs <- bs[stats::complete.cases(bs), , drop = FALSE]
    
    ci <- lapply(c("ab", "a", "b", "cp"), function(k) med_ci(bs[, k]))
    names(ci) <- c("ab", "a", "b", "cp")
    
    modmed_rows[[paste(yv, mv)]] <- tibble(
      outcome = yv, mediator = mv,
      contrast = paste(MED_GROUPS[1], "-", MED_GROUPS[2]),
      n_1 = nrow(dg[[1]]), n_2 = nrow(dg[[2]]),
      index_modmed = pt[["ab"]],
      index_lo = ci$ab[["lo"]], index_hi = ci$ab[["hi"]], index_p = ci$ab[["p"]],
      a_diff = pt[["a"]],  a_lo  = ci$a[["lo"]],  a_hi  = ci$a[["hi"]],  a_p  = ci$a[["p"]],
      b_diff = pt[["b"]],  b_lo  = ci$b[["lo"]],  b_hi  = ci$b[["hi"]],  b_p  = ci$b[["p"]],
      cp_diff = pt[["cp"]], cp_lo = ci$cp[["lo"]], cp_hi = ci$cp[["hi"]], cp_p = ci$cp[["p"]],
      n_boot_ok = nrow(bs)
    )
  }
}

modmed_tbl <- bind_rows(modmed_rows)

if (nrow(modmed_tbl)) {
  cat("\n===== INDEX OF MODERATED MEDIATION (bootstrap) =====\n")
  print(as.data.frame(modmed_tbl %>%
                        select(outcome, mediator, index_modmed, index_lo, index_hi, index_p,
                               b_diff, b_p, cp_diff, cp_p)),
        row.names = FALSE, digits = 3)
}


## ---- 15e. PARALLEL MEDIATORS: 5-HT4 AND 5-HT6 IN ONE MODEL -------------------
##
## The PI's instruction is to keep the receptors separate, and 15c does that.
## This block is the supplement that answers the obvious follow-up question:
## the two receptors are correlated (r ~ 0.30 hydrocodone, ~0.44 tramadol), so
## do they carry the SAME variance in pain or different variance? Both mediators
## enter the b-model together; each b is then adjusted for the other receptor.
## Report this as a robustness check, not as the primary analysis.

par_rows <- list()

for (yv in MED_Y) {
  for (g in MED_GROUPS) {
    
    keep <- c(MED_X, MED_M, yv, MED_COVARS)
    d <- master[master$Plot_Group == g, , drop = FALSE]
    d <- d[stats::complete.cases(d[, keep, drop = FALSE]), keep, drop = FALSE]
    d <- d[apply(d, 1, function(r) all(is.finite(r))), , drop = FALSE]
    if (nrow(d) < MED_MIN_N) next
    
    par_paths <- function(dd) {
      dd <- med_scale(dd, c(MED_X, MED_M, yv))
      if (is.null(dd)) return(NULL)
      cov_rhs <- if (length(MED_COVARS)) paste("+", paste(MED_COVARS, collapse = " + ")) else ""
      as_ <- vapply(MED_M, function(mv) {
        ca <- coef(lm(as.formula(sprintf("`%s` ~ `%s` %s", mv, MED_X, cov_rhs)),
                      data = dd))
        unname(ca[med_pick(names(ca), MED_X)])
      }, numeric(1))
      
      mb <- lm(as.formula(sprintf("`%s` ~ %s + `%s` %s", yv,
                                  paste0("`", MED_M, "`", collapse = " + "),
                                  MED_X, cov_rhs)), data = dd)
      cb  <- coef(mb)
      bs_ <- vapply(MED_M, function(mv) unname(cb[med_pick(names(cb), mv)]),
                    numeric(1))
      
      c(setNames(as_ * bs_, paste0("ab_", MED_M)),
        setNames(as_,       paste0("a_",  MED_M)),
        setNames(bs_,       paste0("b_",  MED_M)),
        cp       = unname(cb[med_pick(names(cb), MED_X)]),
        ab_total = sum(as_ * bs_))
    }
    
    pt <- par_paths(d)
    if (is.null(pt)) next
    bs <- t(vapply(seq_len(MED_BOOT), function(i) {
      r <- tryCatch(par_paths(d[sample.int(nrow(d), nrow(d), replace = TRUE), , drop = FALSE]),
                    error = function(e) NULL)
      if (is.null(r)) rep(NA_real_, length(pt)) else r
    }, numeric(length(pt))))
    colnames(bs) <- names(pt)
    bs <- bs[stats::complete.cases(bs), , drop = FALSE]
    
    par_rows[[paste(yv, g)]] <- bind_rows(lapply(names(pt), function(k) {
      ci <- med_ci(bs[, k])
      tibble(outcome = yv, group = g, n = nrow(d), quantity = k,
             estimate = unname(pt[[k]]),
             lo = ci[["lo"]], hi = ci[["hi"]], boot_p = ci[["p"]])
    }))
  }
}

parallel_tbl <- bind_rows(par_rows)


## ---- 15f. DIAGNOSTICS --------------------------------------------------------

diag_tbl <- bind_rows(lapply(MED_GROUPS, function(g) {
  keep <- c(MED_X, MED_M, MED_Y)
  d <- master[master$Plot_Group == g, , drop = FALSE]
  d <- d[stats::complete.cases(d[, keep, drop = FALSE]), keep, drop = FALSE]
  if (!nrow(d)) return(NULL)
  cmb <- utils::combn(keep, 2, simplify = FALSE)
  bind_rows(lapply(cmb, function(p) tibble(
    group = g, n = nrow(d), var_1 = p[1], var_2 = p[2],
    pearson_r = stats::cor(d[[p[1]]], d[[p[2]]]),
    p_value   = stats::cor.test(d[[p[1]]], d[[p[2]]])$p.value)))
}))


## ---- 15g. PATH DIAGRAMS ------------------------------------------------------
## One diagram per cell. Solid arrow = p < .05, dashed = not significant.
## The indirect effect and its bootstrap CI sit under the title.

## Every path is annotated as  coefficient (SE), p.
## c' is the DIRECT effect; the TOTAL effect c is printed under the X->Y arrow
## so the decomposition c = c' + ab is readable off the figure.
med_path_diagram <- function(row) {
  sig <- function(p) if (!is.na(p) && p < 0.05) "solid" else "22"
  
  nodes <- tibble(
    x = c(0.10, 0.50, 0.90), y = c(0.28, 0.82, 0.28),
    lab = c(med_lab(MED_X), med_lab(row$mediator), med_lab(row$outcome)))
  
  edges <- tibble(
    x    = c(0.16, 0.56, 0.20), y    = c(0.35, 0.76, 0.28),
    xend = c(0.44, 0.84, 0.80), yend = c(0.76, 0.35, 0.28),
    lty  = c(sig(row$a_p), sig(row$b_p), sig(row$cp_p)),
    lab  = c(sprintf("a = %.2f (SE %.2f)\np = %s",
                     row$a_std,  row$a_se,  fmt_p(row$a_p)),
             sprintf("b = %.2f (SE %.2f)\np = %s",
                     row$b_std,  row$b_se,  fmt_p(row$b_p)),
             sprintf("c' = %.2f (SE %.2f), p = %s   [DIRECT]",
                     row$cprime_std, row$cp_se, fmt_p(row$cp_p))),
    lx = c(0.25, 0.75, 0.50), ly = c(0.60, 0.60, 0.205))
  
  ## total effect, printed below the direct effect
  tot <- sprintf("c  = %.2f (SE %.2f), p = %s   [TOTAL = c' + ab]",
                 row$c_total_std, row$c_se, fmt_p(row$c_p))
  
  ggplot() +
    geom_segment(data = edges,
                 aes(x = x, y = y, xend = xend, yend = yend, linetype = I(lty)),
                 linewidth = 0.8, colour = "grey25",
                 arrow = arrow(length = unit(3.2, "mm"), type = "closed")) +
    geom_label(data = nodes, aes(x, y, label = lab),
               size = 4.6, fontface = "bold", label.padding = unit(4, "mm"),
               fill = "white", colour = "black") +
    geom_text(data = edges, aes(lx, ly, label = lab),
              size = 3.8, lineheight = 0.95, colour = "grey15") +
    annotate("text", x = 0.50, y = 0.135, label = tot,
             size = 3.8, colour = "grey15") +
    labs(title = sprintf("%s  --  %s", row$group, med_lab(row$mediator)),
         subtitle = sprintf(
           "n = %d   |   INDIRECT ab = %.3f, bootstrap 95%% CI [%.3f, %.3f]%s",
           row$n, row$ab_std, row$ab_lo, row$ab_hi,
           ifelse(row$ab_sig, "  *", "  (contains zero)")),
         caption = paste0(
           "Standardized paths (within-group z-scores). Solid arrow = p < .05.\n",
           "SE = model standard error. The ab interval is a bootstrap percentile CI, ",
           "not SE-based:\nab is a product of two estimates and is not normally distributed.")) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0.08, 1), expand = FALSE) +
    theme_void(base_size = 13) +
    theme(plot.title    = element_text(face = "bold", hjust = 0, margin = margin(b = 2)),
          plot.subtitle = element_text(hjust = 0, margin = margin(b = 10)),
          plot.caption  = element_text(hjust = 0, colour = "grey40", size = 8),
          plot.margin   = margin(12, 12, 10, 12))
}

for (i in seq_len(nrow(med_tbl))) {
  row <- med_tbl[i, ]
  p <- med_path_diagram(row)
  ggsave(file.path(MED_DIR, sprintf("Mediation_path_%s_%s_%s.png",
                                    row$outcome, gsub(" ", "", row$group), row$mediator)),
         p, width = 7.5, height = 5.2, dpi = 300)
}

## Effect decomposition: total, direct and indirect side by side, each with its
## bootstrap CI. This is the figure that answers "where did the effect go" --
## when c' sits on top of c and ab sits on zero, nothing is routed through M.
if (nrow(med_tbl)) {
  
  decomp <- bind_rows(
    med_tbl %>% transmute(group, mediator, outcome, n,
                          effect = "Total (c)",     est = c_total_std,
                          lo = c_lo,  hi = c_hi),
    med_tbl %>% transmute(group, mediator, outcome, n,
                          effect = "Direct (c')",   est = cprime_std,
                          lo = cp_lo, hi = cp_hi),
    med_tbl %>% transmute(group, mediator, outcome, n,
                          effect = "Indirect (ab)", est = ab_std,
                          lo = ab_lo, hi = ab_hi)
  ) %>%
    mutate(
      effect = factor(effect, levels = c("Indirect (ab)", "Direct (c')", "Total (c)")),
      cell   = paste0(group, "\n", vapply(mediator, med_lab, ""),
                      "  ->  ", vapply(outcome, med_lab, ""), "   (n = ", n, ")")
    )
  
  fp <- ggplot(decomp, aes(x = est, y = effect, colour = group)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.16, linewidth = 0.8) +
    geom_point(size = 3) +
    facet_wrap(~ cell, ncol = 2) +
    scale_colour_manual(values = if (exists("custom_colors")) custom_colors else NULL,
                        guide = "none") +
    labs(title = "Effect decomposition: total = direct + indirect",
         subtitle = sprintf(
           "%s -> receptor -> outcome. Standardized within group; bootstrap percentile 95%% CIs.",
           med_lab(MED_X)),
         x = "Effect (SD units)", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold", size = 9, lineheight = 1.1),
          panel.spacing = unit(1.1, "lines"))
  
  ggsave(file.path(MED_DIR, "Mediation_effect_decomposition.png"),
         fp, width = 10, height = 7, dpi = 300)
  
  ## Indirect effects alone, all cells on one axis (the original forest plot).
  fp2 <- med_tbl %>%
    mutate(cell = paste0(group, "\n", vapply(mediator, med_lab, ""),
                         "  ->  ", vapply(outcome, med_lab, ""))) %>%
    ggplot(aes(x = ab_std, y = cell, colour = group)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmin = ab_lo, xmax = ab_hi), height = 0.18, linewidth = 0.8) +
    geom_point(size = 3.2) +
    scale_colour_manual(values = if (exists("custom_colors")) custom_colors else NULL) +
    labs(title = "Indirect effects (ab) with bootstrap 95% CIs",
         subtitle = sprintf("%s -> receptor -> outcome, standardized within group",
                            med_lab(MED_X)),
         x = "Indirect effect (SD units)", y = NULL, colour = NULL) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
  
  ggsave(file.path(MED_DIR, "Mediation_indirect_effects_forest.png"),
         fp2, width = 9, height = 5.5, dpi = 300)
}


## ---- 15h. EXPORT -------------------------------------------------------------

write_xlsx(
  list(
    Single_mediator   = med_tbl,
    Moderated_med     = modmed_tbl,
    Parallel_mediator = parallel_tbl,
    Diagnostics       = diag_tbl,
    Data_used         = bind_rows(med_data)
  ),
  file.path(MED_DIR, "Stats_Mediation_ROE_receptor_NRS.xlsx")
)

message("\nMediation results, path diagrams and forest plot written to: ",
        normalizePath(MED_DIR))


## ---- 15i. OPTIONAL SENSITIVITY CHECKS ----------------------------------------
## Uncomment individually. Each one re-runs section 15 with a different setting;
## easiest is to change the config at 15a and re-source from there.
##
## 1) Adjust for sex (the main confound: 87% vs 52% female):
##      MED_COVARS <- c("female")
##
## 2) Use prescribed rather than measured exposure:
##      MED_X <- "log_MME"        # larger n, but MME is what was prescribed
##
## 3) Add pain quality as a second outcome:
##      MED_Y <- c("nrs", "PC2")
##
## 4) Restrict to the overlapping ROE range (addresses the limited-overlap
##    caveat raised in the lab meeting):
##      rng <- range(intersect(
##        range(master$log_ROE[master$Plot_Group == MED_GROUPS[1]], na.rm = TRUE),
##        range(master$log_ROE[master$Plot_Group == MED_GROUPS[2]], na.rm = TRUE)))
##      master_overlap <- master[is.na(master$log_ROE) |
##                               (master$log_ROE >= rng[1] & master$log_ROE <= rng[2]), ]
##    then swap `master` for `master_overlap` in 15c-15f.
################################################################################
