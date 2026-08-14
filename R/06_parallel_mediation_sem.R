################################################################################
## 17. PARALLEL MULTIPLE-MEDIATOR MODEL (lavaan):
##     X = log ROE  ->  {several receptors, simultaneously}  ->  Y = NRS
##
## WHAT THIS ADDS OVER SECTION 15
##   Section 15 fits one mediator at a time, so each receptor's indirect effect
##   is estimated ignoring the others. Receptor maps are correlated, so those
##   single-mediator ab's double-count shared variance. Here all mediators enter
##   the outcome equation together: b_k is receptor k's unique contribution
##   holding the other mediators (and ROE) constant, and the specific indirect
##   effects a_k*b_k are therefore comparable to each other and can be
##   differenced (ind_R1 - ind_R2) as a formal test of which receptor carries
##   more of the ROE-pain path.
##
## WHY THE EARLIER VERSION ERRORED WITH "no variance in ROE"
##   Two separate causes, both fixed below:
##     1. log_ROE = log10(ROE) and ROE == 0 for 4 subjects in each group, so
##        log_ROE contains -Inf. na.omit() does NOT drop -Inf, lavaan then
##        computes a non-finite sample variance for that column and aborts with
##        the "no variance" error. Sections 15 and 16 never hit this because
##        they filter rows with an explicit is.finite() test -- 17 now does the
##        same (see pmed_prep).
##     2. Raw ROE is on the order of 1e-4 (tramadol) to 1e-2 (hydrocodone);
##        its variance (2e-8 in tramadol) is numerically indistinguishable from
##        zero next to NRS (var ~ 5) and trips the same check even when finite.
##        Everything is z-scored within cell before fitting, which removes the
##        scaling problem entirely and makes the paths standardized -- required
##        anyway for the ind_Rj - ind_Rk contrasts to mean anything, since the
##        receptors have different raw SDs.
##
## MODEL, FIT PER GROUP
##   M_k  = a_k * X + e_k                      (k = 1..K mediators)
##   Y    = cprime * X + sum_k b_k * M_k + e_y
##   mediator residuals left free to covary (they are correlated maps)
##   ind_k := a_k * b_k ; ind_total := sum ind_k ; total := cprime + ind_total
##
## FITS PRODUCED
##   17d  one model per group (the primary result -- the ROE-pain slope
##        reverses sign between drugs, so a pooled model is not interpretable)
##   17e  pooled model with a group dummy as covariate, for reference
##   17f  two-group model with group differences in each specific indirect
##        effect -- the parallel-mediator analogue of the index of moderated
##        mediation in 15d
##
## SAMPLE SIZE -- READ THIS BEFORE INTERPRETING
##   After dropping ROE == 0 rows there are ~19 (hydrocodone) and ~17
##   (tramadol) complete cases. The outcome equation with K = 4 mediators
##   spends 5 predictors on that. This is a descriptive decomposition, not a
##   powered test; the CIs will be wide and should be reported as such. Keep
##   PMED_M short and pharmacologically motivated rather than throwing the
##   19-receptor panel at it -- with K > 4 the model is not identified in any
##   useful sense at this n.
##
## Self-contained: needs `master`, `FOCUS_GROUPS`, `OUT_DIR`, `MED_LABELS`
## (all from 00_config.R/01_load_data.R), plus lavaan and writexl.
##
## Run standalone: Rscript R/06_parallel_mediation_sem.R
## Sourced by: Run_All.R.
################################################################################

.HTK_PROJECT_DIR <- "/Users/kushagraverma/Work/Projects/Hydrocodone+Tramadol_Project"
## Always re-sourced -- see the note above this line in 01_load_data.R.
source(file.path(.HTK_PROJECT_DIR, "R", "00_config.R"))
if (!exists("master"))              source(file.path(.HTK_PROJECT_DIR, "R", "01_load_data.R"))

library(lavaan)

## ---- 17a. CONFIG -------------------------------------------------------------

PMED_X      <- "log_ROE"                    # exposure
PMED_M      <- c("R_MOR", "R_5HT1A", "R_5HT4", "R_5HT6")   # parallel mediators
PMED_Y      <- "nrs"                        # outcome
PMED_GROUPS <- FOCUS_GROUPS                 # Hydrocodone group / Tramadol group
PMED_COVARS <- character(0)                 # e.g. c("age", "female")

PMED_BOOT   <- 5000                         # bootstrap draws (10000 if you can wait)
PMED_SEED   <- 42
PMED_MIN_N  <- 12                           # skip a group below this
PMED_CI     <- "perc"                       # percentile CI; BCa is unreliable at this n

PMED_DIR <- file.path(OUT_DIR, "lavaan_mediation")
dir.create(PMED_DIR, showWarnings = FALSE, recursive = TRUE)

set.seed(PMED_SEED)

## NOTE: the script-level fmt_p (section 2) takes a SINGLE p and returns "--"
## for NA -- calling it on a column raises
##   "the condition has length > 1".
## This section formats whole columns at a time, so it uses its own vectorized
## formatter rather than fmt_p. Keep them separate; do not "simplify" this away.
pmed_fmt_p <- function(p) {
  ifelse(is.na(p), "--",
  ifelse(p < 1e-5, "<1e-5",
  ifelse(p < 0.001, sprintf("%.1e", p), sprintf("%.3f", p))))
}

## MED_LABELS comes from 00_config.R.
pmed_lab <- function(v) {
  hit <- MED_LABELS[v]
  if (length(hit) == 1 && !is.na(hit)) unname(hit) else v
}


## ---- 17b. DATA PREP ----------------------------------------------------------

## Select the modelled columns, keep only rows that are complete AND finite
## (this is what kills the -Inf from log10(0)), then z-score every variable.
## Returns NULL with an explanatory message instead of letting lavaan fail with
## its opaque "no variance" error.
pmed_prep <- function(df, xv, mvs, yv, covs = character(0), label = "") {

  keep <- c(xv, mvs, yv, covs)
  miss <- setdiff(keep, names(df))
  if (length(miss)) {
    message(sprintf("[%s] missing column(s): %s -- skipped",
                    label, paste(miss, collapse = ", ")))
    return(NULL)
  }

  d  <- df[, keep, drop = FALSE]
  n0 <- nrow(d)
  ok <- stats::complete.cases(d) &
        apply(d, 1, function(r) all(is.finite(as.numeric(r))))
  d  <- d[ok, , drop = FALSE]

  n_nonfinite <- sum(stats::complete.cases(df[, keep, drop = FALSE])) - nrow(d)
  if (n_nonfinite > 0)
    message(sprintf("[%s] dropped %d row(s) with non-finite values (e.g. log10(ROE) where ROE == 0)",
                    label, n_nonfinite))
  if (n0 - nrow(d) - n_nonfinite > 0)
    message(sprintf("[%s] dropped %d row(s) with missing values",
                    label, n0 - nrow(d) - n_nonfinite))

  if (nrow(d) < PMED_MIN_N) {
    message(sprintf("[%s] n = %d < PMED_MIN_N (%d) -- skipped",
                    label, nrow(d), PMED_MIN_N))
    return(NULL)
  }

  ## variance diagnostic: report before standardizing so a degenerate variable
  ## is named explicitly rather than surfacing as a lavaan error
  sds <- vapply(d, stats::sd, numeric(1))
  dead <- names(sds)[!is.finite(sds) | sds < 1e-12]
  if (length(dead)) {
    message(sprintf("[%s] no usable variance in: %s (sd = %s) -- skipped",
                    label, paste(dead, collapse = ", "),
                    paste(signif(sds[dead], 3), collapse = ", ")))
    return(NULL)
  }

  ## rename to syntactically safe model names, z-score, keep the mapping
  new <- c("X", paste0("M", seq_along(mvs)), "Y",
           if (length(covs)) paste0("C", seq_along(covs)) else character(0))
  out <- as.data.frame(scale(d))
  names(out) <- new
  attr(out, "map") <- setNames(keep, new)
  attr(out, "sd_raw") <- sds
  out
}

## Build the lavaan syntax for K parallel mediators. Generated rather than
## hard-coded so PMED_M can be any length without editing the model string.
pmed_syntax <- function(K, ncov = 0) {

  cov_rhs <- if (ncov) paste("", paste0("+ d", seq_len(ncov), "*C", seq_len(ncov)),
                             collapse = " ") else ""

  a_paths <- paste0("  M", 1:K, " ~ a", 1:K, "*X",
                    if (ncov) paste0(cov_rhs) else "", collapse = "\n")

  y_path  <- paste0("  Y ~ cprime*X + ",
                    paste0("b", 1:K, "*M", 1:K, collapse = " + "),
                    if (ncov) cov_rhs else "")

  ## free residual covariances among mediators
  cors <- if (K > 1) {
    pr <- utils::combn(K, 2)
    paste0("  M", pr[1, ], " ~~ M", pr[2, ], collapse = "\n")
  } else ""

  inds <- paste0("  ind_M", 1:K, " := a", 1:K, "*b", 1:K, collapse = "\n")

  tot  <- paste0("  ind_total := ", paste0("ind_M", 1:K, collapse = " + "), "\n",
                 "  total := cprime + ind_total")

  ## pairwise contrasts between specific indirect effects
  contr <- if (K > 1) {
    pr <- utils::combn(K, 2)
    paste0("  diff_M", pr[1, ], "_M", pr[2, ],
           " := ind_M", pr[1, ], " - ind_M", pr[2, ], collapse = "\n")
  } else ""

  paste(a_paths, "", y_path, "", cors, "", inds, "", tot, "", contr, sep = "\n")
}

## Tidy the := rows and the labelled paths out of a fit, with readable names.
pmed_tidy <- function(fit, map, mvs, group_lab = NA_character_) {

  pe <- lavaan::parameterEstimates(fit, standardized = TRUE, ci = TRUE,
                                   boot.ci.type = PMED_CI)
  pe <- pe[pe$op %in% c("~", ":=") | (pe$op == "~~" & pe$lhs != pe$rhs), ]

  pretty <- function(s) {
    for (i in seq_along(mvs)) s <- gsub(paste0("M", i), pmed_lab(mvs[i]), s, fixed = TRUE)
    s <- gsub("^X$", pmed_lab(map[["X"]]), s)
    s <- gsub("^Y$", pmed_lab(map[["Y"]]), s)
    s
  }

  data.frame(
    group    = group_lab,
    term     = vapply(pe$lhs, pretty, ""),
    op       = pe$op,
    rhs      = vapply(pe$rhs, pretty, ""),
    label    = pe$label,
    est      = pe$est,
    se       = pe$se,
    z        = pe$z,
    p        = pe$pvalue,
    ci_lo    = pe$ci.lower,
    ci_hi    = pe$ci.upper,
    std_all  = pe$std.all,
    row.names = NULL, stringsAsFactors = FALSE
  )
}

## Fit one parallel-mediation model with bootstrap SEs, surviving the case
## where some bootstrap replicates fail to converge.
pmed_fit <- function(d, K, ncov = 0, label = "", group_var = NULL) {
  mod <- pmed_syntax(K, ncov)
  fit <- tryCatch(
    lavaan::sem(mod, data = d, se = "bootstrap", bootstrap = PMED_BOOT,
                group = group_var, fixed.x = FALSE,
                estimator = "ML", missing = "listwise"),
    error = function(e) { message(sprintf("[%s] lavaan failed: %s", label, conditionMessage(e))); NULL })
  fit
}


## ---- 17c. THE MODEL SYNTAX ACTUALLY BEING FITTED -----------------------------

cat("\n===== PARALLEL MEDIATION MODEL SYNTAX (K =", length(PMED_M), "mediators) =====\n")
cat("mediators: ", paste(sprintf("M%d = %s", seq_along(PMED_M), PMED_M), collapse = ",  "), "\n")
cat(pmed_syntax(length(PMED_M), length(PMED_COVARS)), "\n")


## ---- 17d. PRIMARY: ONE MODEL PER GROUP ---------------------------------------

pmed_rows <- list()
pmed_fits <- list()
pmed_par  <- list()   # raw labelled estimates per group, used by the figures
pmed_n    <- c()

for (g in PMED_GROUPS) {

  lab <- sprintf("17d %s", g)
  dg  <- pmed_prep(master[master$Plot_Group == g, , drop = FALSE],
                   PMED_X, PMED_M, PMED_Y, PMED_COVARS, label = lab)
  if (is.null(dg)) next

  message(sprintf("[%s] fitting on n = %d, %d bootstrap draws ...",
                  lab, nrow(dg), PMED_BOOT))

  fit <- pmed_fit(dg, length(PMED_M), length(PMED_COVARS), label = lab)
  if (is.null(fit)) next

  pmed_fits[[g]] <- fit

  cat("\n=====", g, "-- parallel mediation, n =", nrow(dg),
      "(all variables z-scored within group) =====\n")
  print(summary(fit, standardized = TRUE, ci = TRUE, rsquare = TRUE))

  tid <- pmed_tidy(fit, attr(dg, "map"), PMED_M, group_lab = g)
  tid$n <- nrow(dg)
  pmed_rows[[g]] <- tid

  ## keep the untouched labelled estimates -- the figures index them by the
  ## lavaan labels (a1..aK, b1..bK, cprime, ind_Mk, ...), not by pretty names
  pmed_par[[g]] <- lavaan::parameterEstimates(fit, standardized = TRUE, ci = TRUE,
                                              boot.ci.type = PMED_CI)
  pmed_n[g] <- nrow(dg)
}

pmed_tbl <- if (length(pmed_rows)) do.call(rbind, pmed_rows) else NULL

if (!is.null(pmed_tbl)) {
  cat("\n===== SPECIFIC INDIRECT EFFECTS AND CONTRASTS (standardized) =====\n")
  show <- pmed_tbl[pmed_tbl$op == ":=", ]
  show$p <- pmed_fmt_p(show$p)
  print(show[, c("group", "n", "term", "est", "se", "ci_lo", "ci_hi", "p")],
        row.names = FALSE, digits = 3)
}


## ---- 17e. REFERENCE: POOLED MODEL WITH GROUP AS COVARIATE --------------------
## Reported for completeness only. If the ROE-pain slope reverses sign between
## drugs (as sections 15/16 indicate), the pooled cprime is an average of two
## opposite effects and means little on its own.

pmed_pooled_tbl <- NULL
pool_src <- master[master$Plot_Group %in% PMED_GROUPS, , drop = FALSE]
pool_src$grp_dummy <- as.numeric(pool_src$Plot_Group == PMED_GROUPS[2])

d_pool <- pmed_prep(pool_src, PMED_X, PMED_M, PMED_Y,
                    c(PMED_COVARS, "grp_dummy"), label = "17e pooled")

if (!is.null(d_pool)) {
  fit_pool <- pmed_fit(d_pool, length(PMED_M), length(PMED_COVARS) + 1,
                       label = "17e pooled")
  if (!is.null(fit_pool)) {
    cat("\n===== POOLED (both groups, group entered as covariate), n =",
        nrow(d_pool), "=====\n")
    print(summary(fit_pool, standardized = TRUE, ci = TRUE, rsquare = TRUE))
    pmed_pooled_tbl <- pmed_tidy(fit_pool, attr(d_pool, "map"), PMED_M,
                                 group_lab = "Pooled")
    pmed_pooled_tbl$n <- nrow(d_pool)
    pmed_fits[["Pooled"]] <- fit_pool
  }
}


## ---- 17f. GROUP DIFFERENCES IN EACH SPECIFIC INDIRECT EFFECT -----------------
## Two-group model, all paths free, with the between-group difference in each
## specific indirect effect defined and bootstrapped. This is the parallel-
## mediator analogue of section 15d's index of moderated mediation: a
## significant diff_* means that receptor carries a different amount of the
## ROE-pain path in one drug than the other. Note each group is standardized
## separately in 17d but jointly here, so 17f estimates are on a common scale
## and will not equal the 17d numbers exactly.

pmed_modmed_tbl <- NULL

d_mg <- pmed_prep(pool_src, PMED_X, PMED_M, PMED_Y, PMED_COVARS,
                  label = "17f multigroup")

if (!is.null(d_mg)) {
  d_mg$grp <- pool_src$Plot_Group[
    stats::complete.cases(pool_src[, c(PMED_X, PMED_M, PMED_Y, PMED_COVARS)]) &
    apply(pool_src[, c(PMED_X, PMED_M, PMED_Y, PMED_COVARS)], 1,
          function(r) all(is.finite(as.numeric(r))))]

  K   <- length(PMED_M)
  mod_mg <- paste(
    paste0("  M", 1:K, " ~ c(a", 1:K, "_g1, a", 1:K, "_g2)*X", collapse = "\n"),
    paste0("  Y ~ c(cp_g1, cp_g2)*X + ",
           paste0("c(b", 1:K, "_g1, b", 1:K, "_g2)*M", 1:K, collapse = " + ")),
    if (K > 1) { pr <- utils::combn(K, 2)
                 paste0("  M", pr[1, ], " ~~ M", pr[2, ], collapse = "\n") } else "",
    paste0("  ind_M", 1:K, "_g1 := a", 1:K, "_g1*b", 1:K, "_g1", collapse = "\n"),
    paste0("  ind_M", 1:K, "_g2 := a", 1:K, "_g2*b", 1:K, "_g2", collapse = "\n"),
    paste0("  diff_ind_M", 1:K, " := ind_M", 1:K, "_g1 - ind_M", 1:K, "_g2",
           collapse = "\n"),
    "  diff_cprime := cp_g1 - cp_g2",
    sep = "\n\n")

  fit_mg <- tryCatch(
    lavaan::sem(mod_mg, data = d_mg, group = "grp",
                se = "bootstrap", bootstrap = PMED_BOOT, fixed.x = FALSE),
    error = function(e) { message("[17f] lavaan failed: ", conditionMessage(e)); NULL })

  if (!is.null(fit_mg)) {
    cat("\n===== GROUP DIFFERENCES IN THE SPECIFIC INDIRECT EFFECTS =====\n")
    cat("group 1 =", levels(factor(d_mg$grp))[1],
        "| group 2 =", levels(factor(d_mg$grp))[2], "\n")
    pe <- lavaan::parameterEstimates(fit_mg, ci = TRUE, boot.ci.type = PMED_CI)
    pe <- pe[pe$op == ":=", c("lhs", "est", "se", "ci.lower", "ci.upper", "pvalue")]
    names(pe) <- c("term", "est", "se", "ci_lo", "ci_hi", "p")
    for (i in seq_along(PMED_M))
      pe$term <- gsub(paste0("M", i), pmed_lab(PMED_M[i]), pe$term, fixed = TRUE)
    pe$p_fmt <- pmed_fmt_p(pe$p)
    print(pe, row.names = FALSE, digits = 3)
    pmed_modmed_tbl <- pe
    pmed_fits[["Multigroup"]] <- fit_mg
  }
}


## ---- 17g. FIGURES ------------------------------------------------------------
## Same visual language as section 15g: solid arrow = p < .05, dashed = not,
## every path annotated as coefficient (SE), p. All figures go to PMED_DIR.

## Pull one labelled parameter row out of a lavaan parameterEstimates frame.
pmed_get <- function(pe, lbl) {
  r <- pe[pe$label == lbl, ]
  if (!nrow(r)) return(list(est = NA_real_, se = NA_real_, p = NA_real_,
                            lo = NA_real_, hi = NA_real_))
  list(est = r$est[1], se = r$se[1], p = r$pvalue[1],
       lo = r$ci.lower[1], hi = r$ci.upper[1])
}

## Path diagram for one group: X on the left, the K mediators stacked in the
## middle, Y on the right. Laid out from K so it does not need editing when
## PMED_M changes length.
pmed_path_diagram <- function(g, pe, n) {

  K   <- length(PMED_M)
  sig <- function(p) if (!is.na(p) && p < 0.05) "solid" else "22"

  ## Mediators stacked in the upper half, X and Y level with the middle of that
  ## stack, and the direct X -> Y path curved underneath so it never crosses a
  ## mediator box. All positions derive from K, so nothing needs editing when
  ## PMED_M changes length.
  my   <- if (K == 1) 0.79 else seq(0.98, 0.60, length.out = K)
  ymid <- mean(range(my))

  nodes <- data.frame(
    x   = c(0.08, rep(0.50, K), 0.92),
    y   = c(ymid, my, ymid),
    lab = c(pmed_lab(PMED_X), vapply(PMED_M, pmed_lab, ""), pmed_lab(PMED_Y)),
    stringsAsFactors = FALSE)

  a  <- lapply(paste0("a", 1:K), function(l) pmed_get(pe, l))
  b  <- lapply(paste0("b", 1:K), function(l) pmed_get(pe, l))
  cp <- pmed_get(pe, "cprime")
  it <- pmed_get(pe, "ind_total")
  tt <- pmed_get(pe, "total")

  ## a-paths (X -> M) and b-paths (M -> Y)
  edges <- data.frame(
    x    = c(rep(0.145, K), rep(0.585, K)),
    y    = c(rep(ymid,  K), my),
    xend = c(rep(0.415, K), rep(0.855, K)),
    yend = c(my,            rep(ymid,  K)),
    lty  = c(vapply(a, function(z) sig(z$p), ""),
             vapply(b, function(z) sig(z$p), "")),
    stringsAsFactors = FALSE)

  ## Path labels sit at the midpoint of their own arrow, on an opaque white
  ## box (linewidth = 0) so a label crossing another arrow stays readable.
  elabs <- data.frame(
    x = c(rep(0.280, K), rep(0.720, K)),
    y = c(ymid + (my - ymid) * 0.50, my + (ymid - my) * 0.50),
    lab = c(sprintf("a%d = %.2f (%.2f), p = %s", 1:K,
                    vapply(a, function(z) z$est, 0),
                    vapply(a, function(z) z$se,  0),
                    pmed_fmt_p(vapply(a, function(z) z$p, 0))),
            sprintf("b%d = %.2f (%.2f), p = %s", 1:K,
                    vapply(b, function(z) z$est, 0),
                    vapply(b, function(z) z$se,  0),
                    pmed_fmt_p(vapply(b, function(z) z$p, 0)))),
    stringsAsFactors = FALSE)

  ## The direct X -> Y path is routed orthogonally below the whole mediator
  ## stack (down, across, up) rather than curved: a curve's belly depth depends
  ## on the panel aspect ratio, so it cuts through the bottom mediator box at
  ## some figure sizes and not others. This route is exact at any size.
  ylow <- min(my) - 0.15
  direct <- data.frame(
    x    = c(0.08,  0.08,  0.92),
    y    = c(ymid - 0.055, ylow, ylow),
    xend = c(0.08,  0.92,  0.92),
    yend = c(ylow,  ylow,  ymid - 0.055),
    head = c(FALSE, FALSE, TRUE))

  ## specific indirect effects, one line each below the diagram
  ind_lines <- data.frame(
    y = ylow - 0.20 - (seq_len(K) - 1) * 0.055,
    lab = vapply(seq_len(K), function(k) {
      z <- pmed_get(pe, paste0("ind_M", k))
      sprintf("%s:   ab = %.3f   95%% CI [%.3f, %.3f]%s",
              pmed_lab(PMED_M[k]), z$est, z$lo, z$hi,
              if (!is.na(z$lo) && !is.na(z$hi) && z$lo * z$hi > 0) "   *" else "")
    }, ""),
    stringsAsFactors = FALSE)

  ggplot() +
    ## direct effect: elbow route under the mediator stack, arrowhead on the
    ## final leg only
    geom_segment(data = direct[!direct$head, ],
                 aes(x = x, y = y, xend = xend, yend = yend),
                 linetype = sig(cp$p), linewidth = 0.75, colour = "grey25") +
    geom_segment(data = direct[direct$head, ],
                 aes(x = x, y = y, xend = xend, yend = yend),
                 linetype = sig(cp$p), linewidth = 0.75, colour = "grey25",
                 arrow = arrow(length = unit(2.8, "mm"), type = "closed")) +
    geom_segment(data = edges,
                 aes(x = x, y = y, xend = xend, yend = yend, linetype = I(lty)),
                 linewidth = 0.75, colour = "grey25",
                 arrow = arrow(length = unit(2.8, "mm"), type = "closed")) +
    geom_label(data = nodes, aes(x, y, label = lab),
               size = 4.0, fontface = "bold", label.padding = unit(2.6, "mm"),
               fill = "white", colour = "black") +
    geom_label(data = elabs, aes(x, y, label = lab),
               size = 3.0, colour = "grey15", fill = "white",
               linewidth = 0, label.padding = unit(1.1, "mm")) +
    ## as a data-frame geom, not annotate(): annotate("label") drops
    ## label.size in ggplot2 4.x and the box comes back with a border
    geom_label(data = data.frame(
                 x = 0.50, y = ylow,
                 lab = sprintf("c' = %.2f (%.2f), p = %s   [DIRECT]",
                               cp$est, cp$se, pmed_fmt_p(cp$p))),
               aes(x, y, label = lab), size = 3.1, colour = "grey15",
               fill = "white", linewidth = 0, label.padding = unit(1.4, "mm")) +
    annotate("text", x = 0.50, y = ylow - 0.10, size = 3.2, colour = "grey15",
             fontface = "bold",
             label = sprintf(
               "TOTAL INDIRECT = %.3f [%.3f, %.3f]          TOTAL c = %.3f [%.3f, %.3f]",
               it$est, it$lo, it$hi, tt$est, tt$lo, tt$hi)) +
    geom_text(data = ind_lines, aes(x = 0.50, y = y, label = lab),
              size = 3.0, colour = "grey25") +
    labs(title = sprintf("%s  --  parallel mediation, %d mediators", g, K),
         subtitle = sprintf(
           "n = %d   |   %s -> receptors -> %s   |   all variables z-scored within group",
           n, pmed_lab(PMED_X), pmed_lab(PMED_Y)),
         caption = paste0(
           "Standardized paths, annotated as coefficient (SE), p. Solid arrow = p < .05.\n",
           "b paths are each receptor's UNIQUE contribution, holding the other mediators and ",
           "ROE constant.\nBrackets are bootstrap percentile 95% CIs; * marks an interval ",
           "excluding zero.")) +
    coord_cartesian(xlim = c(0, 1),
                    ylim = c(min(ind_lines$y) - 0.06, 1.03), expand = FALSE) +
    theme_void(base_size = 13) +
    theme(plot.title      = element_text(face = "bold", hjust = 0, margin = margin(b = 2)),
          plot.subtitle   = element_text(hjust = 0, margin = margin(b = 8)),
          plot.caption    = element_text(hjust = 0, colour = "grey40", size = 8),
          plot.background = element_rect(fill = "white", colour = NA),
          plot.margin     = margin(12, 12, 10, 12))
}

for (g in names(pmed_par)) {
  p <- pmed_path_diagram(g, pmed_par[[g]], pmed_n[[g]])
  ggsave(file.path(PMED_DIR, sprintf("ParallelMediation_path_%s_%s.png",
                                     PMED_Y, gsub(" ", "", g))),
         p, width = 10, height = 7.5, dpi = 300, bg = "white")
}

## Forest plot of the specific indirect effects, one panel per group. This is
## the figure that answers "which receptor carries the ROE-pain path".
if (length(pmed_par)) {

  ind_df <- do.call(rbind, lapply(names(pmed_par), function(g) {
    pe <- pmed_par[[g]]
    do.call(rbind, lapply(seq_along(PMED_M), function(k) {
      z <- pmed_get(pe, paste0("ind_M", k))
      data.frame(group = g, mediator = pmed_lab(PMED_M[k]), n = pmed_n[[g]],
                 est = z$est, lo = z$lo, hi = z$hi, p = z$p,
                 stringsAsFactors = FALSE)
    }))
  }))

  ## keep the receptors in PMED_M order top-to-bottom, not alphabetical
  ind_df$mediator <- factor(ind_df$mediator,
                            levels = rev(vapply(PMED_M, pmed_lab, "")))

  f_ind <- ggplot(ind_df, aes(x = est, y = mediator, colour = group)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.18, linewidth = 0.8) +
    geom_point(size = 3.2) +
    facet_wrap(~ paste0(group, "  (n = ", n, ")"), ncol = 1) +
    scale_colour_manual(values = if (exists("custom_colors")) custom_colors else NULL,
                        guide = "none") +
    labs(title = "Specific indirect effects from the parallel-mediator model",
         subtitle = sprintf(
           "%s -> receptor -> %s, each receptor adjusted for the others. Bootstrap percentile 95%% CIs.",
           pmed_lab(PMED_X), pmed_lab(PMED_Y)),
         x = "Indirect effect a*b (SD units)", y = NULL) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"))

  ggsave(file.path(PMED_DIR, "ParallelMediation_indirect_effects_forest.png"),
         f_ind, width = 9, height = 7, dpi = 300, bg = "white")

  ## Effect decomposition: total = direct + total indirect, per group.
  dec_df <- do.call(rbind, lapply(names(pmed_par), function(g) {
    pe <- pmed_par[[g]]
    do.call(rbind, lapply(c("total", "cprime", "ind_total"), function(l) {
      z <- pmed_get(pe, l)
      data.frame(group = g, n = pmed_n[[g]],
                 effect = c(total = "Total (c)", cprime = "Direct (c')",
                            ind_total = "Indirect (sum ab)")[[l]],
                 est = z$est, lo = z$lo, hi = z$hi, stringsAsFactors = FALSE)
    }))
  }))
  dec_df$effect <- factor(dec_df$effect,
                          levels = c("Indirect (sum ab)", "Direct (c')", "Total (c)"))

  f_dec <- ggplot(dec_df, aes(x = est, y = effect, colour = group)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.16, linewidth = 0.8) +
    geom_point(size = 3) +
    facet_wrap(~ paste0(group, "  (n = ", n, ")"), ncol = 1) +
    scale_colour_manual(values = if (exists("custom_colors")) custom_colors else NULL,
                        guide = "none") +
    labs(title = "Effect decomposition: total = direct + total indirect",
         subtitle = sprintf("%s -> %d receptors in parallel -> %s",
                            pmed_lab(PMED_X), length(PMED_M), pmed_lab(PMED_Y)),
         x = "Effect (SD units)", y = NULL) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"))

  ggsave(file.path(PMED_DIR, "ParallelMediation_effect_decomposition.png"),
         f_dec, width = 9, height = 6.5, dpi = 300, bg = "white")

  ## Pairwise contrasts: is one receptor's indirect effect bigger than another's?
  con_df <- do.call(rbind, lapply(names(pmed_par), function(g) {
    pe <- pmed_par[[g]]
    pr <- utils::combn(length(PMED_M), 2)
    do.call(rbind, lapply(seq_len(ncol(pr)), function(j) {
      z <- pmed_get(pe, sprintf("diff_M%d_M%d", pr[1, j], pr[2, j]))
      data.frame(group = g, n = pmed_n[[g]],
                 contrast = sprintf("%s  -  %s",
                                    pmed_lab(PMED_M[pr[1, j]]),
                                    pmed_lab(PMED_M[pr[2, j]])),
                 est = z$est, lo = z$lo, hi = z$hi, stringsAsFactors = FALSE)
    }))
  }))

  f_con <- ggplot(con_df, aes(x = est, y = contrast, colour = group)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.16, linewidth = 0.8) +
    geom_point(size = 3) +
    facet_wrap(~ paste0(group, "  (n = ", n, ")"), ncol = 1) +
    scale_colour_manual(values = if (exists("custom_colors")) custom_colors else NULL,
                        guide = "none") +
    labs(title = "Pairwise contrasts between specific indirect effects",
         subtitle = "An interval excluding zero means those two receptors carry different amounts of the ROE-pain path",
         x = "Difference in indirect effect (SD units)", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"))

  ggsave(file.path(PMED_DIR, "ParallelMediation_indirect_contrasts.png"),
         f_con, width = 9, height = 7.5, dpi = 300, bg = "white")
}

## Between-group differences from the two-group model (17f).
if (!is.null(pmed_modmed_tbl)) {

  dm <- pmed_modmed_tbl[grepl("^diff_", pmed_modmed_tbl$term), ]

  f_mg <- ggplot(dm, aes(x = est, y = term)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.16,
                   linewidth = 0.8, colour = "grey25") +
    geom_point(size = 3, colour = "grey15") +
    labs(title = "Group differences in the ROE -> receptor -> pain paths",
         subtitle = sprintf("%s minus %s; two-group model, bootstrap percentile 95%% CIs",
                            PMED_GROUPS[1], PMED_GROUPS[2]),
         x = "Difference (SD units)", y = NULL,
         caption = paste0("An interval excluding zero means that path differs by drug. ",
                          "diff_cprime is the direct ROE-pain slope difference.")) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          plot.caption = element_text(hjust = 0, colour = "grey40", size = 8))

  ggsave(file.path(PMED_DIR, "ParallelMediation_group_differences.png"),
         f_mg, width = 9, height = 5.5, dpi = 300, bg = "white")
}


## ---- 17h. EXPORT -------------------------------------------------------------

pmed_sheets <- list(
  Model_syntax = data.frame(
    line = strsplit(pmed_syntax(length(PMED_M), length(PMED_COVARS)), "\n")[[1]],
    stringsAsFactors = FALSE),
  Config = data.frame(
    setting = c("X", "mediators", "Y", "groups", "covariates", "bootstrap",
                "CI type", "seed", "standardized"),
    value   = c(PMED_X, paste(PMED_M, collapse = ", "), PMED_Y,
                paste(PMED_GROUPS, collapse = ", "),
                if (length(PMED_COVARS)) paste(PMED_COVARS, collapse = ", ") else "(none)",
                PMED_BOOT, PMED_CI, PMED_SEED,
                "yes -- all variables z-scored within fitted sample"),
    stringsAsFactors = FALSE))

if (!is.null(pmed_tbl))        pmed_sheets$Per_group      <- pmed_tbl
if (!is.null(pmed_pooled_tbl)) pmed_sheets$Pooled         <- pmed_pooled_tbl
if (!is.null(pmed_modmed_tbl)) pmed_sheets$Group_diffs    <- pmed_modmed_tbl
if (!is.null(pmed_tbl))
  pmed_sheets$Indirect_only <- pmed_tbl[pmed_tbl$op == ":=", ]

writexl::write_xlsx(pmed_sheets,
  file.path(PMED_DIR, "Stats_Parallel_Mediation_ROE_receptors_NRS.xlsx"))

message("\nParallel mediation results written to: ",
        normalizePath(PMED_DIR, mustWork = FALSE))


## ---- 17i. NOTES --------------------------------------------------------------
## - Everything is standardized, so est is in SD units: a 1 SD rise in log ROE
##   moves NRS by `total` SD, of which `ind_total` SD travels through the
##   receptors and `cprime` SD does not.
## - The 4 subjects per group with ROE == 0 are dropped, not imputed. If you
##   want them in, replace log10(ROE) at section 2 with
##   log10(ROE + 0.5*min(ROE[ROE > 0])) and re-source; do it as a sensitivity
##   check and say so, because the offset is arbitrary.
## - Percentile CIs, not BCa: at n < 20 the BCa acceleration estimate is
##   unstable. Same choice as section 15.
## - Cross-check: with K = 1 this reduces to section 15's single-mediator model
##   and ind_M1 should match ab_std for that receptor to within bootstrap noise.
## - The model is saturated in the paths but not in the residual covariances;
##   global fit indices are near-perfect by construction and should not be
##   reported as evidence the model is right.

