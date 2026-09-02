################################################################################
## 20. BETWEEN-DRUG COMPARISON: DO THE TWO DRUGS SHARE THE SAME SLOPES?
##     Pooled models, drug carried as a predictor.
##
## THE SPLIT WITH SECTION 18
##   Section 18 (07_nested_regression.R) fits every model SEPARATELY inside each
##   drug group -- one model in Hydrocodone (n = 19), one in Tramadol (n = 17).
##   That is the stratified estimate, and it is the right thing to report IF the
##   two drugs really do behave differently. But two separate models can never
##   TEST that: a predictor coming out significant in one group and not the
##   other is not evidence that the two groups differ (the difference between
##   "significant" and "non-significant" is not itself significant).
##
##   This section is the test. Both groups go into ONE model with drug as a
##   predictor, so the difference between them becomes an estimable quantity
##   with a standard error and a p-value.
##
##       07 models each drug.   09 compares the drugs.
##
## THE TWO MODELS
##   M_add   y = b0 + b_ROE*log_ROE + b_HT4*R_5HT4 + b_HT6*R_5HT6
##                + b_G*drug + e
##   M_int   M_add + drug x log_ROE + drug x R_5HT4 + drug x R_5HT6
##
##   M_add assumes both drugs share the same slopes and differ only in average
##   pain level -- two parallel lines at different heights. b_G is that height
##   difference and nothing more.
##
##   M_int drops that assumption. EACH INTERACTION COEFFICIENT IS THE SLOPE
##   DIFFERENCE BETWEEN THE DRUGS for that predictor. The nested F test of M_add
##   against M_int, on 3 df, is the omnibus question: do the drugs differ on any
##   slope at all? That test is the justification for section 18's entire
##   stratified design, which until now was asserted in a header comment rather
##   than tested.
##
##   Both models are fit for EVERY model set in MREG_MODELS, not just the
##   two-receptor one, so the comparison is available wherever section 18 has a
##   stratified pair. The HT4_HT6 set is the one reported as the headline pair,
##   because it is the set the marginal-effect and interaction figures use.
##
## SIGN CONVENTION -- READ THIS BEFORE READING ANY COEFFICIENT
##   Drug enters as a FACTOR with Hydrocodone as the reference level (the order
##   of FOCUS_GROUPS in 00_config.R). So:
##
##       every group coefficient reads  TRAMADOL MINUS HYDROCODONE.
##
##   Two consequences worth stating, because both have caused confusion:
##
##     (a) In M_int the plain R_5HT4 row is not "the 5-HT4 effect overall". It
##         is the 5-HT4 slope IN HYDROCODONE, because that is the group the
##         factor puts at zero. The tramadol slope is that row PLUS the
##         interaction row. 20d does this arithmetic and reports both drugs'
##         slopes directly, so nobody has to do it by hand.
##
##     (b) The analysis plan we were sent specified drug coded 1 = Tramadol,
##         2 = Hydrocodone, which runs the comparison the other way round. That
##         was deliberately NOT followed. Which drug is the reference changes
##         nothing except the sign: identical fitted values, identical R2, F,
##         standard errors and p-values. Matching the direction section 18
##         already used everywhere is worth more than matching an arbitrary
##         coding choice, since otherwise the same quantity would appear with
##         opposite signs in two neighbouring outputs.
##
## WHY NOTHING IS REFIT HERE
##   The pooled fits in 20b are the ONLY place the pooled NRS models are fit.
##   Section 18 used to fit them and report only the interaction rows; that
##   loop moved here whole, into preg_fit_all() -- the ONLY fitting logic for
##   a pooled model, for any outcome: 20b calls it once for NRS, 20j calls it
##   again for SOWS, and neither outcome has a second, independently-written
##   copy of the loop. 20c reads the stored `lm` objects and describes them.
##   There is deliberately no second fitting path for a given outcome, so the
##   coefficient table and the interaction table cannot drift apart.
##
## EXPOSURE: log_ROE, n = 36 or 44 depending on ROE_ZERO_IMPUTE (00_config.R)
##   log_ROE = log10(ROE) is -Inf for the subjects with ROE == 0. With
##   ROE_ZERO_IMPUTE FALSE the is.finite() filter drops them -- same rule as
##   sections 15-18 -- so the pooled n (36) is exactly the two stratified n's
##   added together (19 + 17). With it TRUE (the default) those 8 subjects
##   get a finite, imputed log_ROE instead (01_load_data.R) and the pooled n
##   rises to 44 (23 + 21).
##
##   The source spreadsheet also carries its own logROE column, computed as
##   LOG10(100000*ROE). The x100000 is algebraically a +5 shift and is
##   statistically inert (a constant added to a predictor moves only the
##   intercept). Its ROE == 0 rows, however, hold a hand-typed 0 where the
##   formula returned #NUM!, and on that scale 0 means ROE = 1e-5 mg/L -- a real
##   concentration BELOW the smallest measured value. Those 8 points would enter
##   at an invented x carrying ~25% of the leverage, and they move the omnibus
##   test from p = .0003 to p = .096. That column is therefore not used; see the
##   note at 01_load_data.R for both that and the (larger) leverage cost of the
##   ROE_ZERO_IMPUTE default.
##
## Self-contained: needs `master`, `FOCUS_GROUPS`, `OUT_DIR` (00_config.R /
## 01_load_data.R) and section 18's helpers and fitted objects
## (07_nested_regression.R). Sources whatever is missing. No new packages.
##
## Run standalone: Rscript R/09_pooled_regression.R
## Sourced by: Run_All.R, AFTER 07 and BEFORE 02 (02 builds Table 4 from the
## objects created here).
################################################################################

.HTK_PROJECT_DIR <- "/Users/kushagraverma/Work/Projects/Hydrocodone+Tramadol_Project"
## Always re-sourced -- see the note above this line in 01_load_data.R.
source(file.path(.HTK_PROJECT_DIR, "R", "00_config.R"))
if (!exists("master"))       source(file.path(.HTK_PROJECT_DIR, "R", "01_load_data.R"))
## Everything from section 18: the helpers (mreg_theme/save/lab/short/fmt_p/
## stars/model_short/table_figure), the constants (MREG_Y, MREG_MODELS, ...) and
## the stratified fits (mreg_fits, mreg_coef_tbl) that figure (7) draws.
## In a Run_All.R run 07 has already executed and this does not fire.
if (!exists("mreg_fit_one")) source(file.path(.HTK_PROJECT_DIR, "R", "07_nested_regression.R"))


## ---- 20a. CONFIG -------------------------------------------------------------

## Kept separate from MREG_DIR: these outputs are about the comparison between
## drugs, not about either drug's own model.
PREG_DIR <- file.path(OUT_DIR, "between_drug_comparison")
dir.create(PREG_DIR, showWarnings = FALSE, recursive = TRUE)

## The model set whose pooled pair is reported as the headline. Same constant
## the marginal-effect and interaction figures key off, and the same set
## sections 18d/18e treat as the primary two-receptor model.
MREG_FULL <- "HT4_HT6"

## The two pooled models, named by the role they play in the comparison rather
## than by their predictors (their predictors are MREG_MODELS[[MREG_FULL]] plus
## drug, plus the interactions). Term vectors are in the order the analysis plan
## listed them -- receptors, then exposure, then drug, then the interactions --
## which is also the order they appear in every table below.
PREG_MODELS <- list(
  Pooled_main   = c("R_5HT4", "R_5HT6", "log_ROE", "group_TvsH"),
  Pooled_grpint = c("R_5HT4", "R_5HT6", "log_ROE", "group_TvsH",
                    "group_x_5HT4", "group_x_5HT6", "group_x_logROE")
)

PREG_MODEL_LABELS <- c(
  Pooled_main   = "5-HT4 + 5-HT6 + ROE + drug",
  Pooled_grpint = "5-HT4 + 5-HT6 + ROE + drug + drug x each predictor"
)

## Short labels for the table figure's column headers, same hand-wrapped shape
## as MREG_MODEL_SHORT (leading tag, then what the model adds).
PREG_MODEL_SHORT <- c(
  Pooled_main   = "P1\n5-HT4 + 5-HT6\n+ ROE + drug",
  Pooled_grpint = "P2\n+ drug x\neach predictor"
)

## One "group" label, because a pooled model is not fit within a group. Used as
## the `group` key throughout so mreg_table_figure() and build_regression_table()
## -- both of which filter on that column -- work unchanged.
PREG_GROUP_LAB <- "Both groups (pooled)"


## ---- 20b. POOLED FITS --------------------------------------------------------
## Fitting the two drugs separately shows the slopes look different; it does not
## test that they ARE different. Pool the two groups and interact every
## predictor with drug: each interaction coefficient IS the raw-unit slope
## difference (tramadol minus hydrocodone, given the factor's level order), and
## the model-level F test asks whether the drugs differ on any slope at all.
##
## Moved here from section 18e unchanged. This is the only place in the pipeline
## where the pooled NRS models are fit -- and, per outcome, the only place any
## pooled model is fit: preg_fit_all() below is the fitting logic, called once
## here for NRS and again in 20j for SOWS, so neither outcome gets a second,
## drifting copy of the loop.
##
## preg_fit_all() returns four pieces: `pool` (the group-filtered data used for
## every model key), `inter_fits` (named by model key: list(add, int, data,
## preds)), `interaction_tbl` (one row per drug x predictor interaction term)
## and `intomni_tbl` (the nested-F omnibus row per model key).
preg_fit_all <- function(y_var, model_defs) {
  pool <- master %>%
    dplyr::filter(!is.na(Plot_Group), Plot_Group %in% MREG_GROUPS) %>%
    dplyr::mutate(.grp = factor(as.character(Plot_Group), levels = MREG_GROUPS))

  inter_rows   <- list()
  intomni_rows <- list()
  inter_fits   <- list()

  for (mk in names(model_defs)) {
    preds <- model_defs[[mk]]
    vars  <- c(y_var, preds, ".grp")
    if (!all(setdiff(vars, ".grp") %in% names(pool))) next

    dp <- as.data.frame(pool[, vars, drop = FALSE])
    dp <- dp[stats::complete.cases(dp), , drop = FALSE]
    num <- setdiff(vars, ".grp")
    dp  <- dp[apply(dp[, num, drop = FALSE], 1, function(r) all(is.finite(r))), , drop = FALSE]
    if (nrow(dp) < MREG_MIN_N || nlevels(droplevels(dp$.grp)) < 2) next

    rhs   <- paste(sprintf("`%s`", preds), collapse = " + ")
    m_add <- stats::lm(stats::as.formula(sprintf("`%s` ~ .grp + %s", y_var, rhs)), data = dp)
    m_int <- stats::lm(stats::as.formula(sprintf("`%s` ~ .grp * (%s)", y_var, rhs)), data = dp)
    inter_fits[[mk]] <- list(add = m_add, int = m_int, data = dp, preds = preds)

    ci <- suppressWarnings(stats::confint(m_int))
    cf <- summary(m_int)$coefficients
    keep <- grep("^\\.grp.*:", rownames(cf))

    if (length(keep)) {
      trm <- rownames(cf)[keep]
      inter_rows[[mk]] <- tibble::tibble(
        model      = mk,
        term       = gsub("`", "", trm),
        predictor  = gsub("`", "", sub("^.*:", "", trm)),
        comparison = paste0(MREG_GROUPS[2], " minus ", MREG_GROUPS[1]),
        n          = nrow(dp),
        est        = cf[keep, 1],
        se         = cf[keep, 2],
        t          = cf[keep, 3],
        p          = cf[keep, 4],
        ci_lo      = ci[match(trm, rownames(ci)), 1],
        ci_hi      = ci[match(trm, rownames(ci)), 2]
      )
    }

    av <- stats::anova(m_add, m_int)
    intomni_rows[[mk]] <- tibble::tibble(
      model  = mk, n = nrow(dp),
      df_num = av$Df[2], df_den = av$Res.Df[2],
      F_stat = av$F[2],  p = av$`Pr(>F)`[2]
    )
  }

  list(pool = pool, inter_fits = inter_fits,
       interaction_tbl = dplyr::bind_rows(inter_rows),
       intomni_tbl     = dplyr::bind_rows(intomni_rows))
}

preg_nrs_fit <- preg_fit_all(MREG_Y, MREG_MODELS)

mreg_pool             <- preg_nrs_fit$pool
mreg_inter_fits       <- preg_nrs_fit$inter_fits
mreg_interaction_tbl  <- preg_nrs_fit$interaction_tbl
mreg_intomni_tbl      <- preg_nrs_fit$intomni_tbl

cat("\n================================================================\n")
cat("20. BETWEEN-DRUG COMPARISON: ", MREG_Y, " ~ drug x (log ROE + 5-HT4 / 5-HT6)\n", sep = "")
cat("    All coefficients read ", MREG_GROUPS[2], " MINUS ", MREG_GROUPS[1], ".\n", sep = "")
cat("================================================================\n")

if (nrow(mreg_interaction_tbl)) {
  mreg_interaction_tbl <- mreg_interaction_tbl %>%
    dplyr::mutate(model_lab = unname(MREG_MODEL_LABELS[model]),
                  term_lab  = mreg_short(predictor),
                  sig       = !is.na(p) & p < 0.05,
                  p_fmt     = mreg_fmt_p(p))
  cat("\n----- SLOPE DIFFERENCES BETWEEN DRUGS (interaction terms, raw units) -----\n")
  cat("(", MREG_GROUPS[2], " minus ", MREG_GROUPS[1], ")\n\n", sep = "")
  print(as.data.frame(mreg_interaction_tbl %>%
          dplyr::select(model, predictor, n, est, se, ci_lo, ci_hi, t, p)),
        row.names = FALSE, digits = 3)
}
if (nrow(mreg_intomni_tbl)) {
  mreg_intomni_tbl$p_fmt <- mreg_fmt_p(mreg_intomni_tbl$p)
  cat("\n----- OMNIBUS TEST: do any slopes differ between drugs? -----\n\n")
  print(as.data.frame(mreg_intomni_tbl), row.names = FALSE, digits = 3)
}


## ---- 20c. DESCRIBE THE POOLED MODELS -----------------------------------------
## Section 18e kept only the interaction rows out of these fits and threw the
## rest away -- no coefficient table, no R2, no n, nothing a reader could put in
## a paper. This block reads the SAME fitted objects and produces the tables in
## the shape mreg_fit_one() returns, so mreg_table_figure() (18h) and
## build_regression_table() (02) both accept them with no changes.
##
## Nothing here fits a model.

## lm() names a backticked predictor with its backticks, and the factor's main
## effect as ".grp<level>". Both are unusable as table row keys: the backticks
## break name matching and ".grpTramadol group" is not something MED_LABELS can
## be expected to hold. Rename to stable keys once, here, and every downstream
## consumer sees ordinary names.
##
##   .grpTramadol group             -> group_TvsH
##   .grpTramadol group:`R_5HT4`    -> group_x_5HT4
##   .grpTramadol group:`log_ROE`   -> group_x_logROE
##
## Keyed off MREG_GROUPS[2] rather than the literal string, so renaming a group
## in 00_config.R cannot leave a stale prefix here.
preg_tag <- function(v) {
  v <- sub("^R_", "", v)          # R_5HT4  -> 5HT4
  sub("^log_", "log", v)          # log_ROE -> logROE
}
preg_clean_terms <- function(raw) {
  tm   <- gsub("`", "", raw)
  pref <- paste0(".grp", MREG_GROUPS[2])
  out  <- tm
  is_int  <- startsWith(tm, paste0(pref, ":"))
  out[is_int] <- paste0("group_x_",
                        preg_tag(substring(tm[is_int], nchar(pref) + 2)))
  out[tm == pref] <- "group_TvsH"
  out
}

## VIF straight off the model matrix rather than by re-regressing named
## predictors as mreg_fit_one() does: the interacted model's columns include a
## factor dummy and three products, which are not columns of `data`.
##
## Expect LARGE values on the interaction terms. They are literally products of
## the model's own main effects, so they are collinear with them by
## construction. That inflates the MAIN-EFFECT standard errors; it does not
## touch the interaction tests or the nested F, which are what this section is
## for. A high max VIF here is not a broken model.
preg_vif <- function(m) {
  X <- stats::model.matrix(m)
  X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  if (ncol(X) < 2) return(stats::setNames(rep(NA_real_, ncol(X)), colnames(X)))
  v <- vapply(seq_len(ncol(X)), function(j) {
    r2 <- summary(stats::lm(X[, j] ~ X[, -j, drop = FALSE]))$r.squared
    if (is.finite(r2) && r2 < 1) 1 / (1 - r2) else Inf
  }, numeric(1))
  stats::setNames(v, colnames(X))
}

## One fitted lm -> the (coef, fit) pair the table builders expect.
##
## NOTE: no standardized betas. Producing them would mean refitting on z-scored
## data, and this section's whole point is that these models are fit exactly
## once. Worse, a standardized product term is not the product of the
## standardized terms, so the number would not mean what beta_std means in
## section 18. Raw units only -- stated in the Config sheet too.
preg_describe <- function(m, model_key, group_lab = PREG_GROUP_LAB) {
  s   <- summary(m)
  cf  <- s$coefficients
  ci  <- suppressWarnings(stats::confint(m))
  raw <- rownames(cf)

  vif_raw <- preg_vif(m)
  vif_map <- stats::setNames(unname(vif_raw), preg_clean_terms(names(vif_raw)))

  coef_tbl <- tibble::tibble(
    group = group_lab,
    model = model_key,
    term  = preg_clean_terms(raw),
    est   = cf[, 1],
    se    = cf[, 2],
    t     = cf[, 3],
    p     = cf[, 4],
    ci_lo = ci[match(raw, rownames(ci)), 1],
    ci_hi = ci[match(raw, rownames(ci)), 2]
  )
  coef_tbl$beta_std <- NA_real_
  coef_tbl$vif      <- unname(vif_map[coef_tbl$term])

  fst <- s$fstatistic
  res <- stats::residuals(m)
  n   <- length(res)

  ## Same two diagnostics section 18 reports, computed the same way: Shapiro-Wilk
  ## on the residuals, and a hand-rolled Breusch-Pagan (regress squared residuals
  ## on the design, n * R2 is chi-square on #predictors df).
  sw_p <- tryCatch(stats::shapiro.test(res)$p.value, error = function(e) NA_real_)
  bp_p <- tryCatch({
    X <- stats::model.matrix(m)
    r2a <- summary(stats::lm(res^2 ~ X[, colnames(X) != "(Intercept)", drop = FALSE]))$r.squared
    stats::pchisq(n * r2a, df = length(stats::coef(m)) - 1, lower.tail = FALSE)
  }, error = function(e) NA_real_)

  fit_tbl <- tibble::tibble(
    group      = group_lab,
    model      = model_key,
    formula    = paste(deparse(stats::formula(m)), collapse = " "),
    n          = n,
    k          = length(stats::coef(m)) - 1,
    R2         = s$r.squared,
    adj_R2     = s$adj.r.squared,
    F_stat     = if (!is.null(fst)) unname(fst[1]) else NA_real_,
    df1        = if (!is.null(fst)) unname(fst[2]) else NA_real_,
    df2        = if (!is.null(fst)) unname(fst[3]) else NA_real_,
    F_p        = if (!is.null(fst)) stats::pf(fst[1], fst[2], fst[3], lower.tail = FALSE) else NA_real_,
    sigma      = s$sigma,
    RMSE       = sqrt(mean(res^2)),
    AIC        = stats::AIC(m),
    BIC        = stats::BIC(m),
    max_VIF    = if (all(is.na(vif_raw))) NA_real_ else max(vif_raw, na.rm = TRUE),
    max_cooksD = max(stats::cooks.distance(m), na.rm = TRUE),
    shapiro_p  = sw_p,
    bp_p       = bp_p
  )

  list(coef = coef_tbl, fit = fit_tbl)
}

preg_fi <- mreg_inter_fits[[MREG_FULL]]
if (is.null(preg_fi))
  stop("09_pooled_regression.R: no pooled fit for '", MREG_FULL,
       "' -- 20b fit nothing, so there is nothing to report.")

preg_desc <- list(Pooled_main   = preg_describe(preg_fi$add, "Pooled_main"),
                  Pooled_grpint = preg_describe(preg_fi$int, "Pooled_grpint"))

preg_coef_tbl <- dplyr::bind_rows(lapply(preg_desc, `[[`, "coef")) %>%
  dplyr::mutate(model_lab = unname(PREG_MODEL_LABELS[model]),
                term_lab  = ifelse(term == "(Intercept)", "Intercept",
                                   vapply(term, mreg_lab, character(1))),
                p_fmt     = mreg_fmt_p(p))
preg_fit_tbl <- dplyr::bind_rows(lapply(preg_desc, `[[`, "fit"))

## The nested F for THIS pair, pulled from the table 20b already built rather
## than re-running anova() -- one source for the number.
preg_nested_tbl <- mreg_intomni_tbl %>%
  dplyr::filter(model == MREG_FULL) %>%
  dplyr::mutate(reduced = "Pooled_main", full = "Pooled_grpint",
                added_terms = "drug x each predictor") %>%
  dplyr::select(reduced, full, added_terms, n, df_num, df_den, F_stat, p, p_fmt)

cat("\n----- POOLED MODEL COEFFICIENTS (raw units) -----\n")
cat("Model 2's main effects are the ", MREG_GROUPS[1],
    " slopes (the reference group); see 20d.\n\n", sep = "")
print(as.data.frame(preg_coef_tbl %>%
        dplyr::select(model, term, est, se, ci_lo, ci_hi, t, p, vif)),
      row.names = FALSE, digits = 3)

cat("\n----- POOLED MODEL FIT -----\n\n")
print(as.data.frame(preg_fit_tbl %>%
        dplyr::select(model, n, R2, adj_R2, F_stat, df1, df2, F_p,
                      RMSE, AIC, max_VIF, shapiro_p, bp_p)),
      row.names = FALSE, digits = 3)

cat("\n----- DO THE DRUGS DIFFER? (nested F, Model 1 vs Model 2) -----\n\n")
print(as.data.frame(preg_nested_tbl), row.names = FALSE, digits = 3)


## ---- 20d. EACH DRUG'S OWN SLOPE, READ OUT OF THE INTERACTION MODEL -----------
## In Model 2 the plain R_5HT4 row is the 5-HT4 slope in the REFERENCE group
## (Hydrocodone), and the tramadol slope is that row plus its interaction row.
## Nobody should have to do that addition by hand off a coefficient table, so it
## is done here for every predictor and reported as its own table.
##
## The standard error of the sum is NOT se(main) + se(inter) and not the square
## root of the sum of squares either -- the two estimates are correlated, so
##     Var(a + b) = Var(a) + Var(b) + 2*Cov(a, b)
## and the covariance comes off vcov(). Getting this wrong would understate the
## tramadol CIs.
##
## These numbers are also the section's strongest self-check: a model with every
## predictor interacted with group is algebraically identical to fitting the two
## groups separately, so what comes out here MUST equal section 18's stratified
## HT4_HT6 coefficients. 20g prints that comparison.
##
## Pulled out as a function of the fitted (add, int, data, preds) list rather
## than closing over `preg_fi`, so 20j can derive the same two numbers for the
## pooled SOWS interaction model without a second copy of this arithmetic.
preg_group_slopes_for <- function(fit_obj) {
  m   <- fit_obj$int
  b   <- stats::coef(m)
  V   <- stats::vcov(m)
  raw <- names(b)
  cl  <- preg_clean_terms(raw)
  df  <- stats::df.residual(m)
  tcrit <- stats::qt(0.975, df)

  ## match(), not idx[[nm]] -- indexing a named vector by a name it does not
  ## have throws rather than returning NULL, which would turn a renaming bug
  ## into a crash instead of a skipped row.
  gi <- function(nm) { j <- match(nm, cl); if (is.na(j)) NULL else j }

  rows <- list()
  for (pv in fit_obj$preds) {
    i_main <- gi(pv)
    i_int  <- gi(paste0("group_x_", preg_tag(pv)))
    if (is.null(i_main) || is.null(i_int)) next

    for (g in MREG_GROUPS) {
      if (identical(g, MREG_GROUPS[1])) {          # reference group: main effect alone
        est <- b[i_main]
        vr  <- V[i_main, i_main]
      } else {                                      # other group: main + interaction
        est <- b[i_main] + b[i_int]
        vr  <- V[i_main, i_main] + V[i_int, i_int] + 2 * V[i_main, i_int]
      }
      se <- sqrt(vr)
      rows[[paste(pv, g)]] <- tibble::tibble(
        predictor = pv, predictor_lab = mreg_lab(pv), group = g,
        est = unname(est), se = unname(se),
        ci_lo = unname(est - tcrit * se), ci_hi = unname(est + tcrit * se),
        t = unname(est / se),
        p = 2 * stats::pt(abs(unname(est / se)), df, lower.tail = FALSE))
    }
  }
  out <- dplyr::bind_rows(rows)
  if (nrow(out)) out$p_fmt <- mreg_fmt_p(out$p)
  out
}

preg_group_slopes <- preg_group_slopes_for(preg_fi)

if (nrow(preg_group_slopes)) {
  cat("\n----- EACH DRUG'S SLOPE, DERIVED FROM THE INTERACTION MODEL -----\n")
  cat("(not separately fit: main effect for ", MREG_GROUPS[1],
      ", main + interaction for ", MREG_GROUPS[2], ")\n\n", sep = "")
  print(as.data.frame(preg_group_slopes %>%
          dplyr::select(predictor, group, est, se, ci_lo, ci_hi, t, p)),
        row.names = FALSE, digits = 3)
}


## ---- 20e. FIGURES ------------------------------------------------------------
## Moved from section 18i unchanged, except that they write to PREG_DIR. All
## three answer "how do the two drugs compare", which is what makes them belong
## here rather than beside the within-drug figures.

## (7) Marginal (adjusted) effect of each receptor from the model HT4_HT6:
## predicted NRS across the observed range of that receptor with the other
## predictors held at their group means, 95% CI ribbon, raw data overlaid as
## partial residuals so the points sit on the same scale as the line.
## HT4_HT6, the two-receptor main-effects model, is used below for the
## per-receptor marginal-effect plots (7) and the drug-interaction-slopes plot
## (8). Deliberately hardcoded rather than "last model in the list", for two
## separate reasons:
##
##   - M4 and M5 both contain product terms. Both plots hold every OTHER
##     predictor at its mean while sweeping one receptor across the x-axis,
##     which is wrong when "other predictors" includes a product built from
##     the very receptor being swept -- the product would stay frozen at its
##     mean instead of moving with the receptor, so the drawn line would not
##     match what the model actually predicts.
##   - M5 would additionally freeze three covariates at their means to draw one
##     receptor's line. That line is not wrong, but it is a picture of an
##     8-predictor model at n ~ 13-18 and reads as far more precise than it is.
##
## The other models still get their own coefficient forest, regression table,
## nested test and model-comparison entries in section 18; they just aren't the
## model behind THESE two figures.

for (rv in setdiff(MREG_MODELS[[MREG_FULL]], MREG_X)) {
  line_rows <- list(); pt_rows <- list(); lab_rows <- list()

  for (g in names(mreg_fits)) {
    res <- mreg_fits[[g]][[MREG_FULL]]
    if (is.null(res) || !rv %in% res$preds) next
    d <- res$data; m <- res$model

    grid <- data.frame(seq(min(d[[rv]]), max(d[[rv]]), length.out = 120))
    names(grid) <- rv
    for (o in setdiff(res$preds, rv)) grid[[o]] <- mean(d[[o]])
    pr <- stats::predict(m, newdata = grid, interval = "confidence")
    line_rows[[g]] <- tibble::tibble(group = g, x = grid[[rv]],
                                     fit = pr[, 1], lo = pr[, 2], hi = pr[, 3])

    ## Partial residuals: observed minus the fitted contribution of the OTHER
    ## predictors, so each point is "this subject's NRS with exposure and the
    ## other receptor netted out" -- directly comparable to the line.
    b <- stats::coef(m)
    contrib <- rep(0, nrow(d))
    for (o in setdiff(res$preds, rv)) contrib <- contrib + b[[o]] * (d[[o]] - mean(d[[o]]))
    pt_rows[[g]] <- tibble::tibble(group = g, x = d[[rv]], y = d[[res$yv]] - contrib)

    cr <- res$coef[res$coef$term == rv, ]
    lab_rows[[g]] <- tibble::tibble(
      group = g,
      lab = sprintf("%s\nb = %.3f (SE %.3f)\n95%% CI [%.3f, %.3f]\np = %s,  n = %d",
                    g, cr$est[1], cr$se[1], cr$ci_lo[1], cr$ci_hi[1],
                    mreg_fmt_p(cr$p[1]), nrow(d)))
  }

  if (!length(line_rows)) next
  ln <- dplyr::bind_rows(line_rows); pt <- dplyr::bind_rows(pt_rows)
  lb <- dplyr::bind_rows(lab_rows)

  ip <- nrow(mreg_interaction_tbl) &&
    any(mreg_interaction_tbl$model == MREG_FULL & mreg_interaction_tbl$predictor == rv)
  sub_int <- if (ip) {
    r <- mreg_interaction_tbl[mreg_interaction_tbl$model == MREG_FULL &
                                mreg_interaction_tbl$predictor == rv, ]
    sprintf("Slope difference between drugs: %.3f (95%% CI [%.3f, %.3f]), p = %s",
            r$est[1], r$ci_lo[1], r$ci_hi[1], mreg_fmt_p(r$p[1]))
  } else NULL

  p_mar <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(data = ln, ggplot2::aes(x, ymin = lo, ymax = hi, fill = group),
                         alpha = 0.16) +
    ggplot2::geom_line(data = ln, ggplot2::aes(x, fit, colour = group), linewidth = 1.1) +
    ggplot2::geom_point(data = pt, ggplot2::aes(x, y, colour = group), size = 2.7, alpha = 0.85) +
    ggplot2::scale_colour_manual(values = MREG_COLS, name = NULL) +
    ggplot2::scale_fill_manual(values = MREG_COLS, guide = "none") +
    ggplot2::labs(
      title    = sprintf("Adjusted effect of %s on %s", mreg_lab(rv), mreg_lab(MREG_Y)),
      subtitle = paste0(unname(MREG_MODEL_LABELS[MREG_FULL]),
                        ", other predictors held at their group means.",
                        if (!is.null(sub_int)) paste0("\n", sub_int) else ""),
      x = mreg_lab(rv), y = sprintf("%s, adjusted (partial residual)", mreg_lab(MREG_Y)),
      caption = paste0("Points are partial residuals: each subject's ", mreg_lab(MREG_Y),
                       " with the other predictors' contributions removed, so they sit on the ",
                       "same scale as the fitted line.\nRibbon is the 95% CI of the fit.")
    ) +
    mreg_theme()

  ## Stats box in the top-right, one block per group.
  ## NOTE: drawn with geom_label on an inline data.frame, not annotate("label").
  ## ggplot2 4.x drops label.size from annotate() with only a warning, and the
  ## box border comes back -- same trap as section 17's path diagrams.
  yr <- range(c(pt$y, ln$hi, ln$lo))
  box <- data.frame(x = max(ln$x), y = yr[2] + 0.14 * diff(yr),
                    lab = paste(lb$lab, collapse = "\n\n"),
                    stringsAsFactors = FALSE)
  p_mar <- p_mar +
    ggplot2::geom_label(data = box, ggplot2::aes(x, y, label = lab),
                        inherit.aes = FALSE, hjust = 1, vjust = 1, size = 2.9,
                        fill = "white", alpha = 0.85, linewidth = 0,
                        lineheight = 0.98) +
    ggplot2::expand_limits(y = yr[2] + 0.16 * diff(yr))

  mreg_save(p_mar, sprintf("MREG_partial_%s.png", sub("^R_", "", rv)),
            w = 9, h = 7, dir = PREG_DIR)
}

## (8) The interaction picture: each drug's fitted line for each predictor, from
## the full interaction model, with the slope-difference p annotated. This is
## the visual form of the test in 20b.
if (length(mreg_inter_fits) && !is.null(mreg_inter_fits[[MREG_FULL]])) {
  fi <- mreg_inter_fits[[MREG_FULL]]
  m  <- fi$int; dp <- fi$data

  ln_rows <- list(); pt_rows <- list()
  for (pv in fi$preds) {
    for (g in MREG_GROUPS) {
      dg <- dp[dp$.grp == g, , drop = FALSE]
      if (!nrow(dg)) next
      grid <- data.frame(seq(min(dg[[pv]]), max(dg[[pv]]), length.out = 100))
      names(grid) <- pv
      grid$.grp <- factor(g, levels = MREG_GROUPS)
      for (o in setdiff(fi$preds, pv)) grid[[o]] <- mean(dg[[o]])
      pr <- stats::predict(m, newdata = grid, interval = "confidence")
      ln_rows[[paste(pv, g)]] <- tibble::tibble(
        predictor = pv, group = g, x = grid[[pv]],
        fit = pr[, 1], lo = pr[, 2], hi = pr[, 3])

      ## Partial residuals under the interaction model: this group's slope for
      ## the other predictors is the main effect PLUS its interaction term.
      ## The reference level (MREG_GROUPS[1]) has no interaction term at all,
      ## so look it up defensively rather than indexing blind.
      b <- stats::coef(m)
      bget <- function(nm) {
        v <- b[gsub("`", "", names(b)) == nm]
        if (!length(v) || is.na(v[1])) 0 else unname(v[1])
      }
      contrib <- rep(0, nrow(dg))
      for (o in setdiff(fi$preds, pv)) {
        slope_g <- bget(o) + bget(paste0(".grp", g, ":", o))
        contrib <- contrib + slope_g * (dg[[o]] - mean(dg[[o]]))
      }
      pt_rows[[paste(pv, g)]] <- tibble::tibble(
        predictor = pv, group = g, x = dg[[pv]], y = dg[[MREG_Y]] - contrib)
    }
  }

  if (length(ln_rows)) {
    lnI <- dplyr::bind_rows(ln_rows); ptI <- dplyr::bind_rows(pt_rows)

    facet_lab <- vapply(fi$preds, function(pv) {
      r <- mreg_interaction_tbl[mreg_interaction_tbl$model == MREG_FULL &
                                  mreg_interaction_tbl$predictor == pv, ]
      if (!nrow(r)) return(mreg_lab(pv))
      sprintf("%s\nslope difference = %.3f, p = %s",
              mreg_lab(pv), r$est[1], mreg_fmt_p(r$p[1]))
    }, character(1))

    lnI$panel <- factor(facet_lab[lnI$predictor], levels = unname(facet_lab))
    ptI$panel <- factor(facet_lab[ptI$predictor], levels = unname(facet_lab))

    p_int <- ggplot2::ggplot() +
      ggplot2::geom_ribbon(data = lnI, ggplot2::aes(x, ymin = lo, ymax = hi, fill = group),
                           alpha = 0.15) +
      ggplot2::geom_line(data = lnI, ggplot2::aes(x, fit, colour = group), linewidth = 1.1) +
      ggplot2::geom_point(data = ptI, ggplot2::aes(x, y, colour = group),
                          size = 2.5, alpha = 0.85) +
      ggplot2::scale_colour_manual(values = MREG_COLS, name = NULL) +
      ggplot2::scale_fill_manual(values = MREG_COLS, guide = "none") +
      ggplot2::facet_wrap(~ panel, scales = "free_x", nrow = 1) +
      ggplot2::labs(
        title    = "Do the slopes differ between drugs?",
        subtitle = sprintf("Pooled model %s ~ drug x (%s), n = %d.  Slope difference = %s minus %s",
                           mreg_lab(MREG_Y), paste(mreg_short(fi$preds), collapse = " + "),
                           nrow(dp), MREG_GROUPS[2], MREG_GROUPS[1]),
        x = "Predictor (raw units)",
        y = sprintf("%s, adjusted (partial residual)", mreg_lab(MREG_Y)),
        caption = paste0("Each panel's p is the drug x predictor interaction: the test that the two ",
                         "lines are not parallel.\nNote log10 ROE spans a very different range in the ",
                         "two drugs, so its two slopes are estimated over barely overlapping x.")
      ) +
      mreg_theme()

    mreg_save(p_int, "MREG_interaction_slopes.png",
              w = 4.6 * length(fi$preds) + 1.2, h = 6.2, dir = PREG_DIR)
  }
}

## (9) Slope differences as a forest: every interaction coefficient from every
## model, with its CI. A CI clear of zero means the drugs genuinely differ on
## that slope.
if (nrow(mreg_interaction_tbl)) {
  idf <- mreg_interaction_tbl %>%
    dplyr::mutate(
      model_lab = factor(model_lab, levels = unname(MREG_MODEL_LABELS)),
      term_lab  = factor(term_lab,
                         levels = rev(unique(mreg_short(unlist(MREG_MODELS)))))
    )

  p_gd <- ggplot2::ggplot(idf, ggplot2::aes(est, term_lab)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = ci_lo, xmax = ci_hi),
                            orientation = "y", width = 0.16, linewidth = 0.8, colour = "grey30") +
    ggplot2::geom_point(ggplot2::aes(shape = sig), size = 3.2,
                        colour = "grey15", fill = "white") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f, p = %s", est, p_fmt)),
                       vjust = -1.25, size = 3.1, colour = "grey25") +
    ggplot2::scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21),
                                labels = c(`TRUE` = "p < .05", `FALSE` = "n.s."),
                                name = NULL, drop = FALSE) +
    ggplot2::facet_wrap(~ model_lab, ncol = 1, scales = "free_y") +
    ## Per-panel padding, same reasoning as figure (3).
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = c(0.55, 0.95))) +
    ggplot2::labs(
      title    = "Slope differences between drugs",
      subtitle = sprintf("Interaction terms: %s minus %s, raw units, with 95%% CI",
                         MREG_GROUPS[2], MREG_GROUPS[1]),
      x = sprintf("Difference in %s change per unit of predictor", mreg_lab(MREG_Y)),
      y = NULL,
      caption = paste0("A CI clear of the dashed line means the two drugs genuinely differ on that slope. ",
                       "Raw units: the log10 ROE difference is not an effect size,\nsince ROE spans a ",
                       "very different range in the two groups.")
    ) +
    mreg_theme()

  ## Sized like figure (3): tall enough for the densest panel, times the
  ## number of panels.
  mreg_save(p_gd, "MREG_group_difference_forest.png",
            w = 9, h = length(unique(idf$model_lab)) *
                       max(1.6, 0.30 * max(table(idf$model_lab))) + 2.6,
            dir = PREG_DIR)
}


## ---- 20f. THE POOLED REGRESSION TABLE AS A FIGURE ----------------------------
## Section 18h's table figure, called with the pooled coefficient/fit tables. It
## filters on the `group` column, which is why every row above carries
## PREG_GROUP_LAB -- no change to that function was needed beyond letting the
## caller override its caption, since the default one advertises standardized
## betas this section deliberately does not produce.

mreg_table_figure(
  PREG_GROUP_LAB, PREG_MODELS, preg_coef_tbl, preg_fit_tbl, MREG_Y,
  dir = PREG_DIR, file = "PREG_table_pooled.png", shorts = PREG_MODEL_SHORT,
  caption = paste0("*** p < .001, ** p < .01, * p < .05.  Drug is a factor with ",
                   MREG_GROUPS[1], " as reference, so every 'drug' row reads ",
                   MREG_GROUPS[2], " minus ", MREG_GROUPS[1], ".\n",
                   "In P2 the plain receptor and ROE rows are that predictor's slope in ",
                   MREG_GROUPS[1], "; add the matching interaction row for ", MREG_GROUPS[2],
                   " (both drugs' slopes are tabulated in the .xlsx).\n",
                   "High VIF on the interaction rows is expected -- they are products of the ",
                   "model's own main effects -- and affects the main-effect SEs, not the ",
                   "interaction tests or the nested F."))


## ---- 20g. SELF-CHECK ---------------------------------------------------------
## A model with every predictor interacted with group is algebraically identical
## to fitting the two groups separately. So the slopes derived in 20d must equal
## section 18's stratified HT4_HT6 coefficients exactly. If they ever stop
## matching, the term renaming in 20c or the vcov arithmetic in 20d is wrong --
## and the failure would otherwise be silent, because both tables would still
## look perfectly reasonable on their own.

if (nrow(preg_group_slopes) && exists("mreg_coef_tbl")) {
  chk <- preg_group_slopes %>%
    dplyr::select(predictor, group, derived = est) %>%
    dplyr::left_join(
      mreg_coef_tbl %>%
        dplyr::filter(model == MREG_FULL) %>%
        dplyr::select(predictor = term, group, stratified = est),
      by = c("predictor", "group")) %>%
    dplyr::mutate(abs_diff = abs(derived - stratified))

  worst <- suppressWarnings(max(chk$abs_diff, na.rm = TRUE))
  cat("\n----- SELF-CHECK: derived slopes vs section 18's stratified fits -----\n\n")
  print(as.data.frame(chk), row.names = FALSE, digits = 6)
  if (!is.finite(worst)) {
    warning("09: self-check could not run -- no matching stratified coefficients.")
  } else if (worst > 1e-8) {
    warning("09: derived group slopes DISAGREE with section 18's stratified fits ",
            "(max |diff| = ", format(worst), "). The term renaming in 20c or the ",
            "vcov arithmetic in 20d is wrong.")
  } else {
    cat("\nOK: max |difference| =", format(worst, digits = 3),
        "-- the pooled interaction model reproduces the stratified fits exactly.\n")
  }
}


## ---- 20h. EXPORT -------------------------------------------------------------

preg_config <- tibble::tibble(
  setting = c("outcome", "exposure", "groups", "drug_coding", "models",
              "n_rule", "scaling", "vif_note", "provenance", "output_dir"),
  value = c(
    MREG_Y, MREG_X, paste(MREG_GROUPS, collapse = " | "),
    paste0("factor, reference = ", MREG_GROUPS[1],
           ". Every 'drug' coefficient reads ", MREG_GROUPS[2], " MINUS ",
           MREG_GROUPS[1], ". The analysis plan specified 1 = ", MREG_GROUPS[2],
           ", 2 = ", MREG_GROUPS[1], ", i.e. the opposite direction; which drug ",
           "is the reference changes only the SIGN (identical fitted values, R2, ",
           "F, SEs and p-values), and matching section 18's direction was ",
           "preferred over matching the coding."),
    paste(sprintf("%s: %s ~ %s", names(PREG_MODELS), MREG_Y,
                  vapply(PREG_MODELS, paste, character(1), collapse = " + ")),
          collapse = "  |  "),
    paste0("complete cases AND all values finite (drops ROE == 0, whose log10 is ",
           "-Inf). Pooled n is the two stratified n's added together."),
    paste0("raw units only. No standardized betas: these models are fit exactly ",
           "once (in 20b) and a z-scored refit would break that, and a ",
           "standardized product term is not the product of the standardized ",
           "terms, so the number would not mean what beta_std means in section 18."),
    paste0("Model 2's interaction columns are products of its own main effects ",
           "and so are collinear with them by construction. High VIF there ",
           "inflates the MAIN-EFFECT standard errors; it does not affect the ",
           "interaction tests or the nested F, which are what this section is for."),
    paste0("Coefficients are read off the models fit in 20b -- the only place in ",
           "the pipeline these are fit. Nothing here is refit, so the ",
           "coefficient table and the interaction table cannot disagree."),
    PREG_DIR
  )
)

preg_sheets <- list(
  Pooled_coefficients = preg_coef_tbl %>%
    dplyr::select(group, model, model_lab, term, term_lab, est, se, t, p, p_fmt,
                  ci_lo, ci_hi, vif),
  Pooled_model_fit    = preg_fit_tbl,
  Nested_test         = preg_nested_tbl,
  Derived_group_slopes = if (nrow(preg_group_slopes)) preg_group_slopes else
    tibble::tibble(note = "no interaction model fit"),
  Interaction_terms   = if (nrow(mreg_interaction_tbl)) mreg_interaction_tbl else
    tibble::tibble(note = "no interaction models fit"),
  Interaction_omnibus = if (nrow(mreg_intomni_tbl)) mreg_intomni_tbl else
    tibble::tibble(note = "none"),
  Config              = preg_config
)

writexl::write_xlsx(preg_sheets,
                    file.path(PREG_DIR, "Stats_Between_Drug_Comparison.xlsx"))

message("\nSection 20 (between-drug comparison) written to: ",
        normalizePath(PREG_DIR, mustWork = FALSE))


## ---- 20j. THE SAME COMPARISON, WITH SOWS AS THE OUTCOME ----------------------
## Mirrors 07's section 18l: everything above (20a-20h) tests NRS pain; this
## block reruns the identical drug-pooled M_add/M_int comparison with SOWS
## (withdrawal) as the outcome. Same fitting function as 20b (preg_fit_all())
## and the same slope arithmetic as 20d (preg_group_slopes_for()) -- SOWS gets
## its own call to each, not a second copy of either, so the two outcomes
## cannot drift apart the way two independently-written fitting loops would.
##
## Scope is deliberately narrower than 20a-20h: the coefficient table, the
## nested F, each drug's derived slope, and the pooled table figure -- the
## same four things 18l reports for the stratified SOWS models, but pooled.
## No marginal-effect / interaction-forest figures (20e) and no self-check
## (20g): both of those exist because NRS's stratified fits (mreg_coef_tbl)
## are already sitting in this session from section 18 and used across
## several figures; SOWS's stratified fits (mreg_sows, from 18l) are not
## reused by anything else, so duplicating that machinery here for a
## consistency check nothing else depends on would be scope the request did
## not ask for.
##
## SOWS cannot be a predictor of itself, so it is stripped from every model
## (a no-op today, since no model here contains SOWS -- kept as a guard, same
## reasoning and same no-op status as 18l's MREG_MODELS_SOWS).
##
## TWO model-definition lists, not one, and they are NOT interchangeable:
##   MREG_MODELS_FOR_SOWS  real column names (log_ROE, R_5HT4, ...) -- what
##                         preg_fit_all() needs to build a formula and index
##                         `master`. Same shape as MREG_MODELS_SOWS in 18l.
##   PREG_MODELS_SOWS      the POST-FIT renamed term vocabulary
##                         (group_TvsH, group_x_5HT4, ...) that
##                         preg_clean_terms() produces once a model is
##                         already fit. These names are not columns of
##                         `master` -- passing this to preg_fit_all() finds no
##                         matching columns and silently fits nothing, which
##                         is what happened before this comment was added.
##                         It is only for mreg_table_figure()'s `models`
##                         argument (20f uses PREG_MODELS the same way),
##                         which reads row labels out of the ALREADY-RENAMED
##                         coefficient table.
PREG_Y_SOWS           <- "SOWS"
MREG_MODELS_FOR_SOWS  <- lapply(MREG_MODELS, setdiff, PREG_Y_SOWS)
PREG_MODELS_SOWS      <- lapply(PREG_MODELS, setdiff, PREG_Y_SOWS)

PREG_DIR_SOWS <- file.path(OUT_DIR, "between_drug_comparison_SOWS")
dir.create(PREG_DIR_SOWS, showWarnings = FALSE, recursive = TRUE)

cat("\n================================================================\n")
cat("20j. BETWEEN-DRUG COMPARISON: ", PREG_Y_SOWS,
    " ~ drug x (log ROE + 5-HT4 / 5-HT6)\n", sep = "")
cat("     All coefficients read ", MREG_GROUPS[2], " MINUS ", MREG_GROUPS[1], ".\n", sep = "")
cat("================================================================\n")

preg_sows_fit <- preg_fit_all(PREG_Y_SOWS, MREG_MODELS_FOR_SOWS)
preg_sows_fi  <- preg_sows_fit$inter_fits[[MREG_FULL]]

if (is.null(preg_sows_fi)) {

  warning("09 (20j): no pooled SOWS fit for '", MREG_FULL, "' -- section skipped.")

} else {

  preg_sows_desc <- list(
    Pooled_main   = preg_describe(preg_sows_fi$add, "Pooled_main"),
    Pooled_grpint = preg_describe(preg_sows_fi$int, "Pooled_grpint")
  )

  preg_sows_coef_tbl <- dplyr::bind_rows(lapply(preg_sows_desc, `[[`, "coef")) %>%
    dplyr::mutate(model_lab = unname(PREG_MODEL_LABELS[model]),
                  term_lab  = ifelse(term == "(Intercept)", "Intercept",
                                     vapply(term, mreg_lab, character(1))),
                  p_fmt     = mreg_fmt_p(p))
  preg_sows_fit_tbl <- dplyr::bind_rows(lapply(preg_sows_desc, `[[`, "fit"))

  preg_sows_nested_tbl <- preg_sows_fit$intomni_tbl %>%
    dplyr::filter(model == MREG_FULL) %>%
    dplyr::mutate(reduced = "Pooled_main", full = "Pooled_grpint",
                  added_terms = "drug x each predictor",
                  p_fmt = mreg_fmt_p(p)) %>%
    dplyr::select(reduced, full, added_terms, n, df_num, df_den, F_stat, p, p_fmt)

  cat("\n----- POOLED SOWS MODEL COEFFICIENTS (raw units) -----\n\n")
  print(as.data.frame(preg_sows_coef_tbl %>%
          dplyr::select(model, term, est, se, ci_lo, ci_hi, t, p, vif)),
        row.names = FALSE, digits = 3)

  cat("\n----- POOLED SOWS MODEL FIT -----\n\n")
  print(as.data.frame(preg_sows_fit_tbl %>%
          dplyr::select(model, n, R2, adj_R2, F_stat, df1, df2, F_p,
                        RMSE, AIC, max_VIF, shapiro_p, bp_p)),
        row.names = FALSE, digits = 3)

  cat("\n----- DO THE DRUGS DIFFER ON SOWS? (nested F, Model 1 vs Model 2) -----\n\n")
  print(as.data.frame(preg_sows_nested_tbl), row.names = FALSE, digits = 3)

  preg_sows_group_slopes <- preg_group_slopes_for(preg_sows_fi)

  if (nrow(preg_sows_group_slopes)) {
    cat("\n----- EACH DRUG'S SOWS SLOPE, DERIVED FROM THE INTERACTION MODEL -----\n")
    cat("(not separately fit: main effect for ", MREG_GROUPS[1],
        ", main + interaction for ", MREG_GROUPS[2], ")\n\n", sep = "")
    print(as.data.frame(preg_sows_group_slopes %>%
            dplyr::select(predictor, group, est, se, ci_lo, ci_hi, t, p)),
          row.names = FALSE, digits = 3)
  }

  mreg_table_figure(
    PREG_GROUP_LAB, PREG_MODELS_SOWS, preg_sows_coef_tbl, preg_sows_fit_tbl,
    PREG_Y_SOWS, dir = PREG_DIR_SOWS, file = "PREG_SOWS_table_pooled.png",
    shorts = PREG_MODEL_SHORT,
    caption = paste0("*** p < .001, ** p < .01, * p < .05.  Drug is a factor with ",
                     MREG_GROUPS[1], " as reference, so every 'drug' row reads ",
                     MREG_GROUPS[2], " minus ", MREG_GROUPS[1], ".\n",
                     "In P2 the plain receptor and ROE rows are that predictor's slope in ",
                     MREG_GROUPS[1], "; add the matching interaction row for ", MREG_GROUPS[2],
                     " (both drugs' slopes are tabulated in the .xlsx)."))

  writexl::write_xlsx(
    list(
      Pooled_coefficients = preg_sows_coef_tbl %>%
        dplyr::select(group, model, model_lab, term, term_lab, est, se, t, p, p_fmt,
                      ci_lo, ci_hi, vif),
      Pooled_model_fit     = preg_sows_fit_tbl,
      Nested_test          = preg_sows_nested_tbl,
      Derived_group_slopes = if (nrow(preg_sows_group_slopes)) preg_sows_group_slopes else
        tibble::tibble(note = "no interaction model fit"),
      Config = tibble::tibble(
        setting = c("outcome", "exposure", "groups", "models", "note", "output_dir"),
        value = c(
          PREG_Y_SOWS, MREG_X, paste(MREG_GROUPS, collapse = " | "),
          paste(sprintf("%s: %s ~ %s", names(PREG_MODELS_SOWS), PREG_Y_SOWS,
                        vapply(PREG_MODELS_SOWS, paste, character(1), collapse = " + ")),
                collapse = "  |  "),
          paste0("Same pooled model set as the NRS comparison (section 20) with SOWS ",
                 "moved from a covariate-free right-hand side to the outcome. Table and ",
                 "numbers only -- see the NRS output folder for the marginal-effect and ",
                 "interaction-forest figures, which are NRS-specific."),
          PREG_DIR_SOWS
        )
      )
    ),
    file.path(PREG_DIR_SOWS, "Stats_Between_Drug_Comparison_SOWS.xlsx"))

  message("Section 20j (SOWS) written to: ",
          normalizePath(PREG_DIR_SOWS, mustWork = FALSE))
}


## ---- 20i. NOTES --------------------------------------------------------------
## - The omnibus F in 20b is the number that licenses section 18's whole design.
##   If it is significant, reporting one pooled slope per predictor would be
##   wrong and the stratified tables are the correct presentation. If it is not,
##   the pooled Model 1 is the better summary and the stratified split is
##   spending degrees of freedom for nothing. Read it before reading anything
##   else in either section.
##
## - A non-significant omnibus F is NOT evidence the drugs behave identically.
##   Model 2 spends 8 parameters on n = 36, and interaction terms are the
##   lowest-powered thing in a regression. Read a null here as "this sample
##   cannot resolve a difference", not as "there is none".
##
## - Model 2's main effects and Model 1's coefficients answer different
##   questions and should not be compared row by row. Model 1's R_5HT4 is a
##   single slope assumed common to both drugs; Model 2's R_5HT4 is the slope in
##   Hydrocodone alone. They differ for the same reason the two are different
##   models, not because anything is inconsistent.
##
## - Raw coefficients are not comparable ACROSS predictors here any more than in
##   section 18 -- log10 ROE, 5-HT4 and 5-HT6 are on different scales. Unlike
##   section 18 there are no standardized betas to fall back on (see 20c), so
##   compare a predictor only against itself between the two drugs.
##
## - The log10 ROE slope difference in particular is not an effect size. ROE
##   spans ~40x different ranges in the two groups, so the two drugs' ROE slopes
##   are estimated over barely overlapping x. The test of equality is still
##   valid; the magnitude is not interpretable as "how much bigger".
