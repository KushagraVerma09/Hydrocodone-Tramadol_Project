################################################################################
## 18. MULTIPLE REGRESSION: NRS ~ log ROE + 5-HT4 / 5-HT6
##     Nine pre-specified models, fit independently in each drug group.
##     Section 18l repeats the regression table with SOWS as the outcome.
##
## MODELS (each fit separately within Hydrocodone group and Tramadol group).
## Model KEYS are on the left; every model contains log_ROE.
##   HT4                NRS = b0 + b_ROE*log_ROE + b_HT4*R_5HT4                + e
##   HT6                NRS = b0 + b_ROE*log_ROE                + b_HT6*R_5HT6 + e
##   HT4_ROEx           HT4     + b_X4*(log_ROE * R_5HT4)                      + e
##   HT6_ROEx           HT6     + b_X6*(log_ROE * R_5HT6)                      + e
##   HT4_HT6            NRS = b0 + b_ROE*log_ROE + b_HT4*R_5HT4 + b_HT6*R_5HT6 + e
##   HT4_HT6_Int        HT4_HT6 + b_I*(R_5HT4 * R_5HT6)                        + e
##   HT4_HT6_ROEx       HT4_HT6 + b_X4*(log_ROE * R_5HT4)
##                              + b_X6*(log_ROE * R_5HT6)                      + e
##   HT4_HT6_Cov        HT4_HT6 + Sex + MQS non-opioid + log10 MME + SOWS + DOU
##   HT4_HT6_Cov_noMME  HT4_HT6 + Sex + MQS non-opioid +             SOWS + DOU
##
## WHY THE ROE x RECEPTOR MODELS EXIST
##   HT4 and HT6 assume the receptor-pain slope is the SAME at every level of
##   opioid exposure -- that 5-HT4 relates to pain identically in a lightly
##   exposed patient and a heavily exposed one. The ROEx models drop that
##   assumption by adding the product of log ROE and the receptor: a
##   moderation test, asking whether exposure CHANGES the receptor-pain
##   relationship rather than merely sitting alongside it.
##
##   Three of them, deliberately. HT4_ROEx and HT6_ROEx add one product each,
##   so they nest cleanly under HT4 and HT6 as 1-df tests; HT4_HT6_ROEx adds
##   both at once on top of HT4_HT6, a 2-df test of "does exposure moderate
##   EITHER receptor". At n = 17-19 the 1-df tests are the better-powered
##   reading and the 2-df one is the omnibus.
##
##   READ THE PRODUCT TERM AND THE NESTED F, NOT THE MAIN EFFECTS. As in
##   HT4_HT6_Int the products are raw, uncentered (the SCALING note below is
##   why), so in a ROEx model b_HT4 becomes "5-HT4's association with pain
##   when log_ROE = 0", i.e. at ROE = 1 mg/L -- roughly 50x the largest ROE
##   any subject in this data has. That is an extrapolation to a point no
##   patient occupies, and its standard error is correspondingly large.
##   Centering would move those main effects back into the observed range but
##   would NOT change the product coefficient, its p-value, or the nested F --
##   which is the whole result -- so nothing is lost by leaving them raw and
##   reading only what is interpretable. Expect high VIFs on the main effects
##   in these models for the same reason; that is collinearity between a term
##   and its own product, not a data problem.
##
## WHY THE 5-HT4 x 5-HT6 INTERACTION MODEL EXISTS
##   HT4_HT6 asks whether 5-HT4 and 5-HT6 carry SEPARATE information about
##   pain. It does not ask whether their effects depend on EACH OTHER --
##   whether 5-HT4's association with pain is stronger or weaker at different
##   levels of 5-HT6. That second question is a statistical interaction
##   (moderation) term: literally the product of the two receptor values,
##   added to HT4_HT6 as a fourth predictor. Note that in that model b_HT4 and
##   b_HT6 become "this receptor's effect when the OTHER receptor = 0" -- a
##   value that never occurs in this data -- so read the interaction
##   coefficient and the nested F test, not the two main effects.
##
## WHY THERE ARE TWO ADJUSTED MODELS
##   HT4_HT6_Cov asks whether the receptor-pain associations survive
##   adjustment for the obvious confounders: sex, non-opioid medication load
##   (MQS non-opioid), opioid dose (log10 MME), withdrawal symptoms (SOWS) and
##   how long the patient has been on opioids (DOU). Of those five, four are
##   fully observed in both drug groups -- but MME is missing for 5 subjects
##   (1 hydrocodone, 4 tramadol), and complete-case fitting would silently
##   drop them. HT4_HT6_Cov_noMME is the same model without MME, fit on the
##   same rows as every other model here, so the cost of including MME is
##   VISIBLE (compare the two n's) instead of hidden. Because the two are fit
##   on different rows, their R2 and AIC are NOT comparable to each other --
##   see 18k.
##
##   Read both adjusted models as exploratory. HT4_HT6_Cov is 9 parameters on
##   n ~ 18 (hydrocodone) and n ~ 13 (tramadol), and the tramadol group is
##   20/23 female, so its sex coefficient in particular is barely identified.
##   Directions and the model-level fit are worth reading; individual
##   covariate CIs are too wide to be more than a lead.
##
## HOW THIS DIFFERS FROM SECTION 16
##   Section 16 runs the SAME two-predictor model once per receptor across the
##   whole 19-receptor panel -- a screen. At n = 17-19 nothing survives FDR
##   there, and its two figures live entirely in coefficient space (19 rows of
##   overlapping CIs, and a scatter of sign flips) without ever showing a data
##   point. This section is the confirmatory counterpart: two receptors chosen
##   in advance on pharmacological grounds, reported the way a multiple
##   regression is normally reported -- coefficient table, added-variable
##   plots, observed vs predicted, residual diagnostics, model comparison --
##   and then a formal test of whether the slopes differ between drugs.
##   HT4 and HT6 are numerically identical to section 16's R_5HT4 and R_5HT6
##   rows (same formula, same data); that is a deliberate cross-check.
##
## WHY HT4_HT6 IS THE INTERESTING ONE
##   HT4 and HT6 each ask "does this receptor track pain once exposure is
##   netted out". HT4_HT6 asks the harder question: do 5-HT4 and 5-HT6 carry
##   SEPARATE information about pain, or is one of them just standing in for
##   the other? If both coefficients survive in HT4_HT6 at roughly their
##   single-receptor values, they are independent signals. If one collapses
##   toward zero, the two receptors were sharing variance and only one of them
##   is doing work. 18d tests this formally with nested F tests.
##
## WHY SEPARATELY BY GROUP
##   Same reasoning as sections 15-17: the ROE-pain slope reverses sign between
##   drugs, so a pooled model averages the two into something smaller than
##   either. One model per group is the stratified estimate; 18e then tests
##   whether the slopes actually differ rather than leaving it to the eye.
##
## SCALING -- RAW UNITS
##   Unlike sections 15/16/17, coefficients here are reported in RAW units:
##   NRS points per unit of predictor. This is the clinically interpretable
##   scale. Two consequences to keep in mind, both restated in the figure
##   captions:
##     (a) raw coefficients are NOT comparable ACROSS predictors -- log ROE,
##         5-HT4 and 5-HT6 are on different scales, so a larger number does not
##         mean a larger effect. Standardized betas are carried alongside in
##         the xlsx for exactly that comparison.
##     (b) the raw log ROE slope is not directly comparable ACROSS drugs
##         either, since ROE differs ~40x in scale between them. The
##         interaction test in 18e is still a valid test of raw-unit slope
##         equality; just do not read the two raw numbers as effect sizes.
##   Covariate scales in the two adjusted models: MME enters as log10 (it is
##   heavily right-skewed and no subject has MME == 0, so nothing is lost),
##   DOU enters RAW in years (so its coefficient reads as NRS points per extra
##   year on opioids), and sex enters as a 0/1 indicator, sex_female.
##
## THE -Inf TRAP
##   log_ROE = log10(ROE) is -Inf for the 4 subjects per group with ROE == 0,
##   and complete.cases() does NOT drop -Inf. Every model frame below is
##   filtered with is.finite(), matching sections 16 and 17. Expect n = 19
##   (hydrocodone) and n = 17 (tramadol).
##
## Self-contained: needs `master`, `FOCUS_GROUPS`, `OUT_DIR`, `fmt_p`,
## `custom_colors` and `MED_LABELS` (all from 00_config.R/01_load_data.R), and
## the packages loaded there (dplyr, tidyr, tibble, purrr, ggplot2, writexl).
## No new packages -- the 2x2 diagnostic panels use base graphics on a png()
## device precisely so that patchwork/cowplot are not needed.
##
## Run standalone: Rscript R/07_nested_regression.R
## Sourced by: Run_All.R.
################################################################################

.HTK_PROJECT_DIR <- "/Users/kushagraverma/Work/Projects/Hydrocodone+Tramadol_Project"
## Always re-sourced -- see the note above this line in 01_load_data.R.
source(file.path(.HTK_PROJECT_DIR, "R", "00_config.R"))
if (!exists("master"))              source(file.path(.HTK_PROJECT_DIR, "R", "01_load_data.R"))

## ---- 18a. CONFIG -------------------------------------------------------------

MREG_Y      <- "nrs"          # outcome
MREG_X      <- "log_ROE"      # exposure, present in every model
MREG_GROUPS <- FOCUS_GROUPS   # Hydrocodone group / Tramadol group
MREG_MIN_N  <- 12             # skip a group x model cell smaller than this

## Derived predictors, computed once here on `master` so that 18c's per-group
## fits, 18e's pooled fit and 18l's SOWS models all see the identical columns.
##
## (a) The interaction term: the product of 5-HT4 and 5-HT6 activity.
master$R_5HT4_x_5HT6 <- master$R_5HT4 * master$R_5HT6

## (a2) The exposure x receptor interactions: log ROE times each receptor.
## These ask a different question from (a) -- see WHY THE ROE x RECEPTOR MODELS
## EXIST in the header.
##
## log_ROE is -Inf for the ROE == 0 subjects, so these products are -Inf (or
## NaN where the receptor value is exactly 0). That is harmless and needs no
## special handling: log_ROE is itself a predictor in every model, so 18b's
## is.finite() filter has already dropped those same rows before the product is
## ever looked at. The n of these models therefore matches HT4_HT6 exactly,
## which is what makes the nested F tests in 18d valid.
master$log_ROE_x_5HT4 <- master$log_ROE * master$R_5HT4
master$log_ROE_x_5HT6 <- master$log_ROE * master$R_5HT6

## (b) Sex as a numeric 0/1 indicator. `master$female` (01_load_data.R:100) is
## LOGICAL, and handing lm() a logical makes it name the coefficient
## "femaleTRUE" rather than "female" -- which silently breaks every place
## downstream that matches a coefficient name against a predictor name (the VIF
## lookup in 18b, the added-variable plots in 18g, the partial-residual
## contributions in 18i, and the regression table's row ordering). A numeric
## column keeps coefficient name == column name. 1 = female, 0 = male.
master$sex_female <- as.numeric(master$female)

## The nine models. Each entry is the full predictor vector, in the order the
## terms should appear. Change this list (and MREG_Y) and re-source from 18a to
## re-run the whole section on different predictors -- nothing downstream is
## hard-coded to 5-HT4 / 5-HT6 or to how many models there are.
##
## MME enters as log10, DOU raw -- see the SCALING note in the header. The two
## adjusted models differ ONLY in whether log_MME is present: MME is the one
## new predictor with missing values (5 subjects), so the pair makes the cost
## of adjusting for dose visible in the n column rather than hiding it.
MREG_MODELS <- list(
  HT4               = c("log_ROE", "R_5HT4"),
  HT6               = c("log_ROE", "R_5HT6"),
  HT4_ROEx          = c("log_ROE", "R_5HT4", "log_ROE_x_5HT4"),
  HT6_ROEx          = c("log_ROE", "R_5HT6", "log_ROE_x_5HT6"),
  HT4_HT6           = c("log_ROE", "R_5HT4", "R_5HT6"),
  HT4_HT6_Int       = c("log_ROE", "R_5HT4", "R_5HT6", "R_5HT4_x_5HT6"),
  HT4_HT6_ROEx      = c("log_ROE", "R_5HT4", "R_5HT6",
                        "log_ROE_x_5HT4", "log_ROE_x_5HT6"),
  HT4_HT6_Cov       = c("log_ROE", "R_5HT4", "R_5HT6", "sex_female",
                        "MQS_nonopioid", "log_MME", "SOWS", "DOU"),
  HT4_HT6_Cov_noMME = c("log_ROE", "R_5HT4", "R_5HT6", "sex_female",
                        "MQS_nonopioid", "SOWS", "DOU")
)

## Full model labels -- each one names its own predictors, so a reader looking
## at a facet strip or a table header never has to come back up here to find
## out what the model contains.
MREG_MODEL_LABELS <- c(
  HT4               = "ROE + 5-HT4",
  HT6               = "ROE + 5-HT6",
  HT4_ROEx          = "ROE + 5-HT4 + ROEx5-HT4",
  HT6_ROEx          = "ROE + 5-HT6 + ROEx5-HT6",
  HT4_HT6           = "ROE + 5-HT4 + 5-HT6",
  HT4_HT6_Int       = "ROE + 5-HT4 + 5-HT6 + 5-HT4x5-HT6",
  HT4_HT6_ROEx      = "ROE + 5-HT4 + 5-HT6 + ROEx5-HT4 + ROEx5-HT6",
  HT4_HT6_Cov       = "ROE + 5-HT4 + 5-HT6 + Sex + MQS non-opioid + log MME + SOWS + DOU",
  HT4_HT6_Cov_noMME = "ROE + 5-HT4 + 5-HT6 + Sex + MQS non-opioid + SOWS + DOU"
)

## Nested comparisons to run in 18d: c(reduced, full).
##
## HT4_HT6 -> HT4_HT6_Cov is deliberately ABSENT. log_MME has missing values,
## so the adjusted-with-MME model is fit on fewer rows than HT4_HT6, and a
## nested F test across different row sets is not valid -- 18d's own guard
## would refuse it anyway. HT4_HT6 -> HT4_HT6_Cov_noMME IS valid: sex,
## MQS non-opioid, SOWS and DOU are fully observed in both drug groups, so
## that model lands on exactly the same rows as HT4_HT6.
##
## The three ROE x receptor comparisons are all valid: the product columns are
## built from predictors already in the reduced model, so they add no new
## missingness and each pair lands on identical rows. HT4 -> HT4_ROEx and
## HT6 -> HT6_ROEx are 1-df tests; HT4_HT6 -> HT4_HT6_ROEx is 2 df, and at
## n = 17-19 the two 1-df tests are the better-powered way to ask the question.
MREG_NESTED <- list(c("HT4", "HT4_HT6"), c("HT6", "HT4_HT6"),
                    c("HT4", "HT4_ROEx"), c("HT6", "HT6_ROEx"),
                    c("HT4_HT6", "HT4_HT6_Int"),
                    c("HT4_HT6", "HT4_HT6_ROEx"),
                    c("HT4_HT6", "HT4_HT6_Cov_noMME"))

MREG_DIR <- file.path(OUT_DIR, "multiple_regression_HT4_HT6")
dir.create(MREG_DIR, showWarnings = FALSE, recursive = TRUE)

## 18l writes its SOWS tables here, kept separate from the NRS outputs above.
MREG_DIR_SOWS <- file.path(OUT_DIR, "multiple_regression_HT4_HT6_SOWS")
dir.create(MREG_DIR_SOWS, showWarnings = FALSE, recursive = TRUE)

MREG_COLS <- if (exists("custom_colors")) custom_colors else
  setNames(c("#1b9e77", "#d95f02")[seq_along(MREG_GROUPS)], MREG_GROUPS)

## Labels for the predictors this section adds that the shared MED_LABELS
## (00_config.R) does not know about. Kept LOCAL rather than appended to
## MED_LABELS because this file advertises itself as self-contained and
## sections 15-17 have no use for them. Consulted before MED_LABELS below.
MREG_EXTRA_LABELS <- c(
  ## Overrides MED_LABELS' "5-HT4 x 5-HT6 (raw interaction)": with the
  ## mean-centered variant gone there is nothing for "raw" to contrast with.
  R_5HT4_x_5HT6  = "5-HT4 x 5-HT6",
  log_ROE_x_5HT4 = "log10 ROE x 5-HT4",
  log_ROE_x_5HT6 = "log10 ROE x 5-HT6",
  sex_female    = "Sex (female = 1)",
  MQS_nonopioid = "MQS non-opioid",
  log_MME       = "log10 MME",
  SOWS          = "SOWS withdrawal",
  DOU           = "Duration of opioid use (yr)"
)

## Variable label lookup. Robust to keys that are in neither table (unlike
## med_lab at section 15), same shape as section 17's pmed_lab.
mreg_lab <- function(v) {
  own <- MREG_EXTRA_LABELS[v]
  if (length(own) == 1 && !is.na(own)) return(unname(own))
  if (!exists("MED_LABELS")) return(v)
  hit <- MED_LABELS[v]
  if (length(hit) == 1 && !is.na(hit)) unname(hit) else v
}

## NOTE: the script-level fmt_p (section 2, ~line 294) takes a SINGLE p and is
## built out of bare if() calls, so handing it a column raises
##   "the condition has length > 1".
## Everything here formats whole columns, so wrap it. Same situation section 17
## handles with pmed_fmt_p. Do not "simplify" these into one function.
mreg_fmt_p <- function(p) vapply(p, function(x) fmt_p(x), character(1))

## Significance stars for the table figure. No p < .10 tier: the conventional
## "." marker sits right after a decimal number ("-10.613.") and reads as a
## typo.
mreg_stars <- function(p) {
  ifelse(is.na(p), "",
  ifelse(p < 0.001, "***",
  ifelse(p < 0.01,  "**",
  ifelse(p < 0.05,  "*", ""))))
}

## Short axis-friendly labels used on the coefficient forests.
mreg_short <- function(v) {
  out <- sub("^R_", "", v)
  out[out == "5HT4_x_5HT6"]   <- "5-HT4 x 5-HT6"
  out[out == "log_ROE_x_5HT4"] <- "log10 ROE x 5-HT4"
  out[out == "log_ROE_x_5HT6"] <- "log10 ROE x 5-HT6"
  out[out == "log_ROE"]       <- "log10 ROE"
  out[out == "log_MME"]       <- "log10 MME"
  out[out == "sex_female"]    <- "Sex (female)"
  out[out == "MQS_nonopioid"] <- "MQS non-opioid"
  out[out == "DOU"]           <- "DOU (yr)"
  out
}

## Short MODEL labels, for the two spots with only enough room for one model
## per column/bar rather than one per row (the regression-table header (4) and
## the model-comparison x-axis (6)). MREG_MODEL_LABELS's full formula reads
## fine one at a time in a title/subtitle, but nine of them side by side overlap
## into unreadable text. These are keyed the same way as MREG_MODEL_LABELS and
## pre-wrapped by hand -- every model here contains ROE, so the short forms
## only need to say what is added beyond it. The full predictor list is still
## visible in the table's own row labels and in every facet strip.
MREG_MODEL_SHORT <- c(
  HT4               = "ROE +\n5-HT4",
  HT6               = "ROE +\n5-HT6",
  HT4_ROEx          = "ROE + 5-HT4\n+ ROE x 5-HT4",
  HT6_ROEx          = "ROE + 5-HT6\n+ ROE x 5-HT6",
  HT4_HT6           = "ROE + 5-HT4\n+ 5-HT6",
  HT4_HT6_Int       = "+ 5-HT4 x\n5-HT6",
  HT4_HT6_ROEx      = "+ ROE x 5-HT4,\nROE x 5-HT6",
  HT4_HT6_Cov       = "+ Sex, MQS,\nlog MME,\nSOWS, DOU",
  HT4_HT6_Cov_noMME = "+ Sex, MQS,\nSOWS, DOU"
)

## 18l drops SOWS from the two adjusted models, so it needs short labels that
## do not advertise a predictor those models no longer contain. Spelled out
## rather than regex-stripped from the strings above -- they are hand-wrapped,
## and deleting "SOWS, " from a wrapped string leaves the line breaks wrong.
MREG_MODEL_SHORT_SOWS <- MREG_MODEL_SHORT
MREG_MODEL_SHORT_SOWS["HT4_HT6_Cov"]       <- "+ Sex, MQS,\nlog MME, DOU"
MREG_MODEL_SHORT_SOWS["HT4_HT6_Cov_noMME"] <- "+ Sex, MQS,\nDOU"

## Look a short label up BY MODEL KEY, falling back to the key itself so a
## model added to MREG_MODELS without a short-label entry still draws.
mreg_model_short <- function(keys, tbl = MREG_MODEL_SHORT) {
  out <- unname(tbl[keys])
  ifelse(is.na(out), keys, out)
}


## ---- 18b. FIT ONE MODEL ------------------------------------------------------

## Fits `yv ~ preds` on one group's data and returns everything the figures and
## tables downstream need. Returns NULL (with a message) rather than erroring if
## the cell is too small or a predictor is constant.
mreg_fit_one <- function(df, preds, yv, model_key, group_lab) {

  vars <- c(yv, preds)
  miss <- setdiff(vars, names(df))
  if (length(miss)) {
    message("  [", group_lab, " / ", model_key, "] missing column(s): ",
            paste(miss, collapse = ", "), " -- skipped")
    return(NULL)
  }

  ## as.data.frame, not the tibble `master` hands over: tibbles silently drop
  ## row names, and the PIN labels below are what make the outliers that base R
  ## flags in the diagnostic panels name an actual subject instead of a row
  ## number.
  d0 <- as.data.frame(df[, vars, drop = FALSE])
  n0 <- nrow(d0)
  if ("PIN" %in% names(df)) rownames(d0) <- as.character(df$PIN)

  ## Two-stage filter, and report the two reasons separately: complete.cases()
  ## alone would keep the ROE == 0 subjects, whose log_ROE is -Inf.
  ok_cc  <- stats::complete.cases(d0)
  d1     <- d0[ok_cc, , drop = FALSE]
  ok_fin <- apply(d1, 1, function(r) all(is.finite(r)))
  d      <- d1[ok_fin, , drop = FALSE]

  n_miss <- n0 - nrow(d1)
  n_inf  <- nrow(d1) - nrow(d)
  if (n_miss || n_inf)
    message("  [", group_lab, " / ", model_key, "] dropped ", n_miss,
            " row(s) with missing values and ", n_inf,
            " row(s) with non-finite values (log10 of ROE == 0)")

  n <- nrow(d)
  if (n < MREG_MIN_N) {
    message("  [", group_lab, " / ", model_key, "] n = ", n, " < MREG_MIN_N (",
            MREG_MIN_N, ") -- skipped")
    return(NULL)
  }

  sds <- vapply(d, stats::sd, numeric(1))
  dead <- names(sds)[!is.finite(sds) | sds < 1e-12]
  if (length(dead)) {
    message("  [", group_lab, " / ", model_key, "] no variance in: ",
            paste(dead, collapse = ", "), " -- skipped")
    return(NULL)
  }

  f <- stats::as.formula(sprintf("`%s` ~ %s", yv,
                                 paste(sprintf("`%s`", preds), collapse = " + ")))
  m <- stats::lm(f, data = d)
  s <- summary(m)

  ## ---- raw-unit coefficient table (the primary result) ----
  cf  <- s$coefficients
  ci  <- suppressWarnings(stats::confint(m))
  nms <- rownames(cf)

  coef_tbl <- tibble::tibble(
    group    = group_lab,
    model    = model_key,
    term     = gsub("`", "", nms),
    est      = cf[, 1],
    se       = cf[, 2],
    t        = cf[, 3],
    p        = cf[, 4],
    ci_lo    = ci[match(nms, rownames(ci)), 1],
    ci_hi    = ci[match(nms, rownames(ci)), 2]
  )

  ## ---- standardized refit, reference only ----
  ## z-scored within this cell's own complete-case sample, so the betas are
  ## directly comparable across predictors and to sections 15/16/17.
  ds <- d
  for (v in names(ds)) ds[[v]] <- (ds[[v]] - mean(ds[[v]])) / stats::sd(ds[[v]])
  m_std <- stats::lm(f, data = ds)
  cs    <- summary(m_std)$coefficients
  std_map <- setNames(cs[, 1], gsub("`", "", rownames(cs)))
  coef_tbl$beta_std <- unname(std_map[coef_tbl$term])

  ## ---- VIF, computed by hand (no `car` dependency) ----
  ## VIF_j = 1 / (1 - R^2_j), where R^2_j comes from regressing predictor j on
  ## every other predictor. With one predictor there is nothing to inflate.
  vifs <- setNames(rep(NA_real_, length(preds)), preds)
  if (length(preds) > 1) {
    for (j in preds) {
      others <- setdiff(preds, j)
      fj  <- stats::as.formula(sprintf("`%s` ~ %s", j,
                                       paste(sprintf("`%s`", others), collapse = " + ")))
      r2j <- summary(stats::lm(fj, data = d))$r.squared
      vifs[j] <- if (is.finite(r2j) && r2j < 1) 1 / (1 - r2j) else Inf
    }
  }
  coef_tbl$vif <- unname(vifs[coef_tbl$term])

  ## ---- model-level fit statistics ----
  fst <- s$fstatistic
  f_p <- if (!is.null(fst)) stats::pf(fst[1], fst[2], fst[3], lower.tail = FALSE) else NA_real_
  res <- stats::residuals(m)

  ## Shapiro-Wilk on residuals (normality) and a hand-rolled Breusch-Pagan
  ## (constant variance): regress squared residuals on the predictors, then
  ## n * R^2 is chi-square with (#predictors) df under homoskedasticity.
  sw_p <- tryCatch(stats::shapiro.test(res)$p.value, error = function(e) NA_real_)
  bp_p <- tryCatch({
    daux <- d; daux$.u2 <- res^2
    faux <- stats::as.formula(paste(".u2 ~", paste(sprintf("`%s`", preds), collapse = " + ")))
    r2a  <- summary(stats::lm(faux, data = daux))$r.squared
    stats::pchisq(n * r2a, df = length(preds), lower.tail = FALSE)
  }, error = function(e) NA_real_)

  fit_tbl <- tibble::tibble(
    group      = group_lab,
    model      = model_key,
    formula    = paste(deparse(f), collapse = " "),
    n          = n,
    k          = length(preds),
    R2         = s$r.squared,
    adj_R2     = s$adj.r.squared,
    F_stat     = if (!is.null(fst)) unname(fst[1]) else NA_real_,
    df1        = if (!is.null(fst)) unname(fst[2]) else NA_real_,
    df2        = if (!is.null(fst)) unname(fst[3]) else NA_real_,
    F_p        = unname(f_p),
    sigma      = s$sigma,
    RMSE       = sqrt(mean(res^2)),
    AIC        = stats::AIC(m),
    BIC        = stats::BIC(m),
    max_VIF    = if (all(is.na(vifs))) NA_real_ else max(vifs, na.rm = TRUE),
    max_cooksD = max(stats::cooks.distance(m), na.rm = TRUE),
    shapiro_p  = sw_p,
    bp_p       = bp_p
  )

  list(model = m, model_std = m_std, data = d, preds = preds, yv = yv,
       coef = coef_tbl, fit = fit_tbl, key = model_key, group = group_lab)
}


## ---- 18c. FIT EVERY GROUP x MODEL CELL ---------------------------------------

## Fit an entire model set in every group and assemble the two tidy tables the
## rest of the section reads. Factored out of the loop it used to be so that
## 18l can re-run the same machinery with SOWS as the outcome without
## duplicating any of it. `labels` is the key -> display-label vector to stamp
## onto model_lab (18l passes a SOWS-specific one, since dropping SOWS from the
## predictor list changes what the adjusted models are called).
mreg_fit_all <- function(models, yv, groups = MREG_GROUPS,
                         labels = MREG_MODEL_LABELS) {
  fits <- list(); coef_rows <- list(); fit_rows <- list()

  for (g in groups) {
    d_g <- master[!is.na(master$Plot_Group) & master$Plot_Group == g, , drop = FALSE]
    message("Group: ", g, "  (", nrow(d_g), " subjects before filtering)")
    fits[[g]] <- list()

    for (mk in names(models)) {
      res <- mreg_fit_one(d_g, models[[mk]], yv, mk, g)
      if (is.null(res)) next
      fits[[g]][[mk]]           <- res
      coef_rows[[paste(g, mk)]] <- res$coef
      fit_rows[[paste(g, mk)]]  <- res$fit
    }
  }

  coef_tbl <- dplyr::bind_rows(coef_rows)
  fit_tbl  <- dplyr::bind_rows(fit_rows)

  if (!nrow(coef_tbl))
    stop("Section 18: no models could be fit for outcome '", yv,
         "' -- check the model list and `master`.")

  coef_tbl <- coef_tbl %>%
    dplyr::mutate(
      model_lab = unname(labels[model]),
      term_lab  = mreg_short(term),
      sig       = !is.na(p) & p < 0.05,
      p_fmt     = mreg_fmt_p(p)
    )

  fit_tbl <- fit_tbl %>%
    dplyr::mutate(model_lab = unname(labels[model]),
                  F_p_fmt   = mreg_fmt_p(F_p))

  list(fits = fits, coef = coef_tbl, fit = fit_tbl)
}

cat("\n================================================================\n")
cat("18. MULTIPLE REGRESSION: ", MREG_Y, " ~ log ROE + 5-HT4 / 5-HT6\n", sep = "")
cat("================================================================\n")

mreg_all      <- mreg_fit_all(MREG_MODELS, MREG_Y)
mreg_fits     <- mreg_all$fits
mreg_coef_tbl <- mreg_all$coef
mreg_fit_tbl  <- mreg_all$fit

cat("\n----- COEFFICIENTS (raw units: NRS points per unit predictor) -----\n")
print(as.data.frame(mreg_coef_tbl %>%
        dplyr::filter(term != "(Intercept)") %>%
        dplyr::select(group, model, term, est, se, ci_lo, ci_hi, t, p, beta_std, vif)),
      row.names = FALSE, digits = 3)

cat("\n----- MODEL FIT -----\n")
print(as.data.frame(mreg_fit_tbl %>%
        dplyr::select(group, model, n, R2, adj_R2, F_stat, df1, df2, F_p,
                      RMSE, AIC, max_VIF, shapiro_p, bp_p)),
      row.names = FALSE, digits = 3)


## ---- 18d. NESTED MODEL TESTS -------------------------------------------------
## Does the second receptor add anything over the first? HT4 and HT6 are each
## nested inside HT4_HT6, so an F test on the residual sums of squares is the right
## comparison -- but ONLY when both models were fit on the same rows. They are
## here (the same is.finite filter, and neither receptor has extra missingness),
## which is checked explicitly below rather than assumed.

mreg_nested_rows <- list()

for (g in names(mreg_fits)) {
  for (pair in MREG_NESTED) {
    red <- mreg_fits[[g]][[pair[1]]]
    ful <- mreg_fits[[g]][[pair[2]]]
    if (is.null(red) || is.null(ful)) next
    if (nrow(red$data) != nrow(ful$data)) {
      message("  [", g, "] ", pair[1], " vs ", pair[2],
              " fit on different n (", nrow(red$data), " vs ", nrow(ful$data),
              ") -- nested F test not valid, skipped")
      next
    }
    av <- stats::anova(red$model, ful$model)
    added <- setdiff(ful$preds, red$preds)

    ## Keyed on group + BOTH model names, not just the reduced one: MREG_NESTED
    ## has two pairs that share the same reduced model (HT4_HT6 -> HT4_HT6_Int
    ## and HT4_HT6 -> HT4_HT6_Cov_noMME), and a key of just paste(g, pair[1]) would
    ## make the second overwrite the first in this list.
    mreg_nested_rows[[paste(g, pair[1], pair[2])]] <- tibble::tibble(
      group        = g,
      reduced      = pair[1],
      full         = pair[2],
      added_term   = paste(added, collapse = " + "),
      n            = nrow(ful$data),
      df_num       = av$Df[2],
      df_den       = av$Res.Df[2],
      F_stat       = av$F[2],
      p            = av$`Pr(>F)`[2],
      R2_reduced   = summary(red$model)$r.squared,
      R2_full      = summary(ful$model)$r.squared,
      delta_R2     = summary(ful$model)$r.squared - summary(red$model)$r.squared,
      adjR2_reduced = summary(red$model)$adj.r.squared,
      adjR2_full    = summary(ful$model)$adj.r.squared,
      delta_adjR2   = summary(ful$model)$adj.r.squared - summary(red$model)$adj.r.squared,
      delta_AIC     = stats::AIC(ful$model) - stats::AIC(red$model)
    )
  }
}

mreg_nested_tbl <- dplyr::bind_rows(mreg_nested_rows)
if (nrow(mreg_nested_tbl)) {
  mreg_nested_tbl$p_fmt <- mreg_fmt_p(mreg_nested_tbl$p)
  cat("\n----- NESTED MODEL TESTS (does the second receptor add anything?) -----\n")
  cat("(negative delta_AIC favours the fuller model)\n\n")
  print(as.data.frame(mreg_nested_tbl %>%
          dplyr::select(group, reduced, full, added_term, n, F_stat, df_num, df_den,
                        p, delta_R2, delta_adjR2, delta_AIC)),
        row.names = FALSE, digits = 3)
}


## ---- 18e. GROUP DIFFERENCE: INTERACTION MODELS -------------------------------
## Fitting the two drugs separately shows the slopes look different; it does not
## test that they ARE different. Pool the two groups and interact every
## predictor with drug: each interaction coefficient IS the raw-unit slope
## difference (tramadol minus hydrocodone, given the factor's level order), and
## the model-level F test asks whether the drugs differ on any slope at all.

mreg_pool <- master %>%
  dplyr::filter(!is.na(Plot_Group), Plot_Group %in% MREG_GROUPS) %>%
  dplyr::mutate(.grp = factor(as.character(Plot_Group), levels = MREG_GROUPS))

mreg_inter_rows  <- list()
mreg_intomni_rows <- list()
mreg_inter_fits  <- list()

for (mk in names(MREG_MODELS)) {
  preds <- MREG_MODELS[[mk]]
  vars  <- c(MREG_Y, preds, ".grp")
  if (!all(setdiff(vars, ".grp") %in% names(mreg_pool))) next

  dp <- as.data.frame(mreg_pool[, vars, drop = FALSE])
  dp <- dp[stats::complete.cases(dp), , drop = FALSE]
  num <- setdiff(vars, ".grp")
  dp  <- dp[apply(dp[, num, drop = FALSE], 1, function(r) all(is.finite(r))), , drop = FALSE]
  if (nrow(dp) < MREG_MIN_N || nlevels(droplevels(dp$.grp)) < 2) next

  rhs  <- paste(sprintf("`%s`", preds), collapse = " + ")
  m_add <- stats::lm(stats::as.formula(sprintf("`%s` ~ .grp + %s", MREG_Y, rhs)), data = dp)
  m_int <- stats::lm(stats::as.formula(sprintf("`%s` ~ .grp * (%s)", MREG_Y, rhs)), data = dp)
  mreg_inter_fits[[mk]] <- list(add = m_add, int = m_int, data = dp, preds = preds)

  ci <- suppressWarnings(stats::confint(m_int))
  cf <- summary(m_int)$coefficients
  keep <- grep("^\\.grp.*:", rownames(cf))

  if (length(keep)) {
    trm <- rownames(cf)[keep]
    mreg_inter_rows[[mk]] <- tibble::tibble(
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
  mreg_intomni_rows[[mk]] <- tibble::tibble(
    model  = mk, n = nrow(dp),
    df_num = av$Df[2], df_den = av$Res.Df[2],
    F_stat = av$F[2],  p = av$`Pr(>F)`[2]
  )
}

mreg_interaction_tbl <- dplyr::bind_rows(mreg_inter_rows)
mreg_intomni_tbl     <- dplyr::bind_rows(mreg_intomni_rows)

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


## ---- 18f. FIGURE THEME + FILE-NAME HELPERS -----------------------------------
## Shared look for every ggplot in this section, so the 20-odd files read as one
## set. Matches the "regime B" style of sections 15-17 (theme_minimal, no minor
## grid, bold strips, grey caveat caption).

mreg_theme <- function(base = 12) {
  ggplot2::theme_minimal(base_size = base) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.text       = ggplot2::element_text(face = "bold"),
      plot.title       = ggplot2::element_text(face = "bold"),
      plot.caption     = ggplot2::element_text(hjust = 0, colour = "grey40", size = 8),
      legend.position  = "bottom"
    )
}

## Filesystem-safe slug: "Hydrocodone group" -> "Hydrocodonegroup", matching the
## convention section 17's path-diagram filenames already use.
mreg_slug <- function(x) gsub("[^A-Za-z0-9]", "", x)

## `dir` defaults to the NRS output folder; 18l passes MREG_DIR_SOWS.
mreg_save <- function(p, file, w, h, dir = MREG_DIR) {
  ggplot2::ggsave(file.path(dir, file), p, width = w, height = h,
                  dpi = 300, bg = "white")
}

## Wrapped by hand: ggplot does not wrap captions, and an over-long line is
## silently clipped at the right edge of the device.
CAP_RAW <- paste0("Coefficients are in raw units (NRS points per unit of predictor).\n",
                  "Raw coefficients are NOT comparable across predictors -- log10 ROE, ",
                  "5-HT4 and 5-HT6 are on different scales;\n",
                  "use the standardized betas in the accompanying .xlsx to compare magnitudes.")


## ---- 18g. FIGURES, ONE PER GROUP x MODEL -------------------------------------

## (1) The classic 2x2 regression diagnostic panel. Deliberately base graphics:
## plot.lm() is the panel every reader already knows how to read, and it needs
## no extra package (patchwork/cowplot are not installed).
for (g in names(mreg_fits)) {
  for (mk in names(mreg_fits[[g]])) {
    res <- mreg_fits[[g]][[mk]]
    f <- file.path(MREG_DIR, sprintf("MREG_diag_%s_%s.png", mk, mreg_slug(g)))
    grDevices::png(f, width = 2400, height = 2000, res = 220, bg = "white")
    op <- graphics::par(mfrow = c(2, 2), mar = c(4.2, 4.2, 3.2, 1.6),
                        oma = c(0, 0, 3.2, 0), cex = 0.85)
    ## sub.caption = "" suppresses plot.lm's own "stats::lm(f)" strapline, which
    ## would otherwise print on top of the mtext title below.
    graphics::plot(res$model, which = 1:4, col = unname(MREG_COLS[g]),
                   pch = 19, cex = 1.1, sub.caption = "")
    graphics::mtext(sprintf("%s  --  %s   (n = %d)",
                            unname(MREG_MODEL_LABELS[mk]), g, nrow(res$data)),
                    outer = TRUE, cex = 1.15, font = 2, line = 0.6)
    graphics::par(op)
    grDevices::dev.off()
  }
}

## (2) Added-variable (partial regression) plots -- the figure that actually
## shows the DATA behind a multiple regression, and the thing section 16's
## figures never had. For predictor j: residuals of y on all OTHER predictors,
## against residuals of x_j on all other predictors. The slope of the line
## through that cloud IS the model's coefficient for x_j, so a coefficient with
## no visible trend here is a coefficient you should not believe.
for (g in names(mreg_fits)) {
  for (mk in names(mreg_fits[[g]])) {
    res   <- mreg_fits[[g]][[mk]]
    d     <- res$data
    preds <- res$preds
    if (length(preds) < 2) next   # with one predictor an AV plot is just the scatter

    av_rows <- list()
    for (j in preds) {
      others <- setdiff(preds, j)
      ry <- stats::resid(stats::lm(
        stats::as.formula(sprintf("`%s` ~ %s", res$yv,
                                  paste(sprintf("`%s`", others), collapse = " + "))), data = d))
      rx <- stats::resid(stats::lm(
        stats::as.formula(sprintf("`%s` ~ %s", j,
                                  paste(sprintf("`%s`", others), collapse = " + "))), data = d))
      cr <- res$coef[res$coef$term == j, ]
      av_rows[[j]] <- tibble::tibble(
        predictor = j,
        panel = sprintf("%s\nb = %.3f (SE %.3f), p = %s",
                        mreg_lab(j), cr$est[1], cr$se[1], mreg_fmt_p(cr$p[1])),
        rx = rx, ry = ry, sig = !is.na(cr$p[1]) && cr$p[1] < 0.05
      )
    }
    av <- dplyr::bind_rows(av_rows)
    av$panel <- factor(av$panel, levels = unique(av$panel))

    p_av <- ggplot2::ggplot(av, ggplot2::aes(rx, ry)) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey85") +
      ggplot2::geom_vline(xintercept = 0, colour = "grey85") +
      ggplot2::geom_smooth(ggplot2::aes(linetype = sig), method = "lm", formula = y ~ x,
                           se = TRUE, colour = unname(MREG_COLS[g]),
                           fill = unname(MREG_COLS[g]), alpha = 0.15, linewidth = 1) +
      ggplot2::geom_point(size = 2.6, colour = unname(MREG_COLS[g]), alpha = 0.85) +
      ggplot2::scale_linetype_manual(values = c(`TRUE` = "solid", `FALSE` = "22"),
                                     guide = "none") +
      ggplot2::facet_wrap(~ panel, scales = "free_x", nrow = 1) +
      ggplot2::labs(
        title    = sprintf("Added-variable plots -- %s", g),
        subtitle = sprintf("%s   (n = %d).  Each panel: %s and the predictor, both with the OTHER predictors removed",
                           unname(MREG_MODEL_LABELS[mk]), nrow(d), mreg_lab(res$yv)),
        x = "Predictor residual (other predictors partialled out)",
        y = sprintf("%s residual", mreg_lab(res$yv)),
        caption = paste0("The slope of each line is that predictor's coefficient in the multiple regression. ",
                         "Solid = p < .05, dashed = n.s.\nBand is the 95% CI of the partial fit.")
      ) +
      mreg_theme()

    mreg_save(p_av, sprintf("MREG_avplot_%s_%s.png", mk, mreg_slug(g)),
              w = 4.6 * length(preds) + 1.2, h = 5.6)
  }
}


## ---- 18h. FIGURES, ONE PER GROUP ---------------------------------------------

## (3) Coefficient forest: every term of every model for one drug, so you can
## watch the log ROE coefficient move (or not) as receptors enter, and see
## whether 5-HT4 / 5-HT6 hold their single-receptor values once they are in
## HT4_HT6 together, and whether either survives the adjusted models.
for (g in names(mreg_fits)) {
  cf <- mreg_coef_tbl %>%
    dplyr::filter(group == g, term != "(Intercept)") %>%
    dplyr::mutate(
      model_lab = factor(model_lab, levels = unname(MREG_MODEL_LABELS)),
      term_lab  = factor(term_lab,
                         levels = rev(unique(mreg_short(unlist(MREG_MODELS)))))
    )
  if (!nrow(cf)) next

  p_cf <- ggplot2::ggplot(cf, ggplot2::aes(x = est, y = term_lab)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = ci_lo, xmax = ci_hi),
                            orientation = "y", width = 0.16, linewidth = 0.8,
                            colour = unname(MREG_COLS[g])) +
    ggplot2::geom_point(ggplot2::aes(shape = sig), size = 3.2,
                        colour = unname(MREG_COLS[g]), fill = "white") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f, p = %s", est, p_fmt)),
                       vjust = -1.25, size = 3.1, colour = "grey25") +
    ggplot2::scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21),
                                labels = c(`TRUE` = "p < .05", `FALSE` = "n.s."),
                                name = NULL, drop = FALSE) +
    ggplot2::facet_wrap(~ model_lab, ncol = 1, scales = "free_y") +
    ## Pad each panel by a fixed amount ABOVE ITS OWN TOP ROW (0.9 of a row, to
    ## clear the value labels, which sit above their point). The previous
    ## expand_limits() reserved room for the FULL term factor in every panel,
    ## which was harmless at 4 terms but now squashes the two-term panels into
    ## the bottom third of their facet while the eight-term ones stay cramped.
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = c(0.55, 0.95))) +
    ggplot2::labs(
      title    = sprintf("Regression coefficients -- %s", g),
      subtitle = sprintf("%s predicted from log10 ROE and receptor activity; estimate with 95%% CI",
                         mreg_lab(MREG_Y)),
      x = sprintf("Change in %s per unit of predictor (raw units)", mreg_lab(MREG_Y)),
      y = NULL, caption = CAP_RAW
    ) +
    mreg_theme()

  ## facet_wrap gives every panel the SAME height regardless of how many rows
  ## it holds (space = "free_y" is a facet_grid feature), so the figure has to
  ## be tall enough for the DENSEST panel -- the 8-term adjusted model -- and
  ## then multiplied by the panel count. The 1.6 floor keeps a set of
  ## two-term models from collapsing into a stripe.
  n_panels <- length(unique(cf$model_lab))
  max_rows <- max(table(cf$model_lab))
  mreg_save(p_cf, sprintf("MREG_coefficients_%s.png", mreg_slug(g)),
            w = 8.6, h = n_panels * max(1.6, 0.30 * max_rows) + 2.6)
}

## (4) The regression table as a figure -- the "Table 2" of a regression paper,
## every model side by side. Drawn as a ggplot text grid rather than gt::gtsave,
## which would need chromote/webshot2 (not installed).
##
## Written as a FUNCTION rather than an inline loop because 18l reproduces this
## one figure (and only this one) with SOWS as the outcome. Arguments: the
## group, the model set it was fit with, that group's coefficient and fit
## tables, the outcome variable name (for the title), and where to write.
mreg_table_figure <- function(g, models, coef_tbl, fit_tbl, yv,
                              dir = MREG_DIR, file = NULL,
                              shorts = MREG_MODEL_SHORT) {
  mods <- intersect(names(models), unique(fit_tbl$model[fit_tbl$group == g]))
  if (!length(mods)) return(invisible(NULL))
  if (is.null(file)) file <- sprintf("MREG_table_%s.png", mreg_slug(g))

  terms_order <- c("(Intercept)", unique(unlist(models)))
  cf <- coef_tbl %>% dplyr::filter(group == g)

  body_rows <- list()
  for (tm in terms_order) {
    if (!any(cf$term == tm)) next
    cells <- vapply(mods, function(mk) {
      r <- cf[cf$group == g & cf$model == mk & cf$term == tm, ]
      if (!nrow(r)) return("--")
      sprintf("%.3f%s\n(%.3f)", r$est[1], mreg_stars(r$p[1]), r$se[1])
    }, character(1))
    lbl <- if (tm == "(Intercept)") "Intercept" else mreg_lab(tm)
    body_rows[[tm]] <- c(lbl, cells)
  }

  ft <- fit_tbl %>% dplyr::filter(group == g)
  stat_row <- function(lbl, fn) c(lbl, vapply(mods, function(mk) {
    r <- ft[ft$model == mk, ]; if (!nrow(r)) "--" else fn(r)
  }, character(1)))

  foot_rows <- list(
    stat_row("n",       function(r) sprintf("%d", r$n[1])),
    stat_row("R2",      function(r) sprintf("%.3f", r$R2[1])),
    stat_row("adj. R2", function(r) sprintf("%.3f", r$adj_R2[1])),
    stat_row("F (df)",  function(r) sprintf("%.2f (%d, %d)", r$F_stat[1], r$df1[1], r$df2[1])),
    stat_row("model p", function(r) mreg_fmt_p(r$F_p[1])),
    stat_row("RMSE",    function(r) sprintf("%.3f", r$RMSE[1])),
    stat_row("AIC",     function(r) sprintf("%.1f", r$AIC[1])),
    stat_row("max VIF", function(r) if (is.na(r$max_VIF[1])) "--" else sprintf("%.2f", r$max_VIF[1]))
  )

  all_rows <- c(body_rows, foot_rows)
  n_body   <- length(body_rows)
  ncol_tot <- length(mods) + 1

  cellsdf <- do.call(rbind, lapply(seq_along(all_rows), function(i) {
    data.frame(row = i, col = seq_len(ncol_tot),
               txt = as.character(all_rows[[i]]),
               is_stat = i > n_body, stringsAsFactors = FALSE)
  }))
  ## y descends down the page; header sits above row 1.
  cellsdf$y <- -cellsdf$row
  cellsdf$x <- cellsdf$col

  ## Header cells are the SHORT model labels, looked up by model key.
  hdr <- data.frame(x = seq_len(ncol_tot), y = 0,
                    txt = c("", mreg_model_short(mods, shorts)),
                    stringsAsFactors = FALSE)

  p_tab <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = c(0.55, -0.45, -(n_body + 0.45),
                                       -(length(all_rows) + 0.55)),
                        colour = "grey35", linewidth = 0.4) +
    ggplot2::geom_text(data = hdr, ggplot2::aes(x, y, label = txt),
                       fontface = "bold", size = 3.5, lineheight = 0.95) +
    ggplot2::geom_text(data = subset(cellsdf, col == 1),
                       ggplot2::aes(x, y, label = txt, fontface = ifelse(is_stat, "italic", "plain")),
                       hjust = 0, nudge_x = -0.55, size = 3.4, lineheight = 0.95) +
    ggplot2::geom_text(data = subset(cellsdf, col > 1),
                       ggplot2::aes(x, y, label = txt),
                       size = 3.4, lineheight = 0.95) +
    ggplot2::scale_x_continuous(limits = c(0.4, ncol_tot + 0.6), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(-(length(all_rows) + 1.1), 1.6), expand = c(0, 0)) +
    ggplot2::labs(
      title    = sprintf("%s regressed on log10 ROE and receptor activity -- %s",
                         mreg_lab(yv), g),
      subtitle = "Raw-unit coefficient, standard error in parentheses",
      caption  = paste0("*** p < .001, ** p < .01, * p < .05.  ",
                        "Raw coefficients are not comparable across predictors (different scales); ",
                        "see the .xlsx for standardized betas.\n",
                        "VIF > 5 would indicate predictors competing for the same variance",
                        " -- EXCEPT in the interaction models, where a product term is",
                        " collinear with its own components by construction and a high VIF is",
                        " expected, not a fault (it inflates the SEs of the main effects, not",
                        " of the product term or the nested F).\n",
                        "Columns differ in n where a predictor has missing values -- compare the n row.")
    ) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", hjust = 0.5, margin = ggplot2::margin(b = 4)),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, colour = "grey30",
                                            margin = ggplot2::margin(b = 10)),
      plot.caption  = ggplot2::element_text(hjust = 0, colour = "grey40", size = 8,
                                            margin = ggplot2::margin(t = 10)),
      plot.margin   = ggplot2::margin(14, 16, 12, 16),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA)
    )

  ## 2.3 rather than the original 2.6 per column: with nine models the table is
  ## 10 columns wide, and the old factor pushed it past 19 inches even at six.
  mreg_save(p_tab, file, dir = dir,
            w = 2.3 * ncol_tot + 1.4, h = 0.38 * length(all_rows) + 2.4)
  invisible(p_tab)
}

for (g in names(mreg_fits))
  mreg_table_figure(g, MREG_MODELS, mreg_coef_tbl, mreg_fit_tbl, MREG_Y)


## ---- 18i. FIGURES ACROSS GROUPS ----------------------------------------------

## (5) Observed vs predicted. The most honest single picture of how much of NRS
## these models actually explain: points on the diagonal = good prediction,
## a flat cloud = the model knows nothing.
op_rows <- list()
for (g in names(mreg_fits)) for (mk in names(mreg_fits[[g]])) {
  res <- mreg_fits[[g]][[mk]]
  ft  <- mreg_fit_tbl[mreg_fit_tbl$group == g & mreg_fit_tbl$model == mk, ]
  ## SHORT model label here, not the full formula: facet_grid gives each model
  ## a ~3-inch column and the adjusted models' full labels run to 60+
  ## characters, which the strip clips without warning. The full formula is on
  ## the coefficient forest and the regression table.
  op_rows[[paste(g, mk)]] <- tibble::tibble(
    group = g, model = mk,
    model_lab = mreg_model_short(mk),
    observed  = res$data[[res$yv]],
    predicted = unname(stats::fitted(res$model)),
    lab = sprintf("R2 = %.3f\nadj R2 = %.3f\nRMSE = %.2f\nn = %d",
                  ft$R2[1], ft$adj_R2[1], ft$RMSE[1], ft$n[1])
  )
}
op_df <- dplyr::bind_rows(op_rows)

if (nrow(op_df)) {
  op_df$model_lab <- factor(op_df$model_lab,
                            levels = mreg_model_short(names(MREG_MODELS)))
  op_df$group     <- factor(op_df$group, levels = MREG_GROUPS)
  ann <- op_df %>% dplyr::distinct(group, model_lab, lab) %>%
    dplyr::mutate(x = min(op_df$predicted, op_df$observed),
                  y = max(op_df$predicted, op_df$observed))

  p_op <- ggplot2::ggplot(op_df, ggplot2::aes(predicted, observed, colour = group)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey55") +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                         linewidth = 0.8, alpha = 0.9) +
    ggplot2::geom_point(size = 2.6, alpha = 0.85) +
    ggplot2::geom_label(data = ann, ggplot2::aes(x, y, label = lab),
                        inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3,
                        fill = "white", alpha = 0.8, linewidth = 0,
                        lineheight = 0.95) +
    ggplot2::scale_colour_manual(values = MREG_COLS, name = NULL) +
    ggplot2::facet_grid(group ~ model_lab) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title    = sprintf("Observed vs predicted %s", mreg_lab(MREG_Y)),
      subtitle = "Dotted line = perfect prediction; coloured line = the model's own fit",
      x = sprintf("Predicted %s", mreg_lab(MREG_Y)),
      y = sprintf("Observed %s", mreg_lab(MREG_Y)),
      caption = paste0("R2 is in-sample and rises mechanically with every predictor added -- ",
                       "compare adj R2 (and the nested tests) across models, not R2.")
    ) +
    mreg_theme()

  ## Width scales with how many model columns facet_grid actually draws.
  ## 3.2 rather than 3.5 per column now that there are nine of them.
  mreg_save(p_op, "MREG_obs_vs_pred.png",
            w = 3.2 * length(unique(op_df$model_lab)) + 2, h = 8.6)
}

## (6) Model comparison: adjusted R2 and AIC side by side, with the nested F
## tests annotated. adj R2 and AIC both penalise the extra parameter, so if the
## two-receptor model does not beat the single-receptor ones on these it is not
## earning its keep.
##
## CAUTION: HT4_HT6_Cov is fit on FEWER rows than every other
## bar here (log_MME is missing for 5 subjects), and R2/AIC computed on
## different row sets are not comparable. Its bar is legitimate as a
## description of that model; it is not a like-for-like contest against the
## others. HT4_HT6_Cov_noMME is the bar to compare against HT4_HT6.
if (nrow(mreg_fit_tbl)) {
  comp <- mreg_fit_tbl %>%
    dplyr::select(group, model, n, R2, adj_R2, AIC) %>%
    tidyr::pivot_longer(c(R2, adj_R2, AIC), names_to = "metric", values_to = "value") %>%
    dplyr::mutate(
      metric = factor(metric, levels = c("R2", "adj_R2", "AIC"),
                      labels = c("R2 (in-sample)", "Adjusted R2", "AIC (lower is better)")),
      ## x is the SHORT label, keyed off the model name -- nine full formulas
      ## side by side on one axis is unreadable.
      model_lab = factor(mreg_model_short(model),
                         levels = mreg_model_short(names(MREG_MODELS))),
      group = factor(group, levels = MREG_GROUPS)
    )

  ## Nested-F results as the subtitle, ONE LINE PER TEST -- with 4 nested
  ## comparisons now (up from 2), grouping every test for a group onto a
  ## single line runs off the right edge of the device and gets clipped.
  sub_txt <- NULL
  if (nrow(mreg_nested_tbl)) {
    per_row <- sprintf("%s:  +%s over %s: F(%d,%d) = %.2f, p = %s",
                       mreg_nested_tbl$group, mreg_short(mreg_nested_tbl$added_term),
                       mreg_nested_tbl$reduced, mreg_nested_tbl$df_num,
                       mreg_nested_tbl$df_den, mreg_nested_tbl$F_stat,
                       mreg_fmt_p(mreg_nested_tbl$p))
    sub_txt <- paste0("Nested F tests (does the added term improve the fit?)\n",
                      paste(per_row, collapse = "\n"))
  }

  p_comp <- ggplot2::ggplot(comp, ggplot2::aes(model_lab, value, fill = group)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.78), width = 0.72,
                      alpha = 0.9) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", value)),
                       position = ggplot2::position_dodge(width = 0.78),
                       vjust = -0.45, size = 3.1, colour = "grey25") +
    ggplot2::scale_fill_manual(values = MREG_COLS, name = NULL) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.14))) +
    ggplot2::facet_wrap(~ metric, scales = "free_y", nrow = 1) +
    ggplot2::labs(
      title    = "Model comparison",
      subtitle = sub_txt,
      x = NULL, y = NULL,
      caption = paste0("R2 can only rise when a predictor is added; adjusted R2 and AIC ",
                       "penalise the extra parameter, so they are the ones to compare.\n",
                       "Only models fit on the SAME rows are comparable here: the adjusted ",
                       "model with log MME loses the subjects whose MME is missing, so read it ",
                       "on its own, not against the others.")
    ) +
    mreg_theme() +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 8.5, colour = "grey30"),
                   axis.text.x = ggplot2::element_text(size = 9))

  ## Width scales with the number of models (each gets its own bar per group
  ## per metric panel); 2.6 per model rather than 3.5, since even six models at
  ## the old factor made a 24-inch-wide figure. Height scales with how many
  ## nested-test lines are in the subtitle.
  mreg_save(p_comp, "MREG_model_comparison.png",
            w = 2.6 * length(MREG_MODELS) + 3,
            h = 5.5 + 0.3 * (nrow(mreg_nested_tbl) + 1))
}

## (7) Marginal (adjusted) effect of each receptor from the model HT4_HT6:
## predicted NRS across the observed range of that receptor with the other
## predictors held at their group means, 95% CI ribbon, raw data overlaid as
## partial residuals so the points sit on the same scale as the line.
## HT4_HT6, the two-receptor main-effects model, is used below for the
## per-receptor marginal-effect plots (7) and the drug-interaction-slopes plot
## (8). Deliberately hardcoded rather than "last model in the list", for two
## separate reasons:
##
##   - HT4_HT6_Int contains a product term. Both plots hold every OTHER
##     predictor at its mean while sweeping one receptor across the x-axis,
##     which is wrong when "other predictors" includes a product built from
##     the very receptor being swept -- the product would stay frozen at its
##     mean instead of moving with the receptor, so the drawn line would not
##     match what the model actually predicts.
##   - The two adjusted models would freeze six or seven covariates at their
##     means to draw one receptor's line. That line is not wrong, but it is a
##     picture of an 8-predictor model at n ~ 13-19 and reads as far more
##     precise than it is.
##
## The other models still get their own coefficient forest, regression table,
## nested test and model-comparison entries; they just aren't the model behind
## THESE two figures.
MREG_FULL <- "HT4_HT6"

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

  mreg_save(p_mar, sprintf("MREG_partial_%s.png", sub("^R_", "", rv)), w = 9, h = 7)
}

## (8) The interaction picture: each drug's fitted line for each predictor, from
## the full interaction model, with the slope-difference p annotated. This is
## the visual form of the test in 18e.
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
              w = 4.6 * length(fi$preds) + 1.2, h = 6.2)
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
                       max(1.6, 0.30 * max(table(idf$model_lab))) + 2.6)
}


## ---- 18j. EXPORT -------------------------------------------------------------

mreg_config <- tibble::tibble(
  setting = c("outcome", "exposure", "groups", "models", "interaction_raw",
              "sex_coding", "covariate_scales", "min_n", "scaling", "n_rule",
              "output_dir"),
  value = c(
    MREG_Y, MREG_X, paste(MREG_GROUPS, collapse = " | "),
    paste(sprintf("%s: %s ~ %s", names(MREG_MODELS), MREG_Y,
                  vapply(MREG_MODELS, paste, character(1), collapse = " + ")),
          collapse = "  |  "),
    "R_5HT4_x_5HT6 = R_5HT4 * R_5HT6  (HT4_HT6_Int)",
    "sex_female = 1 if female (sex == 2), 0 if male; from master$female",
    paste0("MME entered as log10 (log_MME); DOU entered RAW in years; ",
           "MQS_nonopioid and SOWS entered raw. HT4_HT6_Cov is fit on fewer rows ",
           "than every other model because log_MME has missing values -- its R2/AIC ",
           "are NOT comparable to the others; use HT4_HT6_Cov_noMME for that."),
    as.character(MREG_MIN_N),
    "raw units primary; standardized betas (z within group) in the Coefficients sheet",
    "complete cases AND all values finite (drops ROE == 0, whose log10 is -Inf)",
    MREG_DIR
  )
)

mreg_sheets <- list(
  Coefficients      = mreg_coef_tbl %>%
    dplyr::select(group, model, model_lab, term, term_lab, est, se, t, p, p_fmt,
                  ci_lo, ci_hi, beta_std, vif),
  Model_fit         = mreg_fit_tbl,
  Nested_tests      = if (nrow(mreg_nested_tbl)) mreg_nested_tbl else tibble::tibble(note = "no nested comparisons run"),
  Interaction_terms = if (nrow(mreg_interaction_tbl)) mreg_interaction_tbl else tibble::tibble(note = "no interaction models fit"),
  Interaction_omnibus = if (nrow(mreg_intomni_tbl)) mreg_intomni_tbl else tibble::tibble(note = "none"),
  Diagnostics       = mreg_fit_tbl %>%
    dplyr::select(group, model, n, max_VIF, max_cooksD, shapiro_p, bp_p, sigma, RMSE),
  Config            = mreg_config
)

writexl::write_xlsx(mreg_sheets,
                    file.path(MREG_DIR, "Stats_Multiple_Regression_HT4_HT6_NRS.xlsx"))

message("\nSection 18 results and figures written to: ",
        normalizePath(MREG_DIR, mustWork = FALSE))


## ---- 18l. THE SAME REGRESSION TABLE, WITH SOWS AS THE OUTCOME ----------------
## Everything above predicts NRS pain. This block re-runs the identical model
## set predicting SOWS (withdrawal symptoms) instead, and produces ONLY the
## regression-table figure -- one per drug group -- plus the numbers behind it.
## No diagnostics, added-variable plots, forests or comparison bars: the table
## is what was asked for, and the caveats that make those figures worth having
## for NRS apply here unchanged.
##
## SOWS cannot be a predictor of itself, so it is stripped out of the two
## adjusted models. setdiff() leaves the four receptor-only models untouched.

MREG_Y_SOWS      <- "SOWS"
MREG_MODELS_SOWS <- lapply(MREG_MODELS, setdiff, MREG_Y_SOWS)

## Labels have to lose " + SOWS" too, or the table header would advertise a
## predictor the model does not contain.
MREG_MODEL_LABELS_SOWS <- sub(" \\+ SOWS", "", MREG_MODEL_LABELS)
names(MREG_MODEL_LABELS_SOWS) <- names(MREG_MODEL_LABELS)

cat("\n================================================================\n")
cat("18l. MULTIPLE REGRESSION: ", MREG_Y_SOWS, " ~ log ROE + 5-HT4 / 5-HT6\n", sep = "")
cat("================================================================\n")

mreg_sows <- mreg_fit_all(MREG_MODELS_SOWS, MREG_Y_SOWS,
                          labels = MREG_MODEL_LABELS_SOWS)

cat("\n----- COEFFICIENTS (raw units: SOWS points per unit predictor) -----\n")
print(as.data.frame(mreg_sows$coef %>%
        dplyr::filter(term != "(Intercept)") %>%
        dplyr::select(group, model, term, est, se, ci_lo, ci_hi, t, p, beta_std, vif)),
      row.names = FALSE, digits = 3)

cat("\n----- MODEL FIT -----\n")
print(as.data.frame(mreg_sows$fit %>%
        dplyr::select(group, model, n, R2, adj_R2, F_stat, df1, df2, F_p,
                      RMSE, AIC, max_VIF, shapiro_p, bp_p)),
      row.names = FALSE, digits = 3)

for (g in names(mreg_sows$fits))
  mreg_table_figure(g, MREG_MODELS_SOWS, mreg_sows$coef, mreg_sows$fit,
                    MREG_Y_SOWS, dir = MREG_DIR_SOWS,
                    file = sprintf("MREG_SOWS_table_%s.png", mreg_slug(g)),
                    shorts = MREG_MODEL_SHORT_SOWS)

writexl::write_xlsx(
  list(
    Coefficients = mreg_sows$coef %>%
      dplyr::select(group, model, model_lab, term, term_lab, est, se, t, p, p_fmt,
                    ci_lo, ci_hi, beta_std, vif),
    Model_fit    = mreg_sows$fit,
    Config       = tibble::tibble(
      setting = c("outcome", "exposure", "groups", "models", "note", "output_dir"),
      value = c(
        MREG_Y_SOWS, MREG_X, paste(MREG_GROUPS, collapse = " | "),
        paste(sprintf("%s: %s ~ %s", names(MREG_MODELS_SOWS), MREG_Y_SOWS,
                      vapply(MREG_MODELS_SOWS, paste, character(1), collapse = " + ")),
              collapse = "  |  "),
        paste0("Same model set as the NRS analysis with SOWS moved from the right-hand ",
               "side to the left. Regression table only -- see the NRS output folder for ",
               "diagnostics, added-variable plots and the between-drug slope tests."),
        MREG_DIR_SOWS
      )
    )
  ),
  file.path(MREG_DIR_SOWS, "Stats_Multiple_Regression_HT4_HT6_SOWS.xlsx"))

message("Section 18l (SOWS) tables written to: ",
        normalizePath(MREG_DIR_SOWS, mustWork = FALSE))


## ---- 18m. NOTES --------------------------------------------------------------
## - Read the added-variable plots before the coefficient table. A coefficient
##   whose AV panel is a shapeless cloud with one far-out point is a coefficient
##   driven by that point, and the Cook's distance panel in the diagnostics
##   figure will name the subject.
## - n = 19 / 17 with 3-4 parameters in HT4_HT6. These are pre-specified
##   receptors, so no multiplicity correction is applied -- but the CIs are wide
##   and any result here needs replication before it is more than a lead.
##   Section 16's 19-receptor screen is where FDR belongs; this is not that.
## - The two ADJUSTED models are exploratory, full stop. HT4_HT6_Cov is 9
##   parameters on n ~ 18 (hydrocodone) and n ~ 13 (tramadol), leaving 9 and 4
##   residual df. Read the direction and the model-level fit; do not read an
##   individual covariate's CI as an estimate of anything. Sex in particular is
##   near-unidentified in the tramadol group, which is 20/23 female.
## - HT4_HT6_Cov and HT4_HT6_Cov_noMME are fit on DIFFERENT ROWS, because
##   log_MME is missing for 5 subjects (1 hydrocodone, 4 tramadol) and every
##   other predictor in this section is fully observed. Their R2, adj R2, AIC
##   and RMSE are therefore not comparable to each other or to the receptor-only
##   models -- that is exactly why the no-MME twin exists, so the comparison
##   against HT4_HT6 can be made on a common row set. The nested F test that IS
##   valid (HT4_HT6 -> HT4_HT6_Cov_noMME) is in the Nested_tests sheet; 18d's
##   guard would refuse the MME one, which is why it is not in MREG_NESTED.
## - Cross-check: HT4's R_5HT4 coefficient and HT6's R_5HT6 coefficient are the
##   same quantity as section 16's raw b_Receptor_raw for those receptors, and
##   the same as section 15's b path. Any disagreement means one of the three
##   sections is filtering rows differently.
## - The 4 subjects per group with ROE == 0 are dropped, not imputed, exactly as
##   in sections 15-17. If you want them in, change the log transform at section
##   2 to log10(ROE + 0.5*min(ROE[ROE > 0])) and re-source; do it as a labelled
##   sensitivity analysis, because the offset is arbitrary.
## - 18l repeats the regression table with SOWS as the outcome, into its own
##   output folder. It reuses mreg_fit_all() and mreg_table_figure() -- if you
##   change the table's layout, both outcomes move together, which is the point.
##   The SOWS run deliberately has no diagnostics or between-drug slope tests.
## - To re-run on a different outcome or different receptors, edit MREG_Y and
##   MREG_MODELS at 18a (and MREG_MODEL_LABELS, MREG_MODEL_SHORT and MREG_NESTED
##   to match) and re-source from 18a. Nothing downstream is hard-coded to two
##   receptors or to how many models there are -- except MREG_FULL just above
##   18g, which is
##   deliberately pinned to "HT4_HT6" rather than auto-following the list (see
##   its comment).
################################################################################
