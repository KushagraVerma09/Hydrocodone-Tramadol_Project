################################################################################
## 18. MULTIPLE REGRESSION: NRS ~ log ROE + 5-HT4 / 5-HT6
##     Five pre-specified models, fit independently in each drug group.
##     Section 18l repeats the regression table with SOWS as the outcome.
##
## MODELS (each fit separately within Hydrocodone group and Tramadol group).
## Model KEYS are on the left, the M-numbers they are referred to by in the
## figures in the middle; every model contains log_ROE.
##   HT4                M1  y = b0 + b_ROE*log_ROE + b_HT4*R_5HT4                + e
##   HT6                M2  y = b0 + b_ROE*log_ROE                + b_HT6*R_5HT6 + e
##   HT4_HT6            M3  y = b0 + b_ROE*log_ROE + b_HT4*R_5HT4 + b_HT6*R_5HT6 + e
##   HT4_HT6_ROEint     M4  M3 + b_4I*(log_ROE * R_5HT4) + b_6I*(log_ROE * R_5HT6)
##   HT4_HT6_ROEint_Cov M5  M4 + Sex + log10 MME + MQS non-opioid
##
## WHY THE INTERACTION MODEL (M4) EXISTS
##   M3 asks whether 5-HT4 and 5-HT6 carry SEPARATE information about the
##   outcome once exposure is netted out. It does not ask whether either
##   receptor's association DEPENDS ON EXPOSURE -- whether 5-HT4 tracks pain
##   more strongly in heavily exposed patients than in lightly exposed ones.
##   That second question is a statistical interaction (moderation) term:
##   literally the product of log ROE and the receptor value, one per receptor,
##   added to M3. Note that in M4 b_ROE, b_HT4 and b_HT6 become "this
##   predictor's effect when the variable it is multiplied by = 0" -- for
##   log_ROE that means ROE = 1, and for the receptors a value that never
##   occurs in this data -- so read the interaction coefficients and the nested
##   F test, not the three main effects.
##
## WHY THERE IS AN ADJUSTED MODEL (M5)
##   M5 asks whether whatever M4 finds survives adjustment for the obvious
##   confounders: sex, opioid dose (log10 MME) and non-opioid medication load
##   (MQS non-opioid). Sex and MQS non-opioid are fully observed in both drug
##   groups; MME is missing for 5 subjects (1 hydrocodone, 4 tramadol), so M5
##   is fit on FEWER ROWS than M1-M4 and its R2 / adj R2 / AIC / RMSE are not
##   comparable to theirs. Compare the n row before comparing anything else.
##   For the same reason there is no valid nested F test from M4 to M5 -- see
##   MREG_NESTED at 18a.
##
##   Read M5 as exploratory. It is 9 parameters on n ~ 18 (hydrocodone) and
##   n ~ 13 (tramadol), and the tramadol group is 20/23 female, so its sex
##   coefficient in particular is barely identified. Directions and the
##   model-level fit are worth reading; individual covariate CIs are too wide
##   to be more than a lead.
##
## HOW THIS DIFFERS FROM SECTION 16
##   Section 16 runs the SAME two-predictor model once per receptor across the
##   whole 19-receptor panel -- a screen. At n = 17-19 nothing survives FDR
##   there, and its two figures live entirely in coefficient space (19 rows of
##   overlapping CIs, and a scatter of sign flips) without ever showing a data
##   point. This section is the confirmatory counterpart: two receptors chosen
##   in advance on pharmacological grounds, reported the way a multiple
##   regression is normally reported -- coefficient table, added-variable
##   plots, observed vs predicted, residual diagnostics, model comparison.
##   The formal test of whether the slopes differ BETWEEN drugs is section 20
##   (R/09_pooled_regression.R), which pools the groups instead.
##   HT4 and HT6 are numerically identical to section 16's R_5HT4 and R_5HT6
##   rows (same formula, same data); that is a deliberate cross-check.
##
## WHY M3 IS THE INTERESTING ONE
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
##   either. One model per group is the stratified estimate; section 20
##   (R/09_pooled_regression.R) then tests whether the slopes actually differ
##   rather than leaving it to the eye. Read that test first: if the drugs do
##   NOT differ, splitting them here is spending degrees of freedom for
##   nothing, and the pooled model is the better summary.
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
##         interaction test in section 20 is still a valid test of raw-unit slope
##         equality; just do not read the two raw numbers as effect sizes.
##   Covariate scales in M5: MME enters as log10 (it is heavily right-skewed
##   and no subject has MME == 0, so nothing is lost), MQS non-opioid enters
##   raw, and sex enters as a 0/1 indicator, sex_female. The two interaction
##   terms in M4/M5 are the products of log10 ROE with the RAW receptor value,
##   so their coefficients read as "change in the receptor's slope per one-unit
##   increase in log10 ROE".
##
## THE -Inf TRAP
##   log_ROE = log10(ROE) is -Inf for the 4 subjects per group with ROE == 0
##   when ROE_ZERO_IMPUTE (00_config.R) is FALSE, and complete.cases() does
##   NOT drop -Inf. Every model frame below is filtered with is.finite(),
##   matching sections 16 and 17. Expect n = 19 (hydrocodone) and n = 17
##   (tramadol) with the switch off; with it on (the default -- see
##   01_load_data.R) those subjects get a finite, imputed log_ROE instead and
##   n rises to 23 / 21.
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
## fits, 18l's SOWS models and section 20's pooled fits all see the identical
## columns.
##
## (a) The two interaction terms: exposure x receptor, one per receptor. Built
## on log_ROE (not raw ROE), because log_ROE is what enters the models as a
## main effect and an interaction has to be the product of the terms actually
## in the model. log_ROE is -Inf for the ROE == 0 subjects, so these products
## are non-finite for exactly those rows -- the is.finite() filter in
## mreg_fit_one() drops them, which is the same set of rows log_ROE itself
## would have cost.
master$R_5HT4_x_ROE <- master$R_5HT4 * master$log_ROE
master$R_5HT6_x_ROE <- master$R_5HT6 * master$log_ROE

## (b) Sex as a numeric 0/1 indicator. `master$female` (01_load_data.R:100) is
## LOGICAL, and handing lm() a logical makes it name the coefficient
## "femaleTRUE" rather than "female" -- which silently breaks every place
## downstream that matches a coefficient name against a predictor name (the VIF
## lookup in 18b, the added-variable plots in 18g, the partial-residual
## contributions in 18i, and the regression table's row ordering). A numeric
## column keeps coefficient name == column name. 1 = female, 0 = male.
master$sex_female <- as.numeric(master$female)

## The five models, M1 - M5 in the order they appear everywhere downstream.
## Each entry is the full predictor vector, in the order the terms should
## appear. Change this list (and MREG_Y) and re-source from 18a to re-run the
## whole section on different predictors -- nothing downstream is hard-coded to
## 5-HT4 / 5-HT6 or to five models.
##
## M4 adds BOTH exposure x receptor products at once rather than one at a time:
## they are asked as one question ("does exposure moderate the receptor-outcome
## association at all?"), and the nested F test in 18d tests them jointly on 2
## df, which is the right test for that question.
##
## M5 = M4 + sex + log10 MME + MQS non-opioid. MME is the only predictor in
## this section with missing values (5 subjects), so M5 is fit on fewer rows
## than M1-M4 -- see the header, 18d's guard, and the n row of every table.
MREG_MODELS <- list(
  HT4                = c("log_ROE", "R_5HT4"),
  HT6                = c("log_ROE", "R_5HT6"),
  HT4_HT6            = c("log_ROE", "R_5HT4", "R_5HT6"),
  HT4_HT6_ROEint     = c("log_ROE", "R_5HT4", "R_5HT6",
                         "R_5HT4_x_ROE", "R_5HT6_x_ROE"),
  HT4_HT6_ROEint_Cov = c("log_ROE", "R_5HT4", "R_5HT6",
                         "R_5HT4_x_ROE", "R_5HT6_x_ROE",
                         "sex_female", "log_MME", "MQS_nonopioid")
)

## Full model labels -- each one names its own predictors, so a reader looking
## at a facet strip or a title never has to come back up here to find out what
## the model contains. No M-number here: these are used one at a time in
## titles and subtitles, where the formula is the whole point. The M-numbers
## live in MREG_MODEL_SHORT below, which is what labels the columns and bars
## where several models sit side by side.
MREG_MODEL_LABELS <- c(
  HT4                = "ROE + 5-HT4",
  HT6                = "ROE + 5-HT6",
  HT4_HT6            = "ROE + 5-HT4 + 5-HT6",
  HT4_HT6_ROEint     = "ROE + 5-HT4 + 5-HT6 + ROEx5-HT4 + ROEx5-HT6",
  HT4_HT6_ROEint_Cov = paste("ROE + 5-HT4 + 5-HT6 + ROEx5-HT4 + ROEx5-HT6",
                             "+ Sex + log MME + MQS non-opioid")
)

## Nested comparisons to run in 18d: c(reduced, full).
##
## M4 -> M5 is deliberately ABSENT. log_MME has missing values, so M5 is fit on
## fewer rows than M4, and a nested F test across different row sets is not
## valid -- 18d's own guard would refuse it anyway. There is no no-MME twin of
## M5 to compare instead: the model set was specified with MME in it, and
## adding an unrequested sixth model to make one comparison work would put a
## column in the figures that nobody asked for. The cost of MME is visible in
## the n row of every table instead.
MREG_NESTED <- list(c("HT4", "HT4_HT6"), c("HT6", "HT4_HT6"),
                    c("HT4_HT6", "HT4_HT6_ROEint"))

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
  R_5HT4_x_ROE  = "5-HT4 x log10 ROE",
  R_5HT6_x_ROE  = "5-HT6 x log10 ROE",
  sex_female    = "Sex (female = 1)",
  MQS_nonopioid = "MQS non-opioid",
  log_MME       = "log10 MME",
  ## Not a predictor in any model here any more -- kept because 18l uses SOWS
  ## as the OUTCOME, and mreg_lab() is what titles that run's figures.
  SOWS          = "SOWS withdrawal"
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
  out[out == "5HT4_x_ROE"]    <- "5-HT4 x log10 ROE"
  out[out == "5HT6_x_ROE"]    <- "5-HT6 x log10 ROE"
  out[out == "log_ROE"]       <- "log10 ROE"
  out[out == "log_MME"]       <- "log10 MME"
  out[out == "sex_female"]    <- "Sex (female)"
  out[out == "MQS_nonopioid"] <- "MQS non-opioid"
  out[out == "DOU"]           <- "DOU (yr)"   # no longer a predictor; harmless
  out
}

## Short MODEL labels, for the spots with only enough room for one model per
## column/bar rather than one per row (the regression-table header (4), the
## observed-vs-predicted facet strips (5) and the model-comparison x-axis (6)).
## MREG_MODEL_LABELS's full formula reads fine one at a time in a
## title/subtitle, but five of them side by side overlap into unreadable text.
##
## Each one leads with its M-number, because that is how the models are
## referred to when they are discussed as a set, and then says what the model
## adds -- every model here contains ROE, so M4/M5 only need to name the new
## terms. The full predictor list is still visible in the table's own row
## labels and in every facet strip. Keyed the same way as MREG_MODEL_LABELS and
## pre-wrapped by hand.
MREG_MODEL_SHORT <- c(
  HT4                = "M1\nROE +\n5-HT4",
  HT6                = "M2\nROE +\n5-HT6",
  HT4_HT6            = "M3\nROE + 5-HT4\n+ 5-HT6",
  HT4_HT6_ROEint     = "M4\n+ ROE x 5-HT4\n+ ROE x 5-HT6",
  HT4_HT6_ROEint_Cov = "M5\n+ Sex, log MME,\nMQS non-opioid"
)

## 18l used to need its own short labels, because the adjusted models then
## contained SOWS and SOWS is that run's OUTCOME. No model contains SOWS any
## more, so the two label sets are identical -- the alias is kept because
## 02_table1.R passes MREG_MODEL_SHORT_SOWS by name when it builds Table 3, and
## because a future model set could reintroduce the clash.
MREG_MODEL_SHORT_SOWS <- MREG_MODEL_SHORT

## Look a short label up BY MODEL KEY, falling back to the key itself so a
## model added to MREG_MODELS without a short-label entry still draws.
mreg_model_short <- function(keys, tbl = MREG_MODEL_SHORT) {
  out <- unname(tbl[keys])
  ifelse(is.na(out), keys, out)
}

## Just the M-number ("M3"), for the spots that refer to a model in running
## text rather than labelling a column -- the nested-test lines under the
## model-comparison figure. Read off the FIRST LINE of the short label so the
## numbering has exactly one source: renumber MREG_MODEL_SHORT and every
## mention follows. Falls back to the key for a model with no short label.
mreg_model_num <- function(keys, tbl = MREG_MODEL_SHORT) {
  out <- sub("\n.*$", "", unname(tbl[keys]))
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
    ## has two pairs that share the same FULL model (HT4 -> HT4_HT6 and
    ## HT6 -> HT4_HT6), and has held pairs sharing a reduced model too, so a key
    ## of just one half of the pair would make one overwrite the other here.
    mreg_nested_rows[[paste(g, pair[1], pair[2])]] <- tibble::tibble(
      group        = g,
      reduced      = pair[1],
      full         = pair[2],
      added_term   = paste(added, collapse = " + "),
      ## Pretty version of the same thing, for figure (6)'s subtitle. Has to be
      ## built by mapping mreg_short() over the terms BEFORE they are pasted:
      ## run on the joined string it would only rewrite the first one.
      added_lab    = paste(mreg_short(added), collapse = " + "),
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


## ---- 18e. MOVED ---------------------------------------------------------------
## The pooled group-difference models that used to live here are now section 20
## (R/09_pooled_regression.R), along with the three figures that were built from
## them (the per-receptor marginal-effect plots, the interaction-slopes figure
## and the slope-difference forest) and their two .xlsx sheets.
##
## The split is by question, not by convenience: everything left in this file
## models ONE drug at a time, and everything in 09 compares the two. 09 runs
## after this file and reuses its helpers and its stratified fits (`mreg_fits`),
## so nothing here needs to know about it.


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
## `caption` is overridable because section 20 reuses this function for the
## pooled models, where the default text is wrong twice over: those models have
## no standardized betas to point at, and their columns do not "differ in n".
## Default unchanged, so sections 18h and 18l are untouched.
mreg_table_figure <- function(g, models, coef_tbl, fit_tbl, yv,
                              dir = MREG_DIR, file = NULL,
                              shorts = MREG_MODEL_SHORT,
                              caption = paste0(
                                "*** p < .001, ** p < .01, * p < .05.  ",
                                "Raw coefficients are not comparable across predictors (different scales); ",
                                "see the .xlsx for standardized betas.\n",
                                "VIF > 5 would indicate predictors competing for the same variance.  ",
                                "Columns differ in n where a predictor has missing values -- compare the n row.")) {
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
  cellsdf$x <- cellsdf$col

  ## Header cells are the SHORT model labels, looked up by model key.
  hdr <- data.frame(x = seq_len(ncol_tot), y = 0,
                    txt = c("", mreg_model_short(mods, shorts)),
                    stringsAsFactors = FALSE)

  ## HEADER GEOMETRY. The header text is CENTERED on y = 0 and is as many lines
  ## tall as the longest short label -- three, now that those labels lead with
  ## an M-number, where they used to be two. The rules were hardcoded at +0.55
  ## and -0.45, sized for exactly two lines, so a third line struck through the
  ## lower rule. Everything below is derived from hdr_lines instead.
  ##
  ## 0.20 per line is calibrated to reproduce the old hardcoded numbers at two
  ## lines (2 * 0.20 + 0.05 = 0.45), so a two-line header still lays out
  ## exactly as it did before this became dynamic.
  hdr_lines <- max(vapply(strsplit(hdr$txt, "\n", fixed = TRUE), length, integer(1)))
  hdr_half  <- 0.20 * hdr_lines
  rule_hdr  <- -(hdr_half + 0.05)

  ## Body rows sit 1 y-unit apart and hold two lines, so row 1's text reaches
  ## up to about -1 + 0.45. Once the header grows past two lines that collides
  ## with rule_hdr, so shift the whole body down by however much it overlaps
  ## (plus 0.10 of clearance). At two lines this is 0 and nothing moves.
  body_off <- max(0, (-rule_hdr + 0.45 + 0.10) - 1)

  ## y descends down the page; header sits above row 1.
  cellsdf$y <- -(cellsdf$row + body_off)

  p_tab <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = c(hdr_half + 0.15, rule_hdr,
                                       -(n_body + 0.45 + body_off),
                                       -(length(all_rows) + 0.55 + body_off)),
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
    ggplot2::scale_y_continuous(limits = c(-(length(all_rows) + 1.1 + body_off),
                                           hdr_half + 1.2),
                                expand = c(0, 0)) +
    ggplot2::labs(
      title    = sprintf("%s regressed on log10 ROE and receptor activity -- %s",
                         mreg_lab(yv), g),
      subtitle = "Raw-unit coefficient, standard error in parentheses",
      caption  = caption
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

  ## 2.3 rather than the original 2.6 per column: with five models the table is
  ## 6 columns wide, and the old factor made it needlessly wide. Height carries
  ## the extra header lines too, so the taller header does not eat into the
  ## body's share of a fixed-height device.
  mreg_save(p_tab, file, dir = dir,
            w = 2.3 * ncol_tot + 1.4,
            h = 0.38 * length(all_rows) + 2.4 + 0.38 * max(0, hdr_lines - 2))
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
  ## 3.2 rather than 3.5 per column now that there are six of them.
  mreg_save(p_op, "MREG_obs_vs_pred.png",
            w = 3.2 * length(unique(op_df$model_lab)) + 2, h = 8.6)
}

## (6) Model comparison: adjusted R2 and AIC side by side, with the nested F
## tests annotated. adj R2 and AIC both penalise the extra parameter, so if the
## two-receptor model does not beat the single-receptor ones on these it is not
## earning its keep.
##
## CAUTION: M5 (HT4_HT6_ROEint_Cov) is fit on FEWER rows than every other bar
## here (log_MME is missing for 5 subjects), and R2/AIC computed on different
## row sets are not comparable. Its bar is legitimate as a description of that
## model; it is not a like-for-like contest against M1-M4, which are all fit on
## the same rows and can be read against each other freely.
if (nrow(mreg_fit_tbl)) {
  comp <- mreg_fit_tbl %>%
    dplyr::select(group, model, n, R2, adj_R2, AIC) %>%
    tidyr::pivot_longer(c(R2, adj_R2, AIC), names_to = "metric", values_to = "value") %>%
    dplyr::mutate(
      metric = factor(metric, levels = c("R2", "adj_R2", "AIC"),
                      labels = c("R2 (in-sample)", "Adjusted R2", "AIC (lower is better)")),
      ## x is the SHORT label, keyed off the model name -- six full formulas
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
    per_row <- sprintf("%s:  %s -> %s, adding %s:  F(%d,%d) = %.2f, p = %s",
                       mreg_nested_tbl$group,
                       mreg_model_num(mreg_nested_tbl$reduced),
                       mreg_model_num(mreg_nested_tbl$full),
                       mreg_nested_tbl$added_lab, mreg_nested_tbl$df_num,
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
  ## per metric panel); 2.6 per model rather than 3.5, since five models at the
  ## old factor made a 20-inch-wide figure. Height scales with how many
  ## nested-test lines are in the subtitle.
  mreg_save(p_comp, "MREG_model_comparison.png",
            w = 2.6 * length(MREG_MODELS) + 3,
            h = 5.5 + 0.3 * (nrow(mreg_nested_tbl) + 1))
}

## (7), (8), (9) MOVED
## The per-receptor marginal-effect plots, the drug-interaction-slopes figure
## and the slope-difference forest are now section 20
## (R/09_pooled_regression.R), with the pooled models they annotate. All three
## are pictures of how the two drugs COMPARE; everything still in this file is
## a picture of one drug at a time.
##
## Figure (7) is drawn from the stratified fits in `mreg_fits`, which stay
## here -- 09 runs afterwards and reads them.


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
    paste0("R_5HT4_x_ROE = R_5HT4 * log_ROE and R_5HT6_x_ROE = R_5HT6 * log_ROE ",
           "(M4, M5)"),
    "sex_female = 1 if female (sex == 2), 0 if male; from master$female",
    paste0("MME entered as log10 (log_MME); MQS_nonopioid entered raw. M5 ",
           "(HT4_HT6_ROEint_Cov) is fit on fewer rows than M1-M4 because log_MME ",
           "has missing values -- its R2/adj R2/AIC/RMSE are NOT comparable to ",
           "theirs, and no nested F test into it is valid. Compare the n row first."),
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
  ## Interaction_terms and Interaction_omnibus moved to section 20's workbook
  ## (New_Outputs/between_drug_comparison/Stats_Between_Drug_Comparison.xlsx)
  ## along with the models that produce them. Same contents, different file.
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
## SOWS cannot be a predictor of itself, so it is stripped out of every model.
## As of the current M1-M5 set no model contains SOWS, so both the setdiff()
## and the label surgery below are no-ops -- they are kept as guards, so that
## putting SOWS back into a model at 18a cannot silently produce a SOWS ~ SOWS
## regression here.

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
## - n = 19 / 17 with 3-4 parameters in M3. These are pre-specified receptors,
##   so no multiplicity correction is applied -- but the CIs are wide and any
##   result here needs replication before it is more than a lead.
##   Section 16's 19-receptor screen is where FDR belongs; this is not that.
## - M4 spends 2 df on the two exposure x receptor products at n = 17-19. Read
##   the joint nested F test (M3 -> M4) before either interaction coefficient:
##   with this much data one product term can look large on its own noise.
##   Remember also that adding them changes what b_HT4 and b_HT6 MEAN -- they
##   are now the receptor slopes at log10 ROE = 0, i.e. at ROE = 1 -- so a main
##   effect that moves between M3 and M4 has not necessarily changed, it is
##   being reported at a different point.
## - M5 is exploratory, full stop. It is 9 parameters on n ~ 18 (hydrocodone)
##   and n ~ 13 (tramadol), leaving 9 and 4 residual df. Read the direction and
##   the model-level fit; do not read an individual covariate's CI as an
##   estimate of anything. Sex in particular is near-unidentified in the
##   tramadol group, which is 20/23 female.
## - M5 is fit on DIFFERENT ROWS from M1-M4, because log_MME is missing for 5
##   subjects (1 hydrocodone, 4 tramadol) and every other predictor in this
##   section is fully observed. Its R2, adj R2, AIC and RMSE are therefore not
##   comparable to theirs, and no nested F test into it is valid -- which is
##   why M4 -> M5 is not in MREG_NESTED (18d's guard would refuse it anyway).
##   The n row of every table is where that cost is visible.
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
##   receptors or five models -- except MREG_FULL at 20a in
##   R/09_pooled_regression.R, which is deliberately pinned to "HT4_HT6"
##   rather than auto-following the list (see its comment), and REGTAB_MODELS
##   in 02_table1.R, which names the three
##   models (M3, M4, M5) that Tables 2 and 3 publish. If you rename a model key
##   here, rename it there too.
## - Any interaction term you add at 18a has to exist as a COLUMN on `master`
##   (built at 18a, above MREG_MODELS), not as an I(a*b) or a:b in the formula.
##   The VIF loop, the added-variable plots, the partial-residual contributions
##   and the table's row ordering all match coefficient names against predictor
##   names, and R's formula notation produces coefficient names those lookups
##   would miss.
################################################################################
