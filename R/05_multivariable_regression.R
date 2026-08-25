################################################################################
## 16. MULTIVARIABLE REGRESSION: NRS ~ ROE + RECEPTOR
##     Fit separately for every receptor and every drug group.
##
## MODEL (per group, per receptor)
##   NRS = b0 + b_ROE * log_ROE + b_R * Receptor + e
##
## This is NOT a mediation model -- ROE and the receptor sit side by side as
## two predictors with no arrow between them. b_ROE is "the ROE-pain slope
## after netting out this receptor" and b_R is "the receptor-pain slope after
## netting out ROE". Compare to section 15: b_ROE here is the same quantity as
## c' there, and b_R here is the same quantity as b there -- section 15 already
## computed both for 5-HT4 and 5-HT6. This section is the general version,
## run once per receptor across the full 19-receptor panel, without the
## mediation machinery (no bootstrapped indirect effect, because there is no
## indirect effect to estimate in a model with no mediating arrow).
##
## WHY SEPARATELY BY GROUP
##   Same reasoning as section 15: the ROE-pain slope reverses sign between
##   drugs (hydrocodone +, tramadol -), so a pooled model would average the
##   two into something smaller than either. One model per group is the
##   stratified estimate; 16e below tests whether the receptor coefficients
##   actually differ between groups rather than just eyeballing two tables.
##
## WHY SEPARATELY BY RECEPTOR
##   Fitting all 19 receptors in one model per group would spend 19 df on
##   predictors with only 17-19 subjects to estimate them from -- the model
##   would be unidentified. One receptor at a time is the only version that
##   fits at this n. FDR correction (16d) is what keeps 19 separate tests
##   from inflating the false-positive rate.
##
## SCALING
##   As in section 15, all coefficients are standardized WITHIN GROUP (z-scored
##   on the group's own complete-case sample for that receptor), because raw
##   ROE differs ~40x in scale between drugs. Raw-unit coefficients are kept
##   alongside for reference.
##
## COLLINEARITY
##   ROE and receptor activity are typically only weakly correlated within
##   group (|r| mostly < 0.3 in this cohort), so the two predictors are not
##   competing for the same variance -- but r_predictors is reported per model
##   so you can see this rather than assume it.
##
## Self-contained: needs `master`, `FOCUS_GROUPS`, `RECEPTORS`, `OUT_DIR` and
## `fmt_p` (all from 00_config.R/01_load_data.R) and the packages loaded
## there (dplyr, purrr, tibble, ggplot2, writexl). No new packages required.
##
## Run standalone: Rscript R/05_multivariable_regression.R
## Sourced by: Run_All.R.
################################################################################

.HTK_PROJECT_DIR <- "/Users/kushagraverma/Work/Projects/Hydrocodone+Tramadol_Project"
## Always re-sourced -- see the note above this line in 01_load_data.R.
source(file.path(.HTK_PROJECT_DIR, "R", "00_config.R"))
if (!exists("master"))              source(file.path(.HTK_PROJECT_DIR, "R", "01_load_data.R"))

## ---- 16a. CONFIG -------------------------------------------------------------

MVR_X         <- "log_ROE"        # exposure predictor
MVR_RECEPTORS <- RECEPTORS        # full 19-receptor panel, defined in section 2
MVR_Y         <- "nrs"            # outcome; change to "PC2" etc. to re-run
MVR_GROUPS    <- FOCUS_GROUPS     # Hydrocodone group / Tramadol group
MVR_MIN_N     <- 12               # skip a group x receptor cell smaller than this

MVR_DIR <- file.path(OUT_DIR, "multivariable_regression")
dir.create(MVR_DIR, showWarnings = FALSE)

## fmt_p() comes from 00_config.R.

## Look up a coefficient by variable name, trying both the plain and
## backticked form -- R strips backticks from coefficient names when the
## variable name is syntactically valid, but not always.
mvr_pick <- function(nms, v) {
  hit <- match(v, nms)
  if (is.na(hit)) hit <- match(paste0("`", v, "`"), nms)
  if (is.na(hit)) stop("No coefficient found for variable: ", v,
                       "  (available: ", paste(nms, collapse = ", "), ")")
  hit
}


## ---- 16b. CORE MODEL FUNCTION -------------------------------------------------

## Fits NRS ~ ROE + Receptor on a single group x receptor cell.
## Returns both standardized (primary) and raw-unit (reference) results.
mvr_fit_one <- function(df, xv, rv, yv) {
  
  d <- df[, c(xv, rv, yv)]
  d <- d[stats::complete.cases(d), ]
  d <- d[apply(d, 1, function(r) all(is.finite(r))), ]
  n <- nrow(d)
  if (n < MVR_MIN_N) return(NULL)
  
  r_pred <- suppressWarnings(stats::cor(d[[xv]], d[[rv]]))
  
  ## raw-unit model, for reference
  m_raw <- lm(as.formula(sprintf("`%s` ~ `%s` + `%s`", yv, xv, rv)), data = d)
  c_raw <- summary(m_raw)$coefficients
  ix <- mvr_pick(rownames(c_raw), xv); ir <- mvr_pick(rownames(c_raw), rv)
  
  ## standardized model (within this cell's own complete-case sample)
  ds <- d
  for (v in c(xv, rv, yv)) {
    s <- stats::sd(ds[[v]])
    if (!is.finite(s) || s == 0) return(NULL)
    ds[[v]] <- (ds[[v]] - mean(ds[[v]])) / s
  }
  m_std <- lm(as.formula(sprintf("`%s` ~ `%s` + `%s`", yv, xv, rv)), data = ds)
  c_std <- summary(m_std)$coefficients
  jx <- mvr_pick(rownames(c_std), xv); jr <- mvr_pick(rownames(c_std), rv)
  
  f  <- summary(m_std)$fstatistic
  fp <- if (!is.null(f)) stats::pf(f[1], f[2], f[3], lower.tail = FALSE) else NA_real_
  vif <- if (is.finite(r_pred)) 1 / (1 - r_pred^2) else NA_real_   # 2-predictor VIF
  
  tibble(
    n = n, r_predictors = r_pred, vif = vif,
    b_ROE_std      = c_std[jx, 1], se_ROE_std      = c_std[jx, 2], p_ROE      = c_std[jx, 4],
    b_Receptor_std = c_std[jr, 1], se_Receptor_std = c_std[jr, 2], p_Receptor = c_std[jr, 4],
    b_ROE_raw      = c_raw[ix, 1], se_ROE_raw      = c_raw[ix, 2],
    b_Receptor_raw = c_raw[ir, 1], se_Receptor_raw = c_raw[ir, 2],
    R2 = summary(m_std)$r.squared, adj_R2 = summary(m_std)$adj.r.squared,
    F_p = fp
  )
}


## ---- 16c. FIT EVERY GROUP x RECEPTOR CELL -------------------------------------

mvr_rows <- list()

for (g in MVR_GROUPS) {
  d_g <- master[master$Plot_Group == g, , drop = FALSE]
  
  for (rv in MVR_RECEPTORS) {
    if (!rv %in% names(d_g)) next
    fit <- mvr_fit_one(d_g, MVR_X, rv, MVR_Y)
    if (is.null(fit)) next
    
    mvr_rows[[paste(g, rv)]] <- bind_cols(
      tibble(group = g, receptor = sub("^R_", "", rv), outcome = MVR_Y),
      fit
    )
  }
}

mvr_tbl <- bind_rows(mvr_rows)

## FDR correction: applied WITHIN each group, across the receptor panel --
## the receptor coefficient is the one being screened across 19 tests, so
## that is the family the correction is applied to. The ROE coefficient is
## the same exposure re-estimated 19 times with different covariates, not an
## independent family of hypotheses, but its p-value is FDR-corrected too for
## completeness -- expect it to barely move, since it does not depend much on
## which receptor happens to be in the model.
mvr_tbl <- mvr_tbl %>%
  group_by(group) %>%
  mutate(p_Receptor_FDR = p.adjust(p_Receptor, method = if (APPLY_FDR_CORRECTION) "fdr" else "none"),
         p_ROE_FDR      = p.adjust(p_ROE,      method = if (APPLY_FDR_CORRECTION) "fdr" else "none")) %>%
  ungroup() %>%
  mutate(sig_Receptor_raw = p_Receptor < 0.05,
         sig_Receptor_FDR = p_Receptor_FDR < 0.05) %>%
  arrange(group, p_Receptor)

cat("\n===== MULTIVARIABLE REGRESSION: NRS ~ log ROE + Receptor =====\n")
cat("(one model per receptor per group; FDR applied within each group across the 19-receptor panel)\n\n")
print(as.data.frame(mvr_tbl %>%
                      select(group, receptor, n, b_ROE_std, p_ROE, b_Receptor_std, p_Receptor,
                             p_Receptor_FDR, adj_R2)),
      row.names = FALSE, digits = 3)


## ---- 16d. HEADLINE HITS -------------------------------------------------------

cat("\n===== RECEPTORS SIGNIFICANT AT p < .05 (uncorrected) =====\n")
print(as.data.frame(mvr_tbl %>% filter(sig_Receptor_raw) %>%
                      select(group, receptor, n, b_Receptor_std, p_Receptor, p_Receptor_FDR, adj_R2)),
      row.names = FALSE, digits = 3)

cat("\n===== RECEPTORS SURVIVING FDR CORRECTION (q < .05) =====\n")
fdr_hits <- mvr_tbl %>% filter(sig_Receptor_FDR)
if (nrow(fdr_hits)) {
  print(as.data.frame(fdr_hits %>%
                        select(group, receptor, n, b_Receptor_std, p_Receptor, p_Receptor_FDR)),
        row.names = FALSE, digits = 3)
} else {
  cat("(none -- expected at n = 17-19 per group; treat raw hits above as leads,\n",
      " not confirmed effects, until replicated or pooled with a larger cohort)\n")
}


## ---- 16e. HYDROCODONE vs TRAMADOL, SIDE BY SIDE -------------------------------
## Same receptor, coefficient in each group next to each other, so a sign
## flip (a candidate "double dissociation", as with 5-HT4/5-HT6 in section 15)
## is visible without cross-referencing two separate tables.

mvr_wide <- mvr_tbl %>%
  select(group, receptor, n, b_Receptor_std, p_Receptor, p_Receptor_FDR) %>%
  pivot_wider(names_from = group,
              values_from = c(n, b_Receptor_std, p_Receptor, p_Receptor_FDR),
              names_sep = "_") %>%
  mutate(
    sign_flip = sign(.data[[paste0("b_Receptor_std_", MVR_GROUPS[1])]]) !=
      sign(.data[[paste0("b_Receptor_std_", MVR_GROUPS[2])]]),
    either_sig = .data[[paste0("p_Receptor_", MVR_GROUPS[1])]] < 0.05 |
      .data[[paste0("p_Receptor_", MVR_GROUPS[2])]] < 0.05
  ) %>%
  arrange(desc(sign_flip & either_sig))

cat("\n===== HYDROCODONE vs TRAMADOL: receptor coefficient, side by side =====\n")
cat("(sign_flip + either_sig flags candidate double dissociations, e.g. 5-HT4 / 5-HT6)\n\n")
print(as.data.frame(mvr_wide), row.names = FALSE, digits = 3)


## ---- 16f. FIGURES --------------------------------------------------------------

## Coefficient plot: every receptor's standardized effect on NRS (adjusted for
## ROE), one row per receptor, colored by group, ordered by hydrocodone effect
## size so the panel reads like Figures 2/5 but with inferential error bars.
if (nrow(mvr_tbl)) {
  
  ord <- mvr_tbl %>% filter(group == MVR_GROUPS[1]) %>%
    arrange(b_Receptor_std) %>% pull(receptor)
  mvr_tbl2 <- mvr_tbl %>%
    mutate(receptor = factor(receptor, levels = unique(c(ord, receptor))))
  
  p_coef <- ggplot(mvr_tbl2, aes(x = b_Receptor_std, y = receptor, colour = group)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmin = b_Receptor_std - 1.96 * se_Receptor_std,
                       xmax = b_Receptor_std + 1.96 * se_Receptor_std),
                   height = 0.2, position = position_dodge(width = 0.6)) +
    geom_point(aes(shape = sig_Receptor_raw), size = 2.6,
               position = position_dodge(width = 0.6)) +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                       labels = c(`TRUE` = "p < .05", `FALSE` = "n.s."),
                       name = NULL) +
    scale_colour_manual(values = if (exists("custom_colors")) custom_colors else NULL,
                        name = NULL) +
    labs(title = "Receptor -> NRS, adjusted for log ROE",
         subtitle = "Standardized coefficient +/- 95% CI, one multivariable model per receptor per group",
         x = "Standardized coefficient (SD units)", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
  
  ggsave(file.path(MVR_DIR, "MVR_receptor_coefficients_by_group.png"),
         p_coef, width = 9, height = 8.5, dpi = 300)
}

## Hydrocodone-vs-tramadol scatter: each receptor plotted at
## (hydrocodone coefficient, tramadol coefficient). Points in the upper-left
## or lower-right quadrant are dissociating -- opposite sign in the two
## groups, same as 5-HT4/5-HT6 sit relative to each other.
if (nrow(mvr_wide)) {
  xv <- paste0("b_Receptor_std_", MVR_GROUPS[1])
  yv <- paste0("b_Receptor_std_", MVR_GROUPS[2])
  
  ## Checked BEFORE the plot is built, not after: geom_text_repel() below is
  ## called at construction time, so a `requireNamespace()` guard sitting under
  ## the ggplot() call could never run -- the missing-package error had already
  ## been raised and taken the whole script down with it.
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    message("ggrepel not installed -- installing for label placement in the dissociation scatter")
    install.packages("ggrepel")
  }

  p_scatter <- ggplot(mvr_wide, aes(x = .data[[xv]], y = .data[[yv]])) +
    geom_hline(yintercept = 0, colour = "grey70") +
    geom_vline(xintercept = 0, colour = "grey70") +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey70") +
    geom_point(aes(colour = sign_flip & either_sig), size = 2.8) +
    ggrepel::geom_text_repel(aes(label = receptor), size = 3.2, max.overlaps = 20) +
    scale_colour_manual(values = c(`TRUE` = "#d95f02", `FALSE` = "grey40"),
                        guide = "none") +
    labs(title = "Receptor coefficient: hydrocodone vs tramadol",
         subtitle = "Orange = opposite sign in the two groups with at least one p < .05 (candidate dissociation)",
         x = paste0(MVR_GROUPS[1], " -- standardized coefficient"),
         y = paste0(MVR_GROUPS[2], " -- standardized coefficient")) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank())
  
  ggsave(file.path(MVR_DIR, "MVR_hydrocodone_vs_tramadol_scatter.png"),
         p_scatter, width = 7.5, height = 7, dpi = 300)
}


## ---- 16g. EXPORT ---------------------------------------------------------------

write_xlsx(
  list(
    Per_receptor_model = mvr_tbl,
    Side_by_side        = mvr_wide,
    FDR_hits             = fdr_hits
  ),
  file.path(MVR_DIR, "Stats_MVR_ROE_Receptor_NRS.xlsx")
)

message("\nMultivariable regression results and figures written to: ",
        normalizePath(MVR_DIR))


## ---- 16h. NOTES -----------------------------------------------------------------
## - At n = 17-19 per group, this is a 19-receptor screen, not a confirmatory
##   test. The FDR-corrected table (16d) is what you can defend without
##   qualification; the raw p < .05 table is a candidate list for receptors
##   already motivated by pharmacology (5-HT4, 5-HT6) or worth a replication
##   cohort, not a set of independent discoveries.
## - This model and section 15's b/c' paths estimate the same two quantities
##   for 5-HT4 and 5-HT6 (same formula, same data) -- the numbers in mvr_tbl
##   for those two receptors should match cprime_std/cp_se and b_std/b_se in
##   med_tbl exactly. Useful as a cross-check that both sections are correct.
## - To run this on a different outcome (e.g. PC2), change MVR_Y at 16a and
##   re-source from there; nothing else needs to change.
################################################################################
