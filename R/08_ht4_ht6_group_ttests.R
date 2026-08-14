################################################################################
## 19. GROUP DIFFERENCE IN 5-HT4 / 5-HT6 ACTIVITY
##     Hydrocodone group vs Tramadol group
##
## THE QUESTION
##   Every other section in this pipeline uses 5-HT4 and 5-HT6 as PREDICTORS of
##   pain (sections 15-18). This section asks the prior, simpler question:
##
##       do the two drug groups differ in the receptor values themselves?
##
##   i.e. for each measure below, test  H0: mean(Hydrocodone) = mean(Tramadol).
##
## WHAT IS TESTED (three measures, same test each time)
##   5-HT4                R_5HT4
##   5-HT6                R_5HT6
##   5-HT4 & 5-HT6 mean   rowMeans(R_5HT4, R_5HT6)   -- the "HT6_plus_HT4"
##                        receptor set, i.e. the x-axis value those figures use
##
##   The receptor COLUMNS are read out of `receptor_sets` (01_load_data.R)
##   rather than hardcoded a second time, same as section 11 does for HT6, so
##   renaming or remapping a set updates this section automatically instead of
##   silently going stale.
##
## THE MODEL
##   Welch two-sample t-test (primary):
##       value = b0 + b_group * 1{group = Tramadol} + e ,  e ~ N(0, sigma_g^2)
##       t = (mean_H - mean_T) / sqrt(sd_H^2/n_H + sd_T^2/n_T)
##   Welch is the primary test because it does NOT assume the two groups share
##   a variance, and the group sizes here are unequal -- the case where the
##   equal-variance (Student) test misstates its own p-value. Student's t and
##   the Wilcoxon rank-sum test are reported alongside as robustness checks,
##   not as alternatives to choose between after seeing the p-values.
##
## WHY THIS SECTION EXISTS -- THE OVERLAPPING-CI-BARS QUESTION
##   Reading significance off two error bars is not the same test as the
##   t-test, and the two can disagree in BOTH directions:
##
##     - each group's bar is  mean +/- t* * sd/sqrt(n)   (uncertainty in ONE mean)
##     - the t-test is about  (mean_H - mean_T) +/- t* * sqrt(se_H^2 + se_T^2)
##                                                     (uncertainty in the GAP)
##
##   Because sqrt(se_H^2 + se_T^2) < se_H + se_T, two bars can OVERLAP while
##   p < .05. The reverse also happens: two bars can just barely clear each
##   other while p > .05. So this section reports, for every measure:
##     (a) the per-group means with their 95% CIs   -- what the bars show
##     (b) the difference with its 95% CI and p     -- what the test says
##     (c) whether the two per-group CIs overlap    -- so (a) vs (b) disagreeing
##                                                     is visible, not confusing
##   Figure 2 below plots (b) directly, which is the plot that actually answers
##   "are these two groups different?".
##
## WHY A NON-SIGNIFICANT RESULT HERE IS NOT "NO DIFFERENCE"
##   n is ~23-25 per group here -- larger than the regression sections, which
##   are complete-case on ROE/NRS and lose subjects to missing predictors,
##   because this test needs only the receptor value. 19c reports the smallest
##   difference this design
##   could detect with 80% power at alpha = .05. Any true difference smaller
##   than that is one this sample cannot rule in or out -- so a large p here
##   means "not resolved at this n", not "equal".
##
## Self-contained: needs `master`, `FOCUS_GROUPS`, `receptor_sets`,
## `compute_x_value`, `custom_colors`, `OUT_DIR` and `fmt_p` (all from
## 00_config.R / 01_load_data.R) and the packages loaded there (dplyr, purrr,
## tibble, ggplot2, writexl). No new packages required.
##
## Run standalone: Rscript R/08_ht4_ht6_group_ttests.R
## Sourced by: Run_All.R.
################################################################################

.HTK_PROJECT_DIR <- "/Users/kushagraverma/Work/Projects/Hydrocodone+Tramadol_Project"
## Always re-sourced -- see the note above this line in 01_load_data.R.
source(file.path(.HTK_PROJECT_DIR, "R", "00_config.R"))
if (!exists("master"))              source(file.path(.HTK_PROJECT_DIR, "R", "01_load_data.R"))

## ---- 19a. CONFIG ---------------------------------------------------------------

HT_DIR <- file.path(OUT_DIR, "ht4_ht6_group_comparison")
dir.create(HT_DIR, showWarnings = FALSE)

## Which measures to test, as receptor-set names from 01_load_data.R.
## Add a set name here and it gets the full treatment below -- table, both
## figures, FDR correction -- with nothing else to edit.
HT_SETS <- c("5-HT4"              = "HT4_Only",
             "5-HT6"              = "HT6_Only",
             "5-HT4 & 5-HT6 mean" = "HT6_plus_HT4")

stopifnot(all(HT_SETS %in% names(receptor_sets)))

HT_CONF <- 0.95   # confidence level for every CI in this section

## ---- 19b. BUILD THE ANALYSIS FRAME ---------------------------------------------
## One column per measure, computed the same way the figures compute their
## x-axis (compute_x_value), restricted to the two focus groups.

ht_data <- master %>%
  filter(Plot_Group %in% FOCUS_GROUPS) %>%
  mutate(Plot_Group = droplevels(factor(Plot_Group, levels = FOCUS_GROUPS)))

for (lab in names(HT_SETS))
  ht_data[[lab]] <- compute_x_value(ht_data, receptor_sets[[HT_SETS[lab]]])

message(sprintf("Comparing %s: %s (n=%d) vs %s (n=%d)",
                paste(names(HT_SETS), collapse = ", "),
                FOCUS_GROUPS[1], sum(ht_data$Plot_Group == FOCUS_GROUPS[1]),
                FOCUS_GROUPS[2], sum(ht_data$Plot_Group == FOCUS_GROUPS[2])))

## ---- 19c. THE TEST -------------------------------------------------------------
## One row per measure. Everything a reader needs to judge the comparison --
## the two group means with their own CIs, the difference with ITS CI, three
## p-values, an effect size with a CI, the assumption checks, and the smallest
## difference this n could have detected -- lives in that one row.

ht_compare_one <- function(d, measure, gvar = "Plot_Group",
                           g1 = FOCUS_GROUPS[1], g2 = FOCUS_GROUPS[2],
                           conf = HT_CONF) {

  x1 <- d[[measure]][d[[gvar]] == g1]; x1 <- x1[!is.na(x1)]
  x2 <- d[[measure]][d[[gvar]] == g2]; x2 <- x2[!is.na(x2)]
  n1 <- length(x1); n2 <- length(x2)

  base <- tibble(measure = measure, group_1 = g1, group_2 = g2,
                 n_1 = n1, n_2 = n2)

  if (n1 < 3 || n2 < 3) {
    message("  ", measure, ": fewer than 3 usable values in a group -- skipped.")
    return(base %>% mutate(p_welch = NA_real_))
  }

  se1 <- sd(x1) / sqrt(n1)
  se2 <- sd(x2) / sqrt(n2)
  ## Per-group mean CIs -- these are the bars on Figure 1.
  ci1 <- mean(x1) + c(-1, 1) * qt(1 - (1 - conf) / 2, n1 - 1) * se1
  ci2 <- mean(x2) + c(-1, 1) * qt(1 - (1 - conf) / 2, n2 - 1) * se2

  welch   <- stats::t.test(x1, x2, var.equal = FALSE, conf.level = conf)
  student <- stats::t.test(x1, x2, var.equal = TRUE,  conf.level = conf)
  ## Ties are near-certain in continuous receptor data only when values repeat,
  ## but exact = FALSE keeps the normal approximation (and no warning) either way.
  wilcox  <- suppressWarnings(stats::wilcox.test(x1, x2, exact = FALSE))

  ## Effect size on the pooled SD, with the small-sample (Hedges) correction
  ## and a large-sample CI. At n ~ 18 per group the uncorrected d is biased
  ## upward by ~2%, so both are reported rather than picking one.
  s_pooled <- sqrt(((n1 - 1) * var(x1) + (n2 - 1) * var(x2)) / (n1 + n2 - 2))
  d_eff    <- (mean(x1) - mean(x2)) / s_pooled
  J        <- 1 - 3 / (4 * (n1 + n2 - 2) - 1)
  se_d     <- sqrt((n1 + n2) / (n1 * n2) + d_eff^2 / (2 * (n1 + n2)))
  d_ci     <- d_eff + c(-1, 1) * qnorm(1 - (1 - conf) / 2) * se_d

  ## Smallest difference detectable at 80% power, alpha = .05, at this n.
  ## Read next to the observed difference: if the observed gap is well below
  ## this, a large p-value is a statement about the sample size, not the brains.
  mde <- tryCatch(
    stats::power.t.test(n = 2 / (1 / n1 + 1 / n2),   # harmonic mean n per group
                        sd = s_pooled, sig.level = 1 - conf,
                        power = 0.80)$delta,
    error = function(e) NA_real_)

  base %>% mutate(
    mean_1 = mean(x1), sd_1 = sd(x1), se_1 = se1,
    mean_1_ci_lo = ci1[1], mean_1_ci_hi = ci1[2],
    mean_2 = mean(x2), sd_2 = sd(x2), se_2 = se2,
    mean_2_ci_lo = ci2[1], mean_2_ci_hi = ci2[2],
    ## Difference and its CI: group 1 minus group 2, the quantity being tested.
    mean_diff  = unname(diff(rev(welch$estimate))),
    diff_ci_lo = welch$conf.int[1],
    diff_ci_hi = welch$conf.int[2],
    t_welch    = unname(welch$statistic),
    df_welch   = unname(welch$parameter),
    p_welch    = welch$p.value,
    t_student  = unname(student$statistic),
    df_student = unname(student$parameter),
    p_student  = student$p.value,
    W_wilcoxon = unname(wilcox$statistic),
    p_wilcoxon = wilcox$p.value,
    cohens_d   = d_eff,
    hedges_g   = d_eff * J,
    d_ci_lo    = d_ci[1], d_ci_hi = d_ci[2],
    ## Assumption checks for the two t-tests, reported rather than acted on:
    ## Welch is already the answer to a failed variance test, and at this n a
    ## Shapiro p > .05 is weak evidence of normality, not a clean pass.
    p_shapiro_1 = stats::shapiro.test(x1)$p.value,
    p_shapiro_2 = stats::shapiro.test(x2)$p.value,
    p_var_equal = stats::var.test(x1, x2)$p.value,
    ## Do the two per-group bars on Figure 1 overlap? Compared against
    ## significant_welch below, this is what shows that eyeballing the bars
    ## and running the test are two different questions.
    ci_bars_overlap = (ci1[1] <= ci2[2]) && (ci2[1] <= ci1[2]),
    significant_welch = welch$p.value < (1 - conf),
    mde_80_power = mde)
}

ht_tbl <- map_dfr(names(HT_SETS), ~ ht_compare_one(ht_data, .x)) %>%
  ## Three pre-specified measures, and the third is a function of the first
  ## two -- so BH across them is mild and reported for completeness. The two
  ## receptor rows are the pre-specified comparisons; read the raw p there.
  mutate(p_welch_FDR = p.adjust(p_welch, "fdr"), .after = p_welch)

## ---- 19d. CONSOLE REPORT -------------------------------------------------------

cat("\n", strrep("=", 79), "\n", sep = "")
cat("5-HT4 / 5-HT6 ACTIVITY: ", FOCUS_GROUPS[1], " vs ", FOCUS_GROUPS[2], "\n", sep = "")
cat(strrep("=", 79), "\n", sep = "")

for (i in seq_len(nrow(ht_tbl))) {
  r <- ht_tbl[i, ]
  if (is.na(r$p_welch)) next
  cat(sprintf("\n%s\n", r$measure))
  cat(sprintf("  %-20s mean = %+.4f  (95%% CI %+.4f to %+.4f)  SD = %.4f  n = %d\n",
              r$group_1, r$mean_1, r$mean_1_ci_lo, r$mean_1_ci_hi, r$sd_1, r$n_1))
  cat(sprintf("  %-20s mean = %+.4f  (95%% CI %+.4f to %+.4f)  SD = %.4f  n = %d\n",
              r$group_2, r$mean_2, r$mean_2_ci_lo, r$mean_2_ci_hi, r$sd_2, r$n_2))
  cat(sprintf("  difference (%s - %s) = %+.4f  (95%% CI %+.4f to %+.4f)\n",
              sub(" group", "", r$group_1), sub(" group", "", r$group_2),
              r$mean_diff, r$diff_ci_lo, r$diff_ci_hi))
  cat(sprintf("  Welch t(%.1f) = %.2f, p = %s%s   Hedges g = %.2f (95%% CI %.2f to %.2f)\n",
              r$df_welch, r$t_welch, fmt_p(r$p_welch),
              if (r$p_welch < 0.001) " ***" else
                if (r$p_welch < 0.01) " **" else
                  if (r$p_welch < 0.05) " *" else "  (n.s.)",
              r$hedges_g, r$d_ci_lo, r$d_ci_hi))
  cat(sprintf("  robustness: Student t p = %s | Wilcoxon p = %s | BH-adjusted Welch p = %s\n",
              fmt_p(r$p_student), fmt_p(r$p_wilcoxon), fmt_p(r$p_welch_FDR)))
  cat(sprintf("  assumptions: Shapiro p = %s / %s, equal-variance F test p = %s\n",
              fmt_p(r$p_shapiro_1), fmt_p(r$p_shapiro_2), fmt_p(r$p_var_equal)))
  cat(sprintf("  the two 95%% CI bars %s.  The test says %s.\n",
              if (r$ci_bars_overlap) "OVERLAP" else "do NOT overlap",
              if (r$significant_welch) "SIGNIFICANT (p < .05)" else "not significant"))
  if (r$ci_bars_overlap && r$significant_welch)
    cat("     -> bars overlap yet p < .05: the CI of the DIFFERENCE is what the\n",
        "        test uses, and it is narrower than the two bars suggest.\n", sep = "")
  if (!is.na(r$mde_80_power))
    cat(sprintf("  at this n, 80%% power detects a difference of %.4f or larger; observed |diff| = %.4f\n",
                r$mde_80_power, abs(r$mean_diff)))
}
cat("\n")

## ---- 19e. FIGURE 1: GROUP MEANS WITH 95% CI BARS --------------------------------
## The plot the question came from: each group's mean with its own 95% CI, the
## subjects behind it drawn as points so the spread is visible, and the Welch
## p-value in the strip so the bars and the test are read together.

ht_long <- ht_data %>%
  select(Plot_Group, all_of(names(HT_SETS))) %>%
  pivot_longer(all_of(names(HT_SETS)), names_to = "measure", values_to = "value") %>%
  filter(!is.na(value))

ht_strip <- ht_tbl %>%
  transmute(measure,
            strip = sprintf("%s\nWelch p = %s, Hedges g = %.2f",
                            measure, map_chr(p_welch, fmt_p), hedges_g))

ht_points <- ht_long %>%
  left_join(ht_strip, by = "measure") %>%
  mutate(strip = factor(strip, levels = ht_strip$strip))

## One row per (measure, group): the table is stored wide (group 1 / group 2
## side by side, which is how it reads as a results table) and stacked here,
## which is what ggplot needs.
ht_wide <- ht_tbl %>%
  filter(!is.na(p_welch)) %>%
  select(measure, group_1, group_2, mean_1, mean_1_ci_lo, mean_1_ci_hi,
         mean_2, mean_2_ci_lo, mean_2_ci_hi)

ht_means <- bind_rows(
    ht_wide %>% select(measure, Plot_Group = group_1,
                       mean = mean_1, ci_lo = mean_1_ci_lo, ci_hi = mean_1_ci_hi),
    ht_wide %>% select(measure, Plot_Group = group_2,
                       mean = mean_2, ci_lo = mean_2_ci_lo, ci_hi = mean_2_ci_hi)
  ) %>%
  left_join(ht_strip, by = "measure") %>%
  mutate(strip       = factor(strip, levels = ht_strip$strip),
         Plot_Group  = factor(Plot_Group, levels = FOCUS_GROUPS))

p_means <- ggplot(ht_means, aes(x = Plot_Group, y = mean, colour = Plot_Group)) +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_jitter(data = ht_points, aes(y = value), width = 0.12, alpha = 0.35,
              size = 1.9, show.legend = FALSE) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.12, linewidth = 0.9) +
  geom_point(size = 3.4) +
  facet_wrap(~ strip, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = custom_colors, name = NULL) +
  labs(title = "5-HT4 / 5-HT6 activity by drug group",
       subtitle = paste0("Group mean +/- 95% CI, with individual subjects. ",
                         "Overlapping bars are NOT the significance test -- see Figure 2."),
       x = NULL, y = "Receptor-related brain activity") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 12, hjust = 1),
        strip.text = element_text(face = "bold", lineheight = 1.15))

ggsave(file.path(HT_DIR, "HT4_HT6_group_means_95CI.png"),
       p_means, width = 10, height = 5.2, dpi = 300)

## ---- 19f. FIGURE 2: THE DIFFERENCE, WITH ITS 95% CI -----------------------------
## The plot that actually answers the question. One row per measure: the
## hydrocodone-minus-tramadol difference and the 95% CI the Welch test is
## built on. A CI crossing the dashed zero line is a non-significant test;
## how far it extends past zero is the size of difference still compatible
## with the data.

ht_diff <- ht_tbl %>%
  filter(!is.na(p_welch)) %>%
  mutate(measure = factor(measure, levels = rev(names(HT_SETS))),
         sig     = p_welch < 0.05,
         lab     = sprintf("%+.4f  [%+.4f, %+.4f]   p = %s",
                           mean_diff, diff_ci_lo, diff_ci_hi,
                           map_chr(p_welch, fmt_p)))

p_diff <- ggplot(ht_diff, aes(x = mean_diff, y = measure)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = diff_ci_lo, xmax = diff_ci_hi),
                orientation = "y", width = 0.16, linewidth = 0.9,
                colour = "grey25") +
  geom_point(aes(shape = sig), size = 3.4, colour = "grey15", fill = "white") +
  geom_text(aes(label = lab), vjust = -1.3, size = 3.2, colour = "grey25") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21),
                     labels = c(`TRUE` = "p < .05", `FALSE` = "n.s."),
                     name = NULL, drop = FALSE) +
  expand_limits(y = c(0.5, nlevels(ht_diff$measure) + 0.75)) +
  labs(title = sprintf("Difference in receptor activity: %s minus %s",
                       sub(" group", "", FOCUS_GROUPS[1]),
                       sub(" group", "", FOCUS_GROUPS[2])),
       subtitle = paste0("Welch two-sample t-test, difference in means with 95% CI. ",
                         "A CI crossing zero = not significant at p < .05."),
       x = "Difference in receptor-related brain activity", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave(file.path(HT_DIR, "HT4_HT6_group_difference_95CI.png"),
       p_diff, width = 8.5, height = 4.2, dpi = 300)

## ---- 19g. EXPORT ---------------------------------------------------------------

## Per-subject values behind the test, so the table can be reproduced or the
## comparison re-run somewhere else without re-deriving the receptor columns.
ht_subject_values <- ht_data %>%
  select(PIN, Plot_Group, all_of(names(HT_SETS)))

write_xlsx(
  list(Group_comparison = ht_tbl,
       Subject_values   = ht_subject_values),
  file.path(HT_DIR, "Stats_HT4_HT6_group_comparison.xlsx")
)
write.csv(ht_tbl, file.path(HT_DIR, "Stats_HT4_HT6_group_comparison.csv"),
          row.names = FALSE)

message("\n5-HT4 / 5-HT6 group comparison written to: ", normalizePath(HT_DIR))

## ---- 19h. NOTES -----------------------------------------------------------------
## - Section 11 of 03_figures.R prints a t-test on the HT6 x-axis values for the
##   same two groups. That test and the 5-HT6 row here are the SAME comparison
##   on the same subjects, so the two t statistics must agree -- a useful
##   cross-check. This section adds 5-HT4, the two-receptor mean, the CIs, the
##   effect sizes, the robustness tests and the saved output that section 11
##   does not produce.
## - The 5-HT4 & 5-HT6 mean row is not an independent third test: it is a
##   function of the other two. It is here because it is the exact quantity the
##   HT6_plus_HT4 figures put on their x-axis, so a reader can check that a
##   difference seen on those figures is or is not real.
## - This is a difference in receptor values between drug groups. It says
##   nothing about whether either receptor predicts pain -- that is sections
##   15-18, and the two questions can have opposite answers. Two groups can
##   have identical mean 5-HT4 while 5-HT4 predicts pain in opposite directions
##   within each of them.
## - To test a different receptor or set, add its receptor-set name to HT_SETS
##   at 19a and re-source. Nothing else needs to change.
################################################################################
