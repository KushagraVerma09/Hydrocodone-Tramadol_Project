###############################################################################
## 02_table1.R
##
## Table 1: Hydrocodone vs Tramadol characteristics table, with a CBP-O
## reference column appended AFTER the p-value. Column order:
##   Characteristic | Hydrocodone | Tramadol | p-value | CBP-O
## The p-value compares the two opioid groups ONLY -- CBP-O is a descriptive
## reference column and never enters that test (build_characteristics_table's
## add_p = FALSE, plus the 2-group guard inside it, are what enforce this).
##
## Cells with no computable value ("--" -- see fmt_msd/fmt_npct/fmt_p) are
## filled gray in the exported .xlsx. write_xlsx (writexl) can't style cells,
## so the export at the bottom of this file uses openxlsx instead, purely for
## that formatting -- the table's values and column order are unchanged from
## table1_hydro_tram.
##
## Also builds the publication-formatted version (gt) and saves it as
## Publication_Table1_Hydrocodone_vs_Tramadol.html AND .docx (falling back to
## .rtf if Word export fails) -- so running this file alone produces
## everything Table 1 needs, not just the raw .xlsx.
##
## Section 11 then adds TABLES 2 AND 3 in the same house style: the
## multivariable regression models from 07_nested_regression.R, Hydrocodone and
## Tramadol side by side. Those need 07's fitted objects, so this file is no
## longer standalone-cheap -- see the dependency guard at section 11.
##
## Run standalone: Rscript R/02_table1.R  (produces Table 1's .xlsx/.html/.docx
##   plus Tables 2 and 3; section 11 sources 07_nested_regression.R if its
##   objects are not already in the session, which is slow -- 07 refits every
##   model and redraws ~30 figures)
## Sourced by: Run_All.R (which runs 07 BEFORE this file, so the guard below is
##   already satisfied and 07 is not run twice).
###############################################################################

.HTK_PROJECT_DIR <- "/Users/kushagraverma/Work/Projects/Hydrocodone+Tramadol_Project"
## Always re-sourced -- see the note above this line in 01_load_data.R.
source(file.path(.HTK_PROJECT_DIR, "R", "00_config.R"))
if (!exists("master"))              source(file.path(.HTK_PROJECT_DIR, "R", "01_load_data.R"))

## openxlsx is only needed here (for the gray-fill styling on export) -- no
## other script in the pipeline uses it, so it's loaded locally rather than
## added to the shared package list in 00_config.R, the same way
## 06_parallel_mediation_sem.R loads lavaan locally for itself alone.
if (!"openxlsx" %in% rownames(installed.packages())) install.packages("openxlsx")
library(openxlsx)

################################################################################
## 6. TABLE 1 — Hydrocodone vs Tramadol Characteristics (+ CBP-O ref)
################################################################################

fmt_msd <- function(x, d = 2) {
  x <- x[!is.na(x)]
  if (!length(x)) return("--")
  sprintf(paste0("%.", d, "f \u00b1 %.", d, "f"), mean(x), sd(x))
}
fmt_npct <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return("--")
  sprintf("%d (%.2f%%)", sum(x), 100 * mean(x))
}
p_cont <- function(df, var, gvar, groups) {
  d <- df %>% filter(.data[[gvar]] %in% groups, !is.na(.data[[var]])) %>%
    mutate(g = droplevels(factor(.data[[gvar]])))
  if (nlevels(d$g) < 2 || any(table(d$g) < 2)) return(NA_real_)
  if (nlevels(d$g) == 2) stats::t.test(d[[var]] ~ d$g)$p.value
  else summary(stats::aov(d[[var]] ~ d$g))[[1]][["Pr(>F)"]][1]
}
p_cat <- function(df, var, gvar, groups) {
  d <- df %>% filter(.data[[gvar]] %in% groups, !is.na(.data[[var]])) %>%
    mutate(g = droplevels(factor(.data[[gvar]])))
  tab <- table(d[[var]], d$g)
  if (nrow(tab) < 2 || ncol(tab) < 2) return(NA_real_)
  sup <- suppressWarnings(stats::chisq.test(tab))
  if (any(sup$expected < 5)) stats::fisher.test(tab, simulate.p.value = TRUE, B = 1e5)$p.value
  else sup$p.value
}
## fmt_p() comes from 00_config.R -- shared with every other script that
## prints a p-value, so its formatting only ever needs to change in one place.

TABLE_SPEC <- tribble(
  ~var,              ~label,                                  ~type,      ~digits,
  NA,                "Demographic and general health",         "section",  NA,
  "age",             "Age (years)",                            "cont",     2,
  "female",          "Sex/female (%)",                         "bin",      NA,
  "alcohol_yes",     "Alcohol (%)",                            "bin",      NA,
  "smoker_yes",      "Smoking (%)",                            "bin",      NA,
  "race_White",      "  White (%)",                            "bin",      NA,
  "race_Black",      "  Black (%)",                            "bin",      NA,
  "race_Other",      "  Other/Undisclosed (%)",                "bin",      NA,
  "bmi",             "BMI",                                    "cont",     2,
  "MQS_total",       "MQS total",                              "cont",     2,
  "MQS_nonopioid",   "MQS non-opioid",                         "cont",     2,
  NA,                "Concomitant medication",                 "section",  NA,
  "antidepressant",  "Antidepressants — any (%)",              "bin",      NA,
  "ad_ssri_snri",    "  SSRI/SNRI (%)",                        "bin",      NA,
  "ad_tricyclic",    "  Tricyclic/tetracyclic (%)",            "bin",      NA,
  "benzodiazepine",  "Benzodiazepines (%)",                    "bin",      NA,
  "anticonvulsant",  "Anticonvulsants (%)",                    "bin",      NA,
  "nsaid",           "NSAIDs (%)",                             "bin",      NA,
  "muscle_relaxant", "Muscle relaxants (%)",                   "bin",      NA,
  NA,                "Pain characteristics",                   "section",  NA,
  "nrs",             "Pain intensity (NRS)",                   "cont",     2,
  "PainDur",         "Pain duration (years)",                  "cont",     2,
  NA,                "Principal components",                   "section",  NA,
  "PC1",             "PC1 — Functional disability",            "cont",     2,
  "PC2",             "PC2 — Pain quality/severity",            "cont",     2,
  "PC3",             "PC3 — Negative affect",                  "cont",     2,
  NA,                "Opioid consumption and behavior",        "section",  NA,
  "MME",             "MME (mg/day)",                           "cont",     2,
  "ROE",             "ROE (mg/L)",                             "cont",     4,
  "DOU",             "Duration of opioid use (years, DOU)",    "cont",     2,
  "SOWS",            "SOWS",                                   "cont",     2,
  "COMM",            "COMM",                                   "cont",     2
)

build_characteristics_table <- function(df, gvar = "Group", groups = c("CBP-O", "CBP+O"),
                                        add_p = TRUE, display_names = NULL) {
  d <- df %>% filter(.data[[gvar]] %in% groups) %>%
    mutate(race_White = race_cat == "White",
           race_Black = race_cat == "Black",
           race_Other = race_cat == "Other/Undisclosed",
           .g = factor(.data[[gvar]], levels = groups)) %>% droplevels()

  ## The p-value column must be a two-group comparison and nothing else. `d` is
  ## already filtered to `groups`, so p_cont/p_cat below can only ever see those
  ## rows -- reference columns are produced by separate add_p = FALSE calls and
  ## are never pooled into a test. This guard keeps that true if the call sites
  ## are ever edited: with more than two groups p_cont silently switches to a
  ## one-way ANOVA, which is not the requested test.
  if (add_p && nlevels(d$.g) != 2)
    stop("add_p = TRUE requires exactly 2 groups; got ",
         nlevels(d$.g), ": ", paste(levels(d$.g), collapse = ", "))

  ns <- table(d$.g)
  ## display_names lets a raw data value (e.g. Group == "H") show up under a
  ## readable header ("Healthy") without changing what's actually filtered on.
  ## Every other caller leaves this NULL and gets the raw group value back,
  ## same as before.
  disp <- if (is.null(display_names)) names(ns) else unname(display_names[names(ns)])
  hdr <- sprintf("%s (n=%d)", disp, as.integer(ns))

  rows <- pmap(TABLE_SPEC, function(var, label, type, digits) {
    if (type == "section") {
      out <- as.list(c(label, rep("", length(groups)), ""))
    } else {
      cells <- map_chr(groups, function(g) {
        v <- d[[var]][d$.g == g]
        if (type == "cont") fmt_msd(v, digits) else fmt_npct(v)
      })
      p <- if (!add_p) NA_real_ else
        if (type == "cont") p_cont(d, var, ".g", groups) else p_cat(d, var, ".g", groups)
      out <- as.list(c(label, cells, fmt_p(p)))
    }
    names(out) <- c("Characteristic", hdr, "p-value")
    as_tibble(out)
  })
  bind_rows(rows)
}

# 1. Build the target table: Hydrocodone vs Tramadol ONLY, with a p-value.
#    build_characteristics_table filters to FOCUS_GROUPS before p_cont/p_cat
#    ever run, and its add_p = TRUE guard requires exactly 2 groups -- CBP-O
#    (built next) cannot reach this test, by construction, not just by
#    convention. Column order out of this call is already
#    Characteristic | Hydrocodone | Tramadol | p-value, which is exactly
#    where the p-value belongs, so nothing below needs to reorder it.
table1_hydro_tram <- build_characteristics_table(
  master, "Plot_Group", FOCUS_GROUPS
)

# 2. Reference column, descriptive only (add_p = FALSE: no test, so this
#    never influences or is influenced by the Hydrocodone/Tramadol p-value
#    above).
table1_cbpo <- build_characteristics_table(
  master, "Group", c("CBP-O"), add_p = FALSE
)

# 3. Append the reference column AFTER the p-value:
#    Characteristic | Hydrocodone | Tramadol | p-value | CBP-O
table1_hydro_tram <- bind_cols(
  table1_hydro_tram,
  table1_cbpo[, 2]        # "CBP-O (n=...)"
)

cat("\n===== TABLE 1: Hydrocodone vs Tramadol (with CBP-O reference) =====\n")
print(as.data.frame(table1_hydro_tram), row.names = FALSE)

## ---- Export, with "--" (no computable value) cells filled gray --------------
## write_xlsx() can't style cells, so this uses openxlsx purely for that fill --
## the values and column order are exactly table1_hydro_tram above.
tbl_out  <- as.data.frame(table1_hydro_tram)
na_style <- createStyle(fgFill = "#D9D9D9")

wb <- createWorkbook()
addWorksheet(wb, "Hydrocodone_vs_Tramadol")
writeData(wb, "Hydrocodone_vs_Tramadol", tbl_out)

for (col in seq_along(tbl_out)) {
  na_rows <- which(tbl_out[[col]] == "--")
  if (length(na_rows))
    ## +1: row 1 of the sheet is the header written by writeData(); data row i
    ## of tbl_out lands on sheet row i + 1.
    addStyle(wb, "Hydrocodone_vs_Tramadol", style = na_style,
            rows = na_rows + 1, cols = col, gridExpand = TRUE, stack = TRUE)
}

saveWorkbook(wb, file.path(OUT_DIR, "Table1_Hydrocodone_vs_Tramadol.xlsx"), overwrite = TRUE)


################################################################################
## 10. FORMAT PUBLICATION-READY DEMOGRAPHICS TABLE  ->  outputs/
################################################################################
## Moved here from 03_figures.R (where this was originally section 10): the
## publication-formatted .html/.docx is Table 1's own output, not a figure, so
## it belongs with Table 1's .xlsx -- and running this file alone
## (Rscript R/02_table1.R) now produces the .xlsx AND the .html/.docx, instead
## of needing 03_figures.R just for the docx.

table_data <- read_excel(file.path(OUT_DIR, "Table1_Hydrocodone_vs_Tramadol.xlsx"))

pub_table <- table_data %>%
  gt() %>%
  cols_align(align = "left", columns = Characteristic) %>%
  cols_align(align = "center", columns = -Characteristic) %>%
  ## Column order is Characteristic | Hydrocodone | Tramadol | p-value | CBP-O,
  ## matching table1_hydro_tram above. Adding, removing, or reordering a
  ## column there means updating the indices below too.
  tab_spanner(label = "Opioid subgroup comparison", columns = 2:3) %>%
  tab_spanner(label = "Reference group",            columns = 5) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything())
  ) %>%
  tab_style(
    style = cell_text(weight = "bold", style = "italic"),
    locations = cells_body(
      columns = Characteristic,
      rows = `p-value` == "" | is.na(`p-value`)
    )
  ) %>%
  tab_footnote(
    footnote = paste0("p-values compare ", FOCUS_GROUPS[1], " vs ", FOCUS_GROUPS[2],
                      " only. The CBP-O column is a descriptive reference ",
                      "group and is not included in the test."),
    locations = cells_column_labels(columns = `p-value`)
  ) %>%
  tab_options(
    table.border.top.color = "black",
    table.border.top.width = px(2),
    table.border.bottom.color = "black",
    table.border.bottom.width = px(2),
    column_labels.border.top.color = "black",
    column_labels.border.top.width = px(2),
    column_labels.border.bottom.color = "black",
    column_labels.border.bottom.width = px(2),
    table_body.border.bottom.color = "black",
    table_body.border.bottom.width = px(2),
    table_body.hlines.color = "transparent",
    table_body.vlines.color = "transparent"
  )

gtsave(pub_table, file.path(OUT_DIR, "Publication_Table1_Hydrocodone_vs_Tramadol.html"))

tryCatch(
  gtsave(pub_table, file.path(OUT_DIR, "Publication_Table1_Hydrocodone_vs_Tramadol.docx")),
  error = function(e) {
    message("gtsave(.docx) failed (", conditionMessage(e),
            ") — writing .rtf instead; Word opens it natively.")
    gtsave(pub_table, file.path(OUT_DIR, "Publication_Table1_Hydrocodone_vs_Tramadol.rtf"))
  }
)

message("\nPublication table successfully formatted and saved to: ", normalizePath(OUT_DIR))


################################################################################
## 11. TABLES 2 & 3 — Multivariable regression, Hydrocodone vs Tramadol
################################################################################
## Two publication tables built from 07_nested_regression.R's models, in the
## same house style as Table 1 above (openxlsx .xlsx with gray "--" cells, plus
## a gt .html/.docx). Table 2 predicts NRS pain, Table 3 predicts SOWS
## withdrawal.
##
## Only THREE of 07's five models appear, by request -- M3, M4 and M5:
##   HT4_HT6             M3  y ~ log_ROE + R_5HT4 + R_5HT6
##   HT4_HT6_ROEint      M4  + (log_ROE x R_5HT4) + (log_ROE x R_5HT6)
##   HT4_HT6_ROEint_Cov  M5  + sex + log10 MME + MQS non-opioid
## M1 (log_ROE + R_5HT4) and M2 (log_ROE + R_5HT6) are the single-receptor
## models; they are fit and reported in 07's own figures and .xlsx, but they do
## not get a column here.
##
## NOTE on the n row: M5 adds log_MME, which is missing for 5 subjects, so its
## column is fit on fewer rows than M3 and M4 and its adj. R2 cannot be read
## against theirs. The n row is in the table precisely so that is visible --
## see 07_nested_regression.R section 18m.
##
## For SOWS the models come from MREG_MODELS_SOWS, which drops SOWS from the
## right-hand side if it is ever there. With the current M1-M5 set no model
## contains SOWS, so Table 3 simply has no "SOWS withdrawal" predictor row.

## 07's fitted objects, not its .xlsx. Reading
## New_Outputs/multiple_regression_HT4_HT6/Stats_*.xlsx would be cheaper, but it
## would silently serve the PREVIOUS run's numbers any time the models change,
## because that file is written by a script that has not run yet in this
## session. Same self-bootstrapping idiom as the config/data guards at the top.
if (!exists("mreg_coef_tbl") || !exists("mreg_fit_tbl") || !exists("mreg_sows"))
  source(file.path(.HTK_PROJECT_DIR, "R", "07_nested_regression.R"))

## The three models, in column order. Named by KEY, not position, so reordering
## MREG_MODELS in 07 cannot silently change which models this table shows.
REGTAB_MODELS <- c("HT4_HT6", "HT4_HT6_ROEint", "HT4_HT6_ROEint_Cov")

## Builds one wide table: predictor rows, then one column per (group, model),
## then the n / adj. R2 / model p footer rows.
##
## Returns a list rather than a bare data frame because the exporter needs to
## know which columns belong to which group (for the spanners) and what the
## repeated header labels are -- gt column IDs have to be unique, so the frame's
## own names cannot carry that.
build_regression_table <- function(coef_tbl, fit_tbl,
                                   models     = REGTAB_MODELS,
                                   groups     = FOCUS_GROUPS,
                                   shorts     = MREG_MODEL_SHORT,
                                   model_defs = MREG_MODELS) {

  ## One column per group x model, group-major (all of Hydrocodone's models,
  ## then all of Tramadol's). expand.grid varies its FIRST argument fastest,
  ## which is what puts model inside group.
  grid <- expand.grid(model = models, group = groups, stringsAsFactors = FALSE)
  grid$id <- sprintf("G%d_M%d", match(grid$group, groups),
                                match(grid$model, models))

  ## Same term ordering as 07's table figure (07_nested_regression.R:816):
  ## intercept first, then predictors in the order the models introduce them.
  terms_order <- c("(Intercept)", unique(unlist(model_defs[models])))

  ## "--" for a term the model does not contain, matching Table 1's convention
  ## for a cell with no computable value (and picked up by the gray fill on
  ## export). Coefficient + stars only, no standard error -- the .xlsx written
  ## by 07 has se, CIs and standardized betas for anyone who needs them.
  cell <- function(g, mk, tm) {
    r <- coef_tbl[coef_tbl$group == g & coef_tbl$model == mk &
                  coef_tbl$term == tm, ]
    if (!nrow(r)) "--" else sprintf("%.3f%s", r$est[1], mreg_stars(r$p[1]))
  }
  fcell <- function(g, mk, fn) {
    r <- fit_tbl[fit_tbl$group == g & fit_tbl$model == mk, ]
    if (!nrow(r)) "--" else fn(r)
  }

  row_of <- function(label, f) {
    out <- as.list(c(label, unlist(Map(f, grid$group, grid$model))))
    names(out) <- c("Predictor", grid$id)
    as_tibble(out)
  }

  body <- map(terms_order, function(tm)
    ## mreg_lab() gives the pretty name ("5-HT4 activity"); it has no entry for
    ## the intercept, so that one is spelled out, exactly as 07's figure does.
    row_of(if (tm == "(Intercept)") "Intercept" else mreg_lab(tm),
           function(g, mk) cell(g, mk, tm)))
  body <- bind_rows(body)

  ## Drop any predictor row that is "--" in all six columns. This is what
  ## removes the SOWS row from Table 3 without a second term list.
  keep <- apply(body[, grid$id, drop = FALSE], 1, function(r) any(r != "--"))
  body <- body[keep, , drop = FALSE]

  foot <- bind_rows(
    row_of("n",       function(g, mk) fcell(g, mk, function(r) sprintf("%d", r$n[1]))),
    row_of("adj. R2", function(g, mk) fcell(g, mk, function(r) sprintf("%.3f", r$adj_R2[1]))),
    ## Model-level p: the overall F test, not a per-coefficient p.
    row_of("model p", function(g, mk) fcell(g, mk, function(r) mreg_fmt_p(r$F_p[1])))
  )

  list(tbl    = bind_rows(body, foot),
       grid   = grid,
       groups = groups,
       n_body = nrow(body),
       ## MREG_MODEL_SHORT is hand-wrapped with "\n" for ggplot; neither gt nor
       ## an .xlsx header cell wants a literal newline there.
       hdr    = gsub("\n", " ", mreg_model_short(models, shorts)))
}

## Exports one built table both ways, step for step as Table 1 above: openxlsx
## for the gray-filled .xlsx, gt for the .html/.docx with .rtf fallback.
export_regression_table <- function(rt, base, title, subtitle) {
  out <- as.data.frame(rt$tbl)
  n_g <- length(rt$groups)
  n_m <- length(rt$hdr)

  ## ---- .xlsx ----------------------------------------------------------------
  ## The model header labels repeat across the two groups, and a data frame
  ## cannot carry duplicate names -- so the .xlsx qualifies each one with its
  ## group instead of leaning on a spanner row it has no way to draw.
  xl <- out
  names(xl) <- c("Predictor", sprintf("%s: %s", rt$grid$group, rep(rt$hdr, n_g)))

  wb2 <- createWorkbook()
  addWorksheet(wb2, "Regression")
  writeData(wb2, "Regression", xl)
  for (col in seq_along(xl)) {
    na_rows <- which(xl[[col]] == "--")
    if (length(na_rows))
      ## +1 for the header row written by writeData(), as above.
      addStyle(wb2, "Regression", style = na_style,
               rows = na_rows + 1, cols = col, gridExpand = TRUE, stack = TRUE)
  }
  saveWorkbook(wb2, file.path(OUT_DIR, paste0(base, ".xlsx")), overwrite = TRUE)

  ## ---- publication .html / .docx --------------------------------------------
  ## Column 1 is the predictor label; each group's block is the n_m columns
  ## after it. Computed rather than hardcoded as 2:4 / 5:7 so adding a fourth
  ## model to REGTAB_MODELS does not need an edit here.
  gt_tbl <- gt(out) %>%
    cols_align(align = "left",   columns = Predictor) %>%
    cols_align(align = "center", columns = -Predictor) %>%
    cols_label(.list = setNames(as.list(rep(rt$hdr, n_g)), rt$grid$id)) %>%
    tab_header(title = title, subtitle = subtitle)

  for (i in seq_len(n_g))
    gt_tbl <- gt_tbl %>%
      tab_spanner(label   = rt$groups[i],
                  columns = 1 + (i - 1) * n_m + seq_len(n_m))

  gt_tbl <- gt_tbl %>%
    tab_style(style     = cell_text(weight = "bold"),
              locations = cells_column_labels(everything())) %>%
    ## The n / adj. R2 / model p block is model-level, not a coefficient --
    ## italicised to read as a footer, the same way 07's table figure sets
    ## those rows in italic.
    tab_style(style     = cell_text(style = "italic"),
              locations = cells_body(rows = seq_len(nrow(out)) > rt$n_body)) %>%
    tab_footnote(
      footnote  = paste0("Raw-unit regression coefficients; *** p < .001, ",
                         "** p < .01, * p < .05. \"--\" marks a predictor the ",
                         "model does not contain. Raw coefficients are not ",
                         "comparable across predictors on different scales -- ",
                         "see the 07 output folder for standardized betas and ",
                         "confidence intervals."),
      locations = cells_column_labels(columns = Predictor)) %>%
    tab_footnote(
      footnote  = paste0("Models are fit separately within each drug group, ",
                         "not pooled with a group term. \"model p\" is the ",
                         "overall F test for that column's model."),
      locations = cells_column_labels(columns = rt$grid$id[1])) %>%
    tab_options(
      table.border.top.color            = "black",
      table.border.top.width            = px(2),
      table.border.bottom.color         = "black",
      table.border.bottom.width         = px(2),
      column_labels.border.top.color    = "black",
      column_labels.border.top.width    = px(2),
      column_labels.border.bottom.color = "black",
      column_labels.border.bottom.width = px(2),
      table_body.border.bottom.color    = "black",
      table_body.border.bottom.width    = px(2),
      table_body.hlines.color           = "transparent",
      table_body.vlines.color           = "transparent"
    )

  gtsave(gt_tbl, file.path(OUT_DIR, paste0("Publication_", base, ".html")))
  tryCatch(
    gtsave(gt_tbl, file.path(OUT_DIR, paste0("Publication_", base, ".docx"))),
    error = function(e) {
      message("gtsave(.docx) failed (", conditionMessage(e),
              ") — writing .rtf instead; Word opens it natively.")
      gtsave(gt_tbl, file.path(OUT_DIR, paste0("Publication_", base, ".rtf")))
    }
  )
  invisible(out)
}

## ---- Table 2: NRS -----------------------------------------------------------
table2_nrs <- build_regression_table(mreg_coef_tbl, mreg_fit_tbl)

cat("\n===== TABLE 2: NRS pain, multivariable regression =====\n")
print(as.data.frame(table2_nrs$tbl), row.names = FALSE)

export_regression_table(
  table2_nrs, "Table2_Regression_NRS_HT4_HT6",
  title    = "Table 2. NRS pain regressed on log10 ROE and 5-HT4 / 5-HT6 activity",
  subtitle = "Raw-unit coefficients, fit separately within each drug group")

## ---- Table 3: SOWS ----------------------------------------------------------
## MREG_MODELS_SOWS / MREG_MODEL_SHORT_SOWS are the SOWS-specific model list and
## headers from 07 section 18l -- SOWS is the outcome here, so it is not on the
## right-hand side and must not be advertised in a header.
table3_sows <- build_regression_table(
  mreg_sows$coef, mreg_sows$fit,
  shorts     = MREG_MODEL_SHORT_SOWS,
  model_defs = MREG_MODELS_SOWS)

cat("\n===== TABLE 3: SOWS withdrawal, multivariable regression =====\n")
print(as.data.frame(table3_sows$tbl), row.names = FALSE)

export_regression_table(
  table3_sows, "Table3_Regression_SOWS_HT4_HT6",
  title    = "Table 3. SOWS withdrawal regressed on log10 ROE and 5-HT4 / 5-HT6 activity",
  subtitle = "Raw-unit coefficients, fit separately within each drug group")

message("\nTables 2 and 3 (regression) saved to: ", normalizePath(OUT_DIR))
