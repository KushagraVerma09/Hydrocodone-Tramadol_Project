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
## Run standalone: Rscript R/02_table1.R  (produces the .xlsx, .html and
##   .docx together -- nothing else in the pipeline needs to run first)
## Sourced by: Run_All.R.
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
