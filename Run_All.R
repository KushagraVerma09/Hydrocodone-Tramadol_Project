###############################################################################
## Run_All.R  --  MASTER SCRIPT
##
## Runs the whole CBP / long-term opioid use pipeline end to end, in the same
## order as the original single-file script:
##   00  shared config (paths, receptor sets, colors, text sizes, labels)
##   01  load + merge the five spreadsheets into `master`
##   02  Table 1 (Hydrocodone vs Tramadol characteristics)
##   03  Figure 1 / overlay figures, one per receptor set x outcome, plus the
##       demographics table, HT6 t-test, receptor-group-means and
##       exposure-vs-outcome figures
##   04  mediation analysis (single mediator, bootstrapped)
##   05  multivariable regression: NRS ~ ROE + each receptor
##   06  parallel multi-mediator SEM (lavaan)
##   07  nested regression: NRS ~ log ROE + 5-HT4 / 5-HT6, including two
##       adjusted models (sex, MQS non-opioid, log MME, SOWS, DOU); plus the
##       same regression table with SOWS as the outcome
##   08  5-HT4 / 5-HT6 activity compared BETWEEN the two drug groups
##       (Welch t-tests, effect sizes, group-mean and difference CI figures)
##
## Any of R/0N_*.R can also be run on its own (Rscript R/0N_....R, or source()
## it directly) -- each one loads whatever it needs (config/data/table1) the
## first time it's asked for something that doesn't exist yet, so running
## them individually or all together here produces identical results.
##
## To add/rename/remove a receptor, receptor set, group, color, or text size:
## edit R/00_config.R or R/01_load_data.R ONLY. Every script below picks the
## change up automatically -- nothing else needs to be touched.
###############################################################################

PROJECT_DIR <- "/Users/kushagraverma/Work/Projects/Hydrocodone+Tramadol_Project"
R_DIR       <- file.path(PROJECT_DIR, "R")

scripts <- c("00_config.R",
             "01_load_data.R",
             "02_table1.R",
             "03_figures.R",
             "04_mediation_single.R",
             "05_multivariable_regression.R",
             "06_parallel_mediation_sem.R",
             "07_nested_regression.R",
             "08_ht4_ht6_group_ttests.R")

for (s in scripts) {
  message("\n", strrep("=", 79))
  message("Running ", s)
  message(strrep("=", 79))
  source(file.path(R_DIR, s))
}

message("\nAll done. Everything written to: ", normalizePath(OUT_DIR))
message("Optional: run Compose_Figure_Panels.R next to glue figures into ",
        "multi-panel composites.")
