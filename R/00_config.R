###############################################################################
## 00_config.R  --  SHARED CONFIGURATION
##
## This is the ONE place that defines paths, the receptor vocabulary, group
## names, figure styling, and small formatting helpers used across the whole
## pipeline. Every other script in R/ sources this file before doing anything
## else, so a change made HERE (a receptor added/renamed, a color, a text
## size, a group label) is picked up by every figure, table and test the
## pipeline produces. Nothing needs to be edited in more than one place.
## 
## This file produces no output on its own. It is sourced by:
##   01_load_data.R, 02_table1.R, 03_figures.R, 04_med iation_single.R,
##   05_multivariable_regression.R, 06_parallel_mediation_sem.R,
##   07_nested_regression.R, Run_All.R, and ../Compos e_Figure_Panels.R
##
## Split out of the original Hydrocodone_Tramadol_Kush_All_Graphs.R (sections
## 0, 1, 1b) plus a few pieces of shared config that used to live inline in
## later sections (fmt_p, custom_colors/custom_shapes, STATS_MODES,
## MED_LABELS) and have been centralized here so every script sees the same
## values instead of redefining or falling back to defaults.
###############################################################################

## NOTHING in this file is gated on "already loaded" except the one thing that
## genuinely must run once (set.seed, below). That is deliberate, and it has bitten
## this file twice.
##
## Once .HTK_CONFIG_LOADED was TRUE, re-sourcing this file -- directly, or via the
## guard at the top of every other script -- was a no-op, so an edit to a value
## here only took effect after restarting R. Every styling constant is therefore
## always re-evaluated on every source(): edit one and the next script to
## re-source this file picks it up, no restart. The second bite was the package
## list; see the note above the attach loop.

## ---- 0. PACKAGES -------------------------------------------------------------
## ggtext supplies element_markdown(), which is what lets one facet strip carry
## the group name and the stats line at two DIFFERENT sizes. See strip_md() in
## 03_figures.R.
pkgs <- c("readxl", "writexl", "dplyr", "tidyr", "stringr",
          "purrr", "tibble", "ggplot2", "broom", "gt", "ggtext")

## Attaching runs on EVERY source(), deliberately, and is NOT inside the guard
## below.
##
## It used to be, and that was a trap: the guard fires once per R session, so a
## package ADDED to the list above was never attached in any session that had
## already sourced this file -- the whole block was skipped and the new function
## simply did not exist ("could not find function element_markdown"). A fresh
## Rscript never sees it, because there the guard is always cold, so the failure
## only shows up in a long-running session (RStudio) and only for whoever added
## the package. This is the same staleness the comment below describes, and it
## belongs in the same bin.
##
## requireNamespace() is a cheap already-installed check, and library() on an
## attached package is a no-op, so paying this on every source() costs nothing.
## The expensive part -- installed.packages(), which scans every library path --
## is what stays guarded, further down.
invisible(lapply(pkgs, function(p) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}))

## Guard: this ONLY wraps work that must happen exactly once per session.
## set.seed() is in here rather than above because every script re-sources this
## file, and re-seeding at each of those points would restart the random stream
## mid-pipeline and change the bootstrap results.
if (!exists(".HTK_PKGS_LOADED")) {

set.seed(42)

.HTK_PKGS_LOADED <- TRUE

}  ## end once-per-session guard

## ---- 1. PATHS ------------------------------------------------------------
## PROJECT_DIR is the only hardcoded absolute path in the whole pipeline --
## every other script derives its paths from this one via R_DIR/OUT_DIR.
PROJECT_DIR <- "/Users/kushagraverma/Work/Projects/Hydrocodone+Tramadol_Project"
DATA_DIR    <- PROJECT_DIR   # kept as DATA_DIR too: original scripts' name for this
R_DIR       <- file.path(PROJECT_DIR, "R")

## Everything the pipeline writes goes here.
OUT_DIR  <- file.path(DATA_DIR, "New_Outputs")
dir.create(OUT_DIR, showWarnings = FALSE)

F_BEHAVIOR <- file.path(DATA_DIR, "Copy of data_behavior_all_cbp.xlsx")
F_RECEPTOR <- file.path(DATA_DIR, "Copy of data_receptor_AI_all_subjects.xlsx")
F_PCA      <- file.path(DATA_DIR, "PCA_results_CBP_all.xlsx")
F_MQS      <- file.path(DATA_DIR, "mqs final.xlsx")
F_MED      <- file.path(DATA_DIR, "combinedpatients medication.xlsx")

## ---- 2. X-AXIS MODE -----------------------------------------------------------
##   "raw"        -> plot the receptor value as given. No anchoring, no
##                   rescaling, no division by a group difference. For a
##                   single-receptor set this is literally that receptor's
##                   column; for a multi-receptor set it is the mean across
##                   the receptors in the set.
##   "projection" -> old behaviour: 0 = healthy mean, 1 = CBP-O mean
##   "minmax"     -> old behaviour: rescaled to [0,1] across the sample
X_MODE <- "raw"
stopifnot(X_MODE %in% c("raw", "projection", "minmax"))

## Healthy reference line (green dotted). Drawn at the healthy group mean,
## which under X_MODE = "raw" sits at essentially zero.
##
## OFF: the line and its caption sentence were taking vertical space in every
## panel for a reference that sits at ~0 under X_MODE = "raw". Everything that
## draws or describes it is gated on this flag, so setting it back to TRUE
## restores both the line and its caption text with no other edit.
SHOW_HEALTHY_LINE <- FALSE
HEALTHY_LINE_COLOR <- "#3B6D11"

MIN_N_SLOPE <- 5   # groups smaller than this are dropped from slope models

## ---- 3. FIGURE TEXT SIZING -----------------------------------------------------
## Every text element in every figure scales off TEXT_SCALE. The figures are
## written large and then shrunk when several are composed onto one page
## (see Compose_Figure_Panels.R), so the type has to be oversized here to stay
## readable there. Set TEXT_SCALE <- 1 to get the original, smaller styling back.
##
## The sizes are deliberately NOT all equal, and only two things are bold. The
## intended reading order is
##   group name  >  axis titles  >  tick labels / strip stats  >  caption
## and bold is reserved for the group name and the in-panel stats box (the one
## piece of text sitting on top of the plot, where it needs the extra weight to
## stay legible). Everything that used to be bold-and-large -- caption, strip
## stats -- is now plain, which is what stops the figures reading as uniformly
## shouty.
## TEXT_SCALE 4 overshot. At that scale two things break rather than just look
## heavy: the x tick labels of adjacent facet panels run into each other into an
## unreadable smear, and the scatter points merge into blobs that hide where the
## data actually is. TEXT_SCALE 2 (the value before the bump) is the other
## extreme -- hairline fit lines and near-invisible axis rules. 3 is the rung in
## between, and is the last one where the tick labels stay clear. The ladder the
## choice was made from is rendered in New_Outputs/_size_previews/.
TEXT_SCALE       <- 3
FIG_BASE_SIZE    <- 11   * TEXT_SCALE   # theme base size            (33pt)
## The y-axis title is ROTATED, so its length is spent against the figure's
## HEIGHT, and the shortest figure in use is the 6.8 in exposure overlay whose
## panel is only about 4.3 in tall. "Clinical Pain (NRS 0-10)" measures 5.06 in
## at 33pt and 4.14 in at 27pt, so 27 is the largest size that does not get its
## last characters clipped off the top of that figure. (At the old TEXT_SCALE 4
## it was 48pt / 7.36 in -- longer than the entire figure.)
FIG_AXIS_TITLE   <- 9    * TEXT_SCALE   # x / y axis titles          (27pt)
FIG_AXIS_TEXT    <- 8    * TEXT_SCALE   # tick labels                (24pt)
FIG_STRIP_SIZE   <- 10   * TEXT_SCALE   # facet strip: group name, bold (30pt)
## This one is LIVE -- change it and the n / r / slope line really does resize,
## and only that line. It used to be inert: the strip was built with plotmath,
## whose only size commands are the fixed ratios scriptstyle (~0.71x) and
## scriptscriptstyle (~0.5x), so the stats line was pinned to the group name's
## size and no value here could move it. The strip is HTML now (strip_md() in
## 03_figures.R), which resolves a real font-size per line.
##
## Set it equal to FIG_STRIP_SIZE for identical sizes; 9 vs 10 leaves the stats
## line a touch smaller so the group name still leads. Oversized values are
## capped to what fits the strip rather than clipped -- the run prints the size
## it actually used when that happens.
FIG_STRIP_STATS  <- 9    * TEXT_SCALE   # n / r / slope line beneath it (27pt)
## NOT bumped, deliberately. The caption is now a single horizontal line, and at
## 18pt the logROE caption measures 11.7 in against an 11 in overlay figure, so it
## had to fold onto two lines -- which is exactly the vertical space the one-line
## rewrite was meant to save. At 16pt it measures 10.4 in and fits on one line on
## every figure size in use (13/12/11 in). The caption is also the text that least
## needs to be readable from a distance; the ticks and titles got the lift instead.
FIG_CAPTION_SIZE <- 5    * TEXT_SCALE   # t-test / slope caption under the plot
FIG_ANNOT_SIZE   <- 2.7  * TEXT_SCALE   # in-panel stats box (geom_label, mm)
FIG_POINT_SIZE   <- 3.0  * TEXT_SCALE   # scatter point size
FIG_LINE_SIZE    <- 1.2  * TEXT_SCALE   # regression / best-fit line width
FIG_CAPTION_COL  <- "grey25"            # caption colour: present, not shouting

## ---- 3b. FLOATING-AXIS STYLING + IN-PANEL ANNOTATION -------------------------
## The figures use "floating" axes: no gridlines and no panel box, just a bold
## bottom and left axis line with outward-facing ticks, so nothing competes with
## the data. theme_fig() in 03_figures.R is the only consumer.
FIG_FONT          <- "sans"            # ONE family for theme AND annotations
## These two are absolute, NOT scaled by TEXT_SCALE -- ggplot linewidths are in
## mm, so they do not need to track the point sizes, and at 20/8 they rendered as
## solid black slabs thicker than the fit line they were framing.
FIG_AXIS_LINE     <- 4.0               # bottom/left axis line width
FIG_TICK_LINE     <- 2.0               # tick mark width
FIG_TICK_LEN      <- unit(0.22, "cm")  # POSITIVE = ticks point OUTWARD
FIG_AXIS_TEXT_COL <- "grey15"          # tick labels: near-black, high contrast

## The in-panel stats box is measured by label_size_in() to decide which corner
## it goes in, and drawn by geom_label(). Those are two different pieces of code
## that MUST agree about the box's size -- if the drawn lineheight/padding/face
## differs from the measured one, box_fracs() computes a box that doesn't match
## what's rendered and the box starts covering data points. They used to be
## hardcoded literals in both places and kept in sync by hand; naming them here
## is what makes that class of bug impossible.
FIG_ANNOT_LINEHEIGHT <- 0.90           # tight block, was 1.05
FIG_ANNOT_PAD        <- 0.25           # label.padding, in lines; was 0.35
FIG_ANNOT_FACE       <- "bold"
FIG_ANNOT_COL        <- "grey15"       # distinct dark grey, not pure black

## ---- 4. GROUPS + COLORS ---------------------------------------------------
## The ONLY two groups that appear in any figure, table or model.
FOCUS_GROUPS <- c("Hydrocodone group", "Tramadol group")

## Both groups are filled circles. Shape therefore carries no information and
## its legend is switched off wherever it would otherwise be drawn; colour alone
## separates the groups.
custom_colors <- setNames(c("#1b9e77", "#d95f02")[seq_along(FOCUS_GROUPS)], FOCUS_GROUPS)
custom_shapes <- setNames(rep(16, length(FOCUS_GROUPS)), FOCUS_GROUPS)

## Group names for FIGURES ONLY.
##
## The DATA keeps the long names: Plot_Group's factor levels (01_load_data.R),
## every `filter(Plot_Group %in% FOCUS_GROUPS)`, and the NAMES on custom_colors /
## custom_shapes above are all keyed on them, and ggplot matches a manual scale's
## values by those names. So this map is applied at RENDER time in 03_figures.R
## (facet strips, legend keys, captions) and nowhere else -- which is why
## Table 1 (02) and the mediation/regression/t-test scripts (04-08) keep printing
## the full names and cannot break.
GROUP_LABELS <- setNames(c("CBP+Hydrocodone", "CBP+Tramadol")[seq_along(FOCUS_GROUPS)], FOCUS_GROUPS)

## Vectorized: long group name (or factor) -> short label, anything unrecognised
## passes through unchanged.
glab <- function(g) {
  hit <- unname(GROUP_LABELS[as.character(g)])
  ifelse(is.na(hit), as.character(g), hit)
}

## ---- 5. WHERE THE n / r / slope / p TEXT GOES ON A FIGURE ---------------------
## Two placements, both produced on every run into parallel folders that hold
## identical filenames, so the same figure can be compared side by side:
##
##   "header" -> the stats become a second line on the strip above each panel
##               (and the legend key text on the overlay figures). The text is
##               outside the panel, so it cannot touch the data at all.
##   "corner" -> the stats stay inside the panel, in a box placed in whichever
##               corner is emptiest.
##
## Compose_Figure_Panels.R reads these same folder names when it glues figures
## into multi-panel composites -- defined once here so the two scripts can
## never drift out of sync.
STATS_MODES <- c(header = "figs_stats_header", corner = "figs_stats_corner")

## ---- 6. RECEPTOR LABELS FOR FIGURES/TABLES -----------------------------------
## Pretty labels used by the mediation / regression / SEM scripts (04-07) when
## printing axis titles, path diagrams and tables. Add a receptor or predictor
## here and every script that formats it picks the new label up automatically.
MED_LABELS <- c(log_ROE = "log10 ROE",
                log_MME = "log10 MME",
                R_MOR   = "MOR activity",
                R_5HT1A = "5-HT1A activity",
                R_5HT4  = "5-HT4 activity",
                R_5HT6  = "5-HT6 activity",
                nrs     = "NRS pain",
                PC2     = "PC2 pain quality",
                ## Exposure x receptor interaction terms (07_nested_regression.R
                ## models M4/M5) -- added here so their axis titles/figure
                ## titles read as English instead of the raw column name.
                R_5HT4_x_ROE = "5-HT4 x log10 ROE",
                R_5HT6_x_ROE = "5-HT6 x log10 ROE")
med_lab <- function(v) if (!is.na(MED_LABELS[v])) unname(MED_LABELS[v]) else v

## ---- 7. SHARED FORMATTING HELPERS ---------------------------------------------
## p-value formatter used everywhere a single p-value is printed (tables,
## captions, in-panel stats boxes). Vectorized formatters that some sections
## need instead (e.g. formatting a whole column at once) stay local to those
## sections -- see pmed_fmt_p in 06_parallel_mediation_sem.R for why.
fmt_p <- function(p) {
  if (is.na(p)) return("--")
  if (p < 1e-5) return("<10⁻⁵")
  if (p < 0.001) return(sprintf("%.1e", p))
  sprintf("%.3f", p)
}

## Kept only as a flag OTHER code can check to confirm this file has been
## sourced (nothing left in this file gates on it -- see the note at the
## top). The per-script guards that used to read `if (!exists(".HTK_CONFIG_
## LOADED")) source("00_config.R")` were removed for the same reason: that
## guard was what made an edit here invisible until R was restarted, because
## it skipped re-sourcing entirely once this flag was set. 01_load_data.R,
## 02_table1.R, 03_figures.R, 04_mediation_single.R,
## 05_multivariable_regression.R, 06_parallel_mediation_sem.R,
## 07_nested_regression.R, and 08_ht4_ht6_group_ttests.R now source this file
## unconditionally -- cheap, since the only expensive part (packages, above)
## still guards itself.
.HTK_CONFIG_LOADED <- TRUE
