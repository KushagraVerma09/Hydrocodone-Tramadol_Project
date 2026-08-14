###############################################################################
## CBP / LONG-TERM OPIOID USE — MASTER ANALYSIS SCRIPT
##
## What this script does, in order:
##   1. Loads all five spreadsheets and MERGES them into ONE master dataset.
##      Every figure and table downstream is built from that single object.
##   2. Derives the medication flags (antidepressant use) and the opioid
##      sub-groups (Plot_Group).
##   3. Builds a Table-S1-style characteristics table comparing
##      Hydrocodone group vs Tramadol group, with a CBP+O reference column
##      appended. The p-value compares the two opioid groups ONLY -- the
##      reference column is descriptive and never enters a test.
##   4. Builds Figure 1 & overlays for multiple RECEPTOR SUBSETS (e.g., MOR,
##      HT1A+HT2A, etc.) and clinical outcomes (NRS, PC1, PC2, PC3).
##      *** ALL FIGURES ARE RESTRICTED TO Hydrocodone group vs Tramadol group ***
##      *** X = RAW receptor-related brain activity (no 0/1 rescaling)      ***
##      *** EVERY FIGURE DISPLAYS THE X-AXIS T-TEST BETWEEN THE TWO GROUPS  ***
##   5. Tests whether the regression slopes actually differ between groups
##      (global interaction test + pairwise slope contrasts, FDR corrected).
##   6. Writes everything to ./outputs/ with clearly labeled filenames.
################################################################################

## ---- 0. PACKAGES -------------------------------------------------------------
pkgs <- c("readxl", "writexl", "dplyr", "tidyr", "stringr",
          "purrr", "tibble", "ggplot2", "broom", "gt")
missing <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(missing)) install.packages(missing)
invisible(lapply(pkgs, library, character.only = TRUE))

set.seed(42)

## ---- 1. CONFIG ---------------------------------------------------------------
DATA_DIR <- "/Users/kushagraverma/Work/Projects/Hydrocodone+Tramadol_Project" # <- folder holding the five .xlsx files
## Everything this script writes goes here. Kept separate from the older
## ./outputs/ folder so the previous run stays intact for side-by-side comparison.
OUT_DIR  <- file.path(DATA_DIR, "New_Outputs")
dir.create(OUT_DIR, showWarnings = FALSE)

F_BEHAVIOR <- file.path(DATA_DIR, "Copy of data_behavior_all_cbp.xlsx")
F_RECEPTOR <- file.path(DATA_DIR, "Copy of data_receptor_AI_all_subjects.xlsx")
F_PCA      <- file.path(DATA_DIR, "PCA_results_CBP_all.xlsx")
F_MQS      <- file.path(DATA_DIR, "mqs final.xlsx")
F_MED      <- file.path(DATA_DIR, "combinedpatients medication.xlsx")

## >>> CHANGE 5 <<<  X-AXIS MODE
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
SHOW_HEALTHY_LINE <- FALSE
HEALTHY_LINE_COLOR <- "#3B6D11"

MIN_N_SLOPE <- 5   # groups smaller than this are dropped from slope models

## ---- 1b. FIGURE TEXT SIZING --------------------------------------------------
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
TEXT_SCALE       <- 2
FIG_BASE_SIZE    <- 11   * TEXT_SCALE   # theme base size
FIG_AXIS_TITLE   <- 11   * TEXT_SCALE   # x / y axis titles
FIG_AXIS_TEXT    <- 8.5  * TEXT_SCALE   # tick labels
FIG_STRIP_SIZE   <- 11   * TEXT_SCALE   # facet strip: group name, bold
FIG_STRIP_STATS  <- 8.5  * TEXT_SCALE   # n / r / slope line beneath it, plain
FIG_CAPTION_SIZE <- 8    * TEXT_SCALE   # t-test / slope caption under the plot
FIG_ANNOT_SIZE   <- 2.6  * TEXT_SCALE   # in-panel stats box (geom_label, mm)
FIG_POINT_SIZE   <- 3.0  * TEXT_SCALE   # scatter point size
FIG_LINE_SIZE    <- 0.75 * TEXT_SCALE   # regression / best-fit line width
FIG_CAPTION_COL  <- "grey25"            # caption colour: present, not shouting


## The ONLY two groups that appear in any figure or slope model.
FOCUS_GROUPS <- c("Hydrocodone group", "Tramadol group")

## ---- 2. LOAD + HARMONISE -----------------------------------------------------

RECEPTORS_RAW <- c("HT1A","HT1B","HT2a","HT4","HT6","HTT","D1","D2","DAT","NET",
                   "H3","A4B2","M1","VAChT","CB1","MOR","NMDA","GluR5","GABA")
RECEPTORS <- paste0("R_", c("5HT1A","5HT1B","5HT2A","5HT4","5HT6","5HTT",
                            "D1","D2","DAT","NET","H3","a4b2","M1","VAChT",
                            "CB1","MOR","NMDA","mGluR5","GABA"))

## -- 2a. receptor-related brain activity (this file contains the healthy controls)
receptor <- read_excel(F_RECEPTOR) %>%
  rename(Group_rec = Group, age_rec = Age) %>%
  rename(sex_rec = starts_with("Sex")) %>%
  rename_with(~ RECEPTORS[match(.x, RECEPTORS_RAW)], all_of(RECEPTORS_RAW)) %>%
  rename(Rsqr_adj = `Rsqr-Adj`) %>%
  select(PIN, Group_rec, age_rec, sex_rec, all_of(RECEPTORS), Rsqr_adj)

## -- 2b. behaviour / clinical (CBP patients only)
behavior <- read_excel(F_BEHAVIOR)
names(behavior) <- make.names(names(behavior))
behavior <- behavior %>%
  select(-any_of(c("log.PainDur", "logMME", "logROE", "Log.OpDur", "X"))) %>%
  rename(DOU = Op.Dur, SOWS = sows, COMM = comm) %>%
  mutate(
    Group_beh   = if_else(Group == 0, "CBP-O", "CBP+O"),
    log_PainDur = log10(PainDur),
    log_MME     = log10(MME),
    log_ROE     = log10(ROE),
    log_DOU     = log10(DOU)
  ) %>%
  select(-Group)

## -- 2c. principal components
pca <- read_excel(F_PCA)
names(pca)[1] <- "PIN"
pca <- pca %>% rename(PC1 = `Func Dis`, PC2 = `Pain Sev`, PC3 = `Neg Aff`)

## -- 2d. medication quantification scale
mqs <- read_excel(F_MQS) %>%
  rename(MQS_total = `MQS total`, MQS_nonopioid = `MQS noopioid`)

## -- 2e. medication list -> antidepressant flags
med_raw <- read_excel(F_MED) %>% rename(PIN = ID)
med_cols <- grep("^Med", names(med_raw), value = TRUE)

clean_med <- function(x) {
  x <- as.character(x)
  x <- gsub("\u00a0", " ", x)
  x <- str_squish(tolower(x))
  x[x == "" | x == "nan"] <- NA
  x
}
med_long <- med_raw %>%
  pivot_longer(all_of(med_cols), names_to = "slot", values_to = "med") %>%
  mutate(med = clean_med(med)) %>%
  filter(!is.na(med))

medication <- med_long %>%
  group_by(PIN) %>%
  summarise(
    n_med_classes    = n(),
    antidepressant   = any(str_detect(med, "antidepress")),
    ad_ssri_snri     = any(str_detect(med, "serotonin reuptake")),
    ad_tricyclic     = any(str_detect(med, "tricyclic")),
    ad_other         = any(str_detect(med, "antidepressants .* other")),
    benzodiazepine   = any(str_detect(med, "benzodiazepine")),
    anticonvulsant   = any(str_detect(med, "anticonvulsant")),
    nsaid            = any(str_detect(med, "nsaid")),
    muscle_relaxant  = any(str_detect(med, "muscle relaxant")),
    .groups = "drop"
  )

## ---- 3. MERGE INTO ONE MASTER DATASET ----------------------------------------
master <- receptor %>%
  left_join(behavior,   by = "PIN") %>%
  left_join(pca,        by = "PIN") %>%
  left_join(mqs,        by = "PIN") %>%
  left_join(medication, by = "PIN") %>%
  mutate(
    Group  = factor(Group_rec, levels = c("H", "CBP-O", "CBP+O")),
    age    = coalesce(age, age_rec),
    sex    = coalesce(gender, sex_rec),
    female = sex == 2,
    alcohol_yes = alcohol == 2,
    smoker_yes  = tobacco == 2,
    race_cat = case_when(race == 5 ~ "White",
                         race == 3 ~ "Black",
                         !is.na(race) ~ "Other/Undisclosed",
                         TRUE ~ NA_character_),
    across(c(antidepressant, ad_ssri_snri, ad_tricyclic, ad_other,
             benzodiazepine, anticonvulsant, nsaid, muscle_relaxant),
           ~ if_else(is.na(.x) & Group != "H", FALSE, .x))
  )

## -- 3a. opioid sub-groups (Plot_Group) ---------------------------------------
norm_drug <- function(x) {
  x <- gsub("\u00a0", " ", as.character(x))
  x <- str_squish(tolower(x))
  x[x == ""] <- NA
  x
}
master <- master %>%
  mutate(
    o1 = norm_drug(Opioid.1),
    o2 = norm_drug(Opioid.2),
    o2 = if_else(!is.na(o2) & !is.na(o1) & o2 == o1, NA_character_, o2),
    Plot_Group = case_when(
      Group == "H"                    ~ "Healthy",
      Group == "CBP-O"                ~ "CBP-O",
      !is.na(o2)                      ~ "Multiple opioids",
      is.na(o1)                       ~ "CBP+O unspecified",
      str_detect(o1, "hydrocodone")   ~ "Hydrocodone group",
      str_detect(o1, "tramadol")      ~ "Tramadol group",
      str_detect(o1, "oxycodone")     ~ "Oxycodone only",
      TRUE                            ~ "Other opioid only"
    ),
    Plot_Group = factor(Plot_Group,
                        levels = c("Healthy", "CBP-O", "Hydrocodone group",
                                   "Tramadol group", "Oxycodone only",
                                   "Other opioid only", "Multiple opioids",
                                   "CBP+O unspecified"))
  ) %>%
  droplevels()

## Sanity check: both focus groups must actually exist
present <- intersect(FOCUS_GROUPS, levels(master$Plot_Group))
if (length(present) < 2)
  stop("Both focus groups are required but only found: ",
       paste(present, collapse = ", "))
message(sprintf("Focus groups: %s (n=%d) and %s (n=%d)",
                FOCUS_GROUPS[1], sum(master$Plot_Group == FOCUS_GROUPS[1]),
                FOCUS_GROUPS[2], sum(master$Plot_Group == FOCUS_GROUPS[2])))

## ---- 4. X-AXIS VALUE UTILITIES ----------------------------------------------
min_max_scale <- function(x) (x - min(x, na.rm = TRUE)) /
  (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))

## >>> CHANGE 5 <<<
## Under X_MODE = "raw" this returns the receptor value itself — no healthy
## centroid subtraction, no division by the healthy-to-CBP-O gap. Healthy and
## CBP-O subjects play no role in defining the scale at all.
compute_x_value <- function(data, receptor_subset, scale_mode = X_MODE) {
  M <- as.matrix(data[, receptor_subset, drop = FALSE])
  
  if (scale_mode == "raw") {
    return(as.numeric(rowMeans(M, na.rm = TRUE)))
  }
  
  muH  <- colMeans(M[data$Group == "H", , drop = FALSE], na.rm = TRUE)
  Mc   <- sweep(M, 2, muH, "-")
  v    <- colMeans(Mc[data$Group == "CBP-O", , drop = FALSE], na.rm = TRUE)
  proj <- as.numeric((Mc %*% v) / sum(v^2))
  if (scale_mode == "minmax") proj <- min_max_scale(proj)
  proj
}
## Back-compatible alias so older calls still work
compute_disease_axis <- compute_x_value

## Receptor-set key -> the receptor name(s) as they should READ on a figure:
## every serotonin receptor gets its "5" back ("HT4" -> "5HT4") and the
## structural words become symbols ("HT6_plus_HT4" -> "5HT6 + 5HT4"). Render
## time only -- receptor_sets, filenames and sheet names keep the raw keys.
receptor_set_label <- function(rs_name) {
  special <- c(All_Receptors = "All Receptors",
               Six_Specified = "Six Receptors")
  if (rs_name %in% names(special)) return(unname(special[rs_name]))
  toks <- strsplit(rs_name, "_", fixed = TRUE)[[1]]
  toks <- toks[!toks %in% c("Only", "plus", "and")]
  toks <- sub("^5?HT", "5HT", toks)
  paste(toks, collapse = " + ")
}

## Human-readable x-axis label for a given receptor set
x_axis_label <- function(rs_name) {
  pretty <- receptor_set_label(rs_name)
  if (X_MODE == "raw")
    paste0(pretty, " Activity")
  else
    paste0("Biological Disease Axis (0 = Healthy, 1 = CBP-O) [Set: ", pretty, "]")
}

## Define subsets as requested
receptor_sets <- list(
  "All_Receptors"  = RECEPTORS,
  "MOR_Only"       = c("R_MOR"),
  "HT1A_Only"      = c("R_5HT1A"),
  "HT2A_Only"      = c("R_5HT2A"),
  "HT6_Only"       = c("R_5HT6"),
  "HT4_Only"       = c("R_5HT4"),
  "D1_Only"        = c("R_D1"),
  ## Exploratory only. CB1 gets the same figures, x-axis t-test and centroid
  ## table as every other single receptor here, but is deliberately NOT added to
  ## the mediation (MED_M, PMED_M) or the fixed HT4/HT6 regression models in
  ## sections 15, 17 and 18.
  "CB1_Only"       = c("R_CB1"),
  "HT1A_plus_HT2A" = c("R_5HT1A", "R_5HT2A"),
  "HT6_plus_HT4"   = c("R_5HT6", "R_5HT4"),
  "MOR_plus_D1"    = c("R_MOR", "R_D1"),
  "MOR_HT1A_HT6"   = c("R_MOR", "R_5HT1A", "R_5HT6"),
  "Six_Specified"  = c("R_MOR", "R_5HT1A", "R_5HT2A", "R_5HT6", "R_5HT4", "R_D1")
)

## Default X computation for general use
master$X_disease <- compute_x_value(master, RECEPTORS)

## ---- 5. EXPORT THE MERGED MASTER DATASET -------------------------------------
write_xlsx(master, file.path(OUT_DIR, "master_dataset_merged.xlsx"))
write.csv(master, file.path(OUT_DIR, "master_dataset_merged.csv"), row.names = FALSE)

################################################################################
## 6. TABLE 1 — Hydrocodone vs Tramadol Characteristics (+ CBP-O / CBP+O refs)
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
fmt_p <- function(p) {
  if (is.na(p)) return("--")
  if (p < 1e-5) return("<10\u207b\u2075")
  if (p < 0.001) return(sprintf("%.1e", p))
  sprintf("%.3f", p)
}

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

build_characteristics_table <- function(df, gvar = "Group", groups = c("CBP-O", "CBP+O"), add_p = TRUE) {
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
  hdr <- sprintf("%s (n=%d)", names(ns), as.integer(ns))
  
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

# 1. Build the target table (Hydrocodone vs Tramadol, with p-values)
table1_hydro_tram <- build_characteristics_table(
  master, "Plot_Group", FOCUS_GROUPS
)

# 2. Build the reference column for the full CBP+O group. add_p = FALSE: this
#    column is descriptive only and takes no part in the comparison.
table1_cbpplus <- build_characteristics_table(
  master, "Group", c("CBP+O"), add_p = FALSE
)

# 3. Bind the CBP+O reference column on, then push p-value to the end
table1_hydro_tram <- bind_cols(
  table1_hydro_tram,
  table1_cbpplus[, 2]     # "CBP+O (n=...)"
) %>%
  relocate(`p-value`, .after = last_col())

cat("\n===== TABLE 1: Hydrocodone vs Tramadol (with CBP+O reference) =====\n")
print(as.data.frame(table1_hydro_tram), row.names = FALSE)

write_xlsx(list(`Hydrocodone_vs_Tramadol` = table1_hydro_tram),
           file.path(OUT_DIR, "Table1_Hydrocodone_vs_Tramadol.xlsx"))


################################################################################
## 7. FIGURES AND OVERLAYS  (Hydrocodone group vs Tramadol group)
################################################################################

OUTCOMES <- tribble(
  ~var,  ~label,                                    ~tag,
  "nrs", "Clinical Pain Score (NRS 0-10)",          "NRS",
  "PC1", "PC1 score — Functional disability",       "PC1",
  "PC2", "PC2 score — Pain quality/severity",       "PC2",
  "PC3", "PC3 score — Negative affect",             "PC3"
)

PLOT_GROUPS <- FOCUS_GROUPS

custom_colors <- setNames(c("#1b9e77", "#d95f02")[seq_along(PLOT_GROUPS)], PLOT_GROUPS)
## Both groups are filled circles. Shape therefore carries no information and
## its legend is switched off wherever it would otherwise be drawn; colour alone
## separates the groups.
custom_shapes <- setNames(rep(16, length(PLOT_GROUPS)), PLOT_GROUPS)

## ---- 7a. WHERE THE n / r / slope / p TEXT GOES --------------------------------
## Two placements, both produced on every run into parallel folders that hold
## identical filenames, so the same figure can be compared side by side:
##
##   "header" -> the stats become a second line on the strip above each panel
##               (and the legend key text on the overlay figures). The text is
##               outside the panel, so it cannot touch the data at all.
##   "corner" -> the stats stay inside the panel, in a box placed in whichever
##               corner is emptiest (see pick_corner below).
STATS_MODES <- c(header = "figs_stats_header", corner = "figs_stats_corner")
for (.d in file.path(OUT_DIR, STATS_MODES)) dir.create(.d, showWarnings = FALSE)

## Path for a figure in one of the two placement folders.
fig_path <- function(stats_mode, file_stub)
  file.path(OUT_DIR, STATS_MODES[[stats_mode]], paste0(file_stub, ".png"))

## Saved figure sizes, in inches. Named because the corner placement has to know
## the size the figure will be SAVED at to work out how much of a panel the stats
## box covers -- so the number ggsave() gets and the number the placement is
## measured against must be the same one.
FIG1_W <- 13; FIG1_H <- 8      # section 7/9  faceted
OVL_W  <- 11; OVL_H  <- 7      # section 7/9  overlay
XY_W   <- 12; XY_H   <- 7.4    # section 14   faceted
XYO_W  <- 11; XYO_H  <- 6.8    # section 14   overlay

## How much of a panel the stats box occupies, as a fraction of the panel's
## width and height. These are only the FALLBACK values, used when the box and
## panel cannot be measured; the real numbers are measured per figure by
## box_fracs() below, because guessing them is exactly how a box ends up sitting
## on a data point.
LAB_W_FRAC <- 0.50
LAB_H_FRAC <- 0.30

## Physical size of one stats box, in inches. ggplot expresses text `size` in mm
## of font height, hence the .pt conversion to points; label.padding is in lines
## of that font and applies on both sides.
label_size_in <- function(txt, size_mm = FIG_ANNOT_SIZE, lineheight = 1.05,
                          pad_lines = 0.35) {
  fs <- size_mm * ggplot2::.pt
  g  <- grid::textGrob(txt, gp = grid::gpar(fontsize = fs, fontface = "bold",
                                            lineheight = lineheight))
  pad <- 2 * pad_lines * fs / 72.27
  c(w = grid::convertWidth(grid::grobWidth(g),   "in", valueOnly = TRUE) + pad,
    h = grid::convertHeight(grid::grobHeight(g), "in", valueOnly = TRUE) + pad)
}

## Size of ONE panel, in inches, for a plot that will be saved at fig_w x fig_h.
## In a gtable the panels are the rows/columns measured in "null" units and
## everything else (axes, strips, caption, margins) has a fixed size, so the
## panels get whatever is left over after the fixed parts are paid for.
panel_size_in <- function(p, fig_w, fig_h) {
  gt <- ggplot2::ggplotGrob(p)
  hn <- grid::unitType(gt$heights) == "null"
  wn <- grid::unitType(gt$widths)  == "null"
  fixed_h <- if (any(!hn)) sum(grid::convertHeight(gt$heights[!hn], "in", TRUE)) else 0
  fixed_w <- if (any(!wn)) sum(grid::convertWidth(gt$widths[!wn],   "in", TRUE)) else 0
  c(w = max((fig_w - fixed_w) / max(1L, sum(wn)), 0.1),
    h = max((fig_h - fixed_h) / max(1L, sum(hn)), 0.1))
}

## Box size as a fraction of the panel, measured. `p` must already carry its
## facet spec (the panel count sets the panel width) but not the label itself.
## Falls back to the constants above if anything about the measurement fails --
## on a device-less or unusual setup a wrong answer is worse than a rough one.
box_fracs <- function(p, fig_w, fig_h, txt) {
  out <- try({
    pm <- panel_size_in(p, fig_w, fig_h)
    ## widest/tallest label wins: one size is used for every panel
    ls <- vapply(txt, label_size_in, numeric(2))
    c(w = max(ls["w", ]) / pm[["w"]], h = max(ls["h", ]) / pm[["h"]])
  }, silent = TRUE)
  if (inherits(out, "try-error") || !all(is.finite(out)))
    return(c(w = LAB_W_FRAC, h = LAB_H_FRAC))
  ## Never claim the box is bigger than the panel, and keep a little slack.
  c(w = min(out[["w"]] * 1.05, 0.9), h = min(out[["h"]] * 1.05, 0.9))
}

## Refit the same model geom_smooth draws and return the fitted line with its
## 95% confidence band over an x grid. Needed because that ribbon regularly
## extends ABOVE the highest data point and BELOW the lowest, so max(y)/min(y)
## are not the visual extremes of the panel.
fit_band <- function(d, xv, yv, n_grid = 100) {
  ## Copied into plain .x / .y columns so the formula never has to survive
  ## backtick-quoting of names like "log_ROE" or "R_5HT4".
  dd <- data.frame(.x = as.numeric(d[[xv]]), .y = as.numeric(d[[yv]]))
  dd <- dd[is.finite(dd$.x) & is.finite(dd$.y), , drop = FALSE]
  if (nrow(dd) < 3 || length(unique(dd$.x)) < 2) return(NULL)
  g  <- data.frame(.x = seq(min(dd$.x), max(dd$.x), length.out = n_grid))
  pr <- try(suppressWarnings(
    predict(lm(.y ~ .x, data = dd), newdata = g, interval = "confidence")),
    silent = TRUE)
  if (inherits(pr, "try-error")) return(NULL)
  out <- cbind(g, as.data.frame(pr))
  names(out)[1] <- xv
  out
}

## The full vertical extent of what will be drawn in one panel: data points and
## the confidence ribbon, whichever reaches further.
panel_span <- function(d, xv, yv) {
  fb <- fit_band(d, xv, yv)
  ys <- c(d[[yv]], if (!is.null(fb)) c(fb$lwr, fb$upr))
  ys <- ys[is.finite(ys)]
  if (!length(ys)) return(c(NA_real_, NA_real_))
  range(ys)
}

## Score all four corners of one panel and return the emptiest.
##
## Penalty = number of that group's points falling inside the candidate box
##           + number of ribbon grid steps whose [lwr, upr] span crosses it
##           + a charge for sitting on the healthy-mean guide (see below).
## Lowest penalty wins; ties break topleft -> topright -> bottomleft ->
## bottomright, which is the order a reader's eye goes.
##
## `xr` / `yr` are the panel limits, which are shared across facets
## (scales = "fixed"), so a box sized as a fraction of them is the same physical
## size in every panel.
##
## `ref_x` is where healthy_line() draws its dotted vertical line, which carries
## a rotated "Healthy mean" label near the TOP of the panel. Hiding that label
## costs the reader real information, so it is charged like several data points;
## crossing the bare line lower down is only cosmetic and charged as one.
pick_corner <- function(d, xv, yv, xr, yr, ref_x = NA_real_,
                        w_frac = LAB_W_FRAC, h_frac = LAB_H_FRAC) {
  w <- w_frac * diff(xr)
  h <- h_frac * diff(yr)
  fb <- fit_band(d, xv, yv)

  ## anchor point, plus the box's extent measured from it
  cand <- tibble(
    corner = c("topleft", "topright", "bottomleft", "bottomright"),
    lab_x  = c(xr[1], xr[2], xr[1], xr[2]),
    lab_y  = c(yr[2], yr[2], yr[1], yr[1]),
    hjust  = c(0, 1, 0, 1),
    vjust  = c(1, 1, 0, 0)
  )

  cand$penalty <- vapply(seq_len(nrow(cand)), function(i) {
    x_lo <- if (cand$hjust[i] == 0) cand$lab_x[i] else cand$lab_x[i] - w
    x_hi <- x_lo + w
    y_lo <- if (cand$vjust[i] == 0) cand$lab_y[i] else cand$lab_y[i] - h
    y_hi <- y_lo + h

    n_pts <- sum(d[[xv]] >= x_lo & d[[xv]] <= x_hi &
                 d[[yv]] >= y_lo & d[[yv]] <= y_hi, na.rm = TRUE)

    n_rib <- if (is.null(fb)) 0L else
      sum(fb[[xv]] >= x_lo & fb[[xv]] <= x_hi &
          fb$upr >= y_lo & fb$lwr <= y_hi, na.rm = TRUE)

    n_ref <- if (is.na(ref_x) || ref_x < x_lo || ref_x > x_hi) 0 else
      if (cand$vjust[i] == 1) 5 else 1

    as.numeric(n_pts + n_rib + n_ref)
  }, numeric(1))

  best <- which.min(cand$penalty)          # which.min takes the first on a tie
  list(lab_x   = cand$lab_x[best],
       lab_y   = cand$lab_y[best],
       hjust   = cand$hjust[best],
       vjust   = cand$vjust[best],
       corner  = cand$corner[best],
       penalty = cand$penalty[best],
       ## Which side the box sits on, for the headroom fallback below.
       side    = if (cand$vjust[best] == 1) "top" else "bottom")
}

## Place the stats box for every group in one figure, and work out the panel
## limits it needs.
##
## On a crowded panel every corner is occupied, and the emptiest one still
## covers something. When that happens the panel limit on the chosen side is
## pushed out by one box height, which turns that corner into empty space and
## guarantees nothing is hidden. Facets share scales, so the expansion is the
## largest one any group asked for and applies to both panels. Panels where a
## corner was genuinely clear are left at their natural limits, so no figure
## carries dead space it did not need.
## `gvar` groups the PLACEMENT: one box per facet panel, or a single pooled box
## on an overlay (where callers pass a constant column). `extent_var` groups the
## EXTENT calculation and must always name the real drawing groups, because an
## overlay draws one regression ribbon per group and a pooled fit through both
## would not reach as far as either of them.
place_corner_labels <- function(plot_data, gvar, xv, yv, ref_x = NA_real_,
                                extent_var = gvar,
                                w_frac = LAB_W_FRAC, h_frac = LAB_H_FRAC) {
  xr <- range(plot_data[[xv]], na.rm = TRUE)

  spans <- lapply(split(plot_data, droplevels(plot_data[[extent_var]])),
                  panel_span, xv, yv)
  yr <- range(unlist(spans), na.rm = TRUE)
  if (!all(is.finite(yr)) || diff(yr) == 0) yr <- range(plot_data[[yv]], na.rm = TRUE)
  if (diff(yr) == 0) yr <- yr + c(-0.5, 0.5)

  ## A reference line outside the plotted x-range is never in the way.
  if (!is.na(ref_x) && (ref_x < xr[1] || ref_x > xr[2])) ref_x <- NA_real_

  groups <- split(plot_data, droplevels(plot_data[[gvar]]))

  ## The box is a fixed PHYSICAL size, so it always covers h_frac of the panel's
  ## height -- which means widening the y-limits shrinks the box in data units,
  ## and pick_corner's verdict depends on the limits it was given. So: pick,
  ## expand, then pick again against the expanded limits, and only stop when
  ## every box is clear.
  ##
  ## Band height, for `sides` padded sides. The box must fit inside the band
  ## with room to spare, so solve
  ##   pad = SLACK * h_frac * (diff(yr) + sides * pad)
  ##     => pad = SLACK * h_frac * diff(yr) / (1 - SLACK * sides * h_frac)
  ## Two details matter. Padding widens the y-range, and the box is a fixed
  ## physical size, so its height in data units grows with the range -- hence
  ## `sides * pad` on the right. And SLACK > 1 keeps the box off the band's own
  ## edge: solved exactly, the box's lower edge lands precisely on the highest
  ## data point, which still counts as covering it.
  SLACK <- 1.15
  band <- function(span, sides) {
    denom <- 1 - SLACK * sides * h_frac
    if (denom <= 0.05) SLACK * h_frac * span * (sides + 1)
    else SLACK * h_frac * span / denom
  }

  ylim    <- yr
  pad_top <- 0
  pad_bot <- 0

  for (iter in 1:4) {
    pos  <- lapply(groups, pick_corner, xv, yv, xr, ylim, ref_x, w_frac, h_frac)
    need_top <- any(vapply(pos, function(p) p$penalty > 0 && p$side == "top",
                           logical(1)))
    need_bot <- any(vapply(pos, function(p) p$penalty > 0 && p$side == "bottom",
                           logical(1)))
    if (!need_top && !need_bot) break

    pad_top <- max(pad_top, if (need_top) 1 else 0)
    pad_bot <- max(pad_bot, if (need_bot) 1 else 0)
    pad     <- band(diff(yr), pad_top + pad_bot)
    ylim    <- c(yr[1] - pad_bot * pad, yr[2] + pad_top * pad)
  }
  pad_top <- ylim[2] - yr[2]
  pad_bot <- yr[1] - ylim[1]

  ## Re-anchor to the expanded limits so a box on a padded side lands in the new
  ## blank band rather than back on top of the data.
  coords <- tibble(
    .group = factor(names(pos), levels = levels(droplevels(plot_data[[gvar]]))),
    lab_x  = vapply(pos, `[[`, numeric(1), "lab_x"),
    lab_y  = ifelse(vapply(pos, `[[`, numeric(1), "vjust") == 1, ylim[2], ylim[1]),
    hjust  = vapply(pos, `[[`, numeric(1), "hjust"),
    vjust  = vapply(pos, `[[`, numeric(1), "vjust"),
    corner = vapply(pos, `[[`, character(1), "corner")
  )
  names(coords)[1] <- gvar

  ## Residual check, against the limits the figure will actually be drawn with:
  ## how much still falls inside each box. This is the number that matters, so it
  ## is computed here and audited rather than assumed to be zero.
  resid <- do.call(rbind, lapply(seq_len(nrow(coords)), function(k) {
    g  <- as.character(coords[[gvar]][k])
    dg <- groups[[g]]
    w  <- w_frac * diff(xr)
    h  <- h_frac * diff(ylim)
    x_lo <- if (coords$hjust[k] == 0) coords$lab_x[k] else coords$lab_x[k] - w
    y_lo <- if (coords$vjust[k] == 0) coords$lab_y[k] else coords$lab_y[k] - h
    x_hi <- x_lo + w; y_hi <- y_lo + h
    fb <- fit_band(dg, xv, yv)
    data.frame(
      panel  = g,
      corner = coords$corner[k],
      pts_covered = sum(dg[[xv]] >= x_lo & dg[[xv]] <= x_hi &
                        dg[[yv]] >= y_lo & dg[[yv]] <= y_hi, na.rm = TRUE),
      rib_covered = if (is.null(fb)) 0L else
        sum(fb[[xv]] >= x_lo & fb[[xv]] <= x_hi &
            fb$upr >= y_lo & fb$lwr <= y_hi, na.rm = TRUE),
      pad_top = pad_top, pad_bot = pad_bot,
      w_frac = w_frac, h_frac = h_frac, stringsAsFactors = FALSE)
  }))

  list(coords = coords, ylim = ylim, xlim = xr, yr = yr,
       w_frac = w_frac, h_frac = h_frac, residual = resid)
}

## Audit trail for the corner placement. Every corner-mode figure appends its
## residual check here; section 14 then reports the total and writes it out, so a
## box that covers something shows up as a loud warning instead of being found
## later by eye in a composite.
CORNER_AUDIT <- list()
audit_corner <- function(figure, lp) {
  if (is.null(lp$residual)) return(invisible(NULL))
  CORNER_AUDIT[[length(CORNER_AUDIT) + 1]] <<-
    cbind(figure = figure, lp$residual, stringsAsFactors = FALSE)
  invisible(NULL)
}

## Width of one line of text, in inches, at a given point size. Used to decide
## whether the stats fit on a single strip line before committing to it.
text_width_in <- function(txt, size_pt, face = "plain") {
  g <- grid::textGrob(txt, gp = grid::gpar(fontsize = size_pt, fontface = face))
  grid::convertWidth(grid::grobWidth(g), "in", valueOnly = TRUE)
}

## Facet strip carrying a bold group name above a smaller, plain stats line.
##
## ggplot styles a strip with ONE element_text, so two sizes/weights in one strip
## normally needs ggtext. plotmath does it without the dependency: atop() stacks
## the lines, bold() weights the first, and scriptstyle() shrinks the second to
## about 0.8x. The catch is that a lookup vector passed to labeller() is never
## parsed, so the panels are faceted on a column whose LEVELS are the plotmath
## strings and label_parsed is applied to those.
## `stats_txt` may contain newlines; each becomes another stacked atop() line.
strip_math <- function(name, stats_txt) {
  ## Quotes would end the plotmath string early; nothing else needs escaping.
  clean <- function(s) gsub('"', "", s, fixed = TRUE)
  vapply(seq_along(name), function(i) {
    lines <- strsplit(clean(stats_txt[i]), "\n", fixed = TRUE)[[1]]
    inner <- sprintf('scriptstyle(plain("%s"))', lines)
    body  <- inner[length(inner)]
    for (k in rev(seq_len(length(inner) - 1)))
      body <- sprintf("atop(%s, %s)", inner[k], body)
    sprintf('atop(bold("%s"), %s)', clean(name[i]), body)
  }, character(1))
}

## Add the plotmath strip column to `plot_data`, ready for
## facet_wrap(~ .strip, labeller = label_parsed).
strip_facet_data <- function(plot_data, gvar, names_in_order, stats_txt) {
  labs_v <- setNames(strip_math(names_in_order, stats_txt), names_in_order)
  plot_data$.strip <- factor(labs_v[as.character(plot_data[[gvar]])],
                             levels = unname(labs_v))
  plot_data
}

## Width one facet strip has to play with, in inches: the figure minus the
## y-axis furniture (title, tick labels, margins), split across the panels.
## Deliberately an estimate -- the exact number needs a built plot, and the only
## thing it decides is whether the stats take one strip line or two, with two
## always safe. Doing it this way keeps the strip text settled BEFORE the plot is
## built, so the data never has to be swapped in afterwards.
##
## Checked against the real thing: this returns 5.55 / 5.05 in for the two
## faceted figure sizes, where ggplotGrob measures 5.67 / 5.17 -- i.e. it
## under-reads by about 0.12 in, which is the margin the 0.97 below allows for.
strip_budget_in <- function(fig_w, n_panels) max((fig_w - 1.9) / n_panels, 1)

## TRUE if every stats line fits on one strip line. scriptstyle() renders at
## about 0.8x, which is the size these are measured at.
stats_fit_one_line <- function(txt, fig_w, n_panels) {
  all(vapply(txt, text_width_in, numeric(1), FIG_STRIP_SIZE * 0.8) <
        0.97 * strip_budget_in(fig_w, n_panels))
}

## Shared theme for every scatter figure. Pulling it out of the four builders is
## what keeps the sizes consistent between them -- they had drifted apart.
theme_fig <- function(faceted = TRUE, legend = "none") {
  th <- theme_minimal(base_size = FIG_BASE_SIZE) +
    theme(
      axis.title      = element_text(size = FIG_AXIS_TITLE),
      axis.text       = element_text(size = FIG_AXIS_TEXT, colour = "grey30"),
      legend.position = legend,
      legend.text     = element_text(size = FIG_STRIP_STATS),
      plot.caption    = element_text(hjust = 0, size = FIG_CAPTION_SIZE,
                                     colour = FIG_CAPTION_COL, face = "plain",
                                     lineheight = 1.25, margin = margin(t = 8)),
      plot.margin     = margin(6, 10, 4, 6)
    )
  if (faceted)
    th <- th + theme(
      ## Bold, because this strip carries the group name. Header-mode figures
      ## override it to "plain" -- there the label is plotmath and supplies its
      ## own bold() for the name and scriptstyle(plain()) for the stats line.
      strip.text   = element_text(size = FIG_STRIP_SIZE, face = "bold",
                                  lineheight = 1.05,
                                  margin = margin(b = 6, t = 2)),
      panel.border = element_rect(colour = "grey80", fill = NA, linewidth = 0.8)
    )
  th
}

## One line, for a facet strip and for the overlay legend keys. Single-spaced
## after the commas: at two double spaces it measured 4.92 in against a 5.17 in
## section-14 panel, which is too close to rely on.
stats_line <- function(n, r, slope, p)
  sprintf("n = %d, r = %.2f, slope = %.2f (p = %s)",
          n, r, slope, vapply(p, fmt_p, character(1)))

## Two lines, for the strip above a faceted panel. A facet strip is only about
## half the figure wide and ggplot does not wrap strip text -- it just clips it
## at the panel edge -- so the one-line form above does not fit here.
stats_strip <- function(n, r, slope, p)
  sprintf("n = %d,  r = %.2f\nslope = %.2f (p = %s)",
          n, r, slope, vapply(p, fmt_p, character(1)))

## Three lines, for the box inside a panel in the "corner" placement.
stats_block <- function(n, r, slope, p)
  sprintf("n = %d\nr = %.2f\nslope = %.2f (p = %s)",
          n, r, slope, vapply(p, fmt_p, character(1)))

## A rotated y-axis title is limited by the panel HEIGHT, and ggplot clips it
## rather than wrapping. The outcome labels are long enough ("Clinical Pain
## Score (NRS 0-10)", "PC1 score — Functional disability") that at
## FIG_BASE_SIZE they overrun a panel of the height these figures use, so they
## are wrapped here instead of being shortened -- the full wording still shows.
wrap_lab <- function(s, width = 18)
  vapply(s, function(x) {
    if (nchar(x) <= width) return(x)
    ## Break at the em-dash when there is one: these labels are "<name> —
    ## <description>", and splitting there keeps both halves whole. Plain
    ## strwrap would cut "PC3 score — Negative affect" into
    ## "PC3 score — Negative" / "affect".
    if (grepl(" — ", x, fixed = TRUE)) {
      parts <- strsplit(x, " — ", fixed = TRUE)[[1]]
      return(paste0(parts[1], " —\n",
                    paste(parts[-1], collapse = " — ")))
    }
    paste(strwrap(x, width = width), collapse = "\n")
  }, character(1), USE.NAMES = FALSE)

## Welch two-sample t-test on the X-axis values of the two groups.
x_axis_ttest <- function(d,
                         xvar    = "X_disease",
                         gvar    = "Plot_Group",
                         g1      = FOCUS_GROUPS[1],
                         g2      = FOCUS_GROUPS[2],
                         outcome = NA_character_,
                         rs_name = NA_character_) {
  
  x1 <- d[[xvar]][d[[gvar]] == g1 & !is.na(d[[xvar]])]
  x2 <- d[[xvar]][d[[gvar]] == g2 & !is.na(d[[xvar]])]
  
  blank_row <- tibble(receptor_set = rs_name, outcome = outcome,
                      group_1 = g1, group_2 = g2,
                      n_1 = length(x1), n_2 = length(x2),
                      mean_X_1 = if (length(x1)) mean(x1) else NA_real_,
                      mean_X_2 = if (length(x2)) mean(x2) else NA_real_,
                      sd_X_1 = if (length(x1) > 1) sd(x1) else NA_real_,
                      sd_X_2 = if (length(x2) > 1) sd(x2) else NA_real_,
                      mean_diff = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
                      t = NA_real_, df = NA_real_, p_value = NA_real_,
                      cohens_d = NA_real_, significant = NA)
  
  if (length(x1) < 2 || length(x2) < 2) {
    return(list(label = "X-axis t-test: not enough data in one group",
                row = blank_row))
  }
  
  tt <- stats::t.test(x1, x2)                       # Welch by default
  s_pooled <- sqrt(((length(x1) - 1) * var(x1) + (length(x2) - 1) * var(x2)) /
                     (length(x1) + length(x2) - 2))
  d_eff <- (mean(x1) - mean(x2)) / s_pooled
  
  row <- blank_row %>%
    mutate(mean_diff   = unname(diff(rev(tt$estimate))),
           ci_lo       = tt$conf.int[1],
           ci_hi       = tt$conf.int[2],
           t           = unname(tt$statistic),
           df          = unname(tt$parameter),
           p_value     = tt$p.value,
           cohens_d    = d_eff,
           significant = tt$p.value < 0.05)
  
  label_raw <- sprintf(
    "X-axis t-test — %s: mean = %.4f \u00b1 %.4f (n=%d) vs %s: mean = %.4f \u00b1 %.4f (n=%d). Welch t(%.1f) = %.2f, p = %s, Cohen's d = %.2f %s",
    g1, mean(x1), sd(x1), length(x1),
    g2, mean(x2), sd(x2), length(x2),
    unname(tt$parameter), unname(tt$statistic), fmt_p(tt$p.value), d_eff,
    if (tt$p.value < 0.001) "***" else
      if (tt$p.value < 0.01) "**" else
        if (tt$p.value < 0.05) "*" else "(n.s.)"
  )
  label <- paste(strwrap(label_raw, width = 115), collapse = "\n")
  
  list(label = label, row = row)
}

## >>> CHANGE 5 <<<
## Healthy reference line ONLY. The CBP-O line at X = 1 is gone, since under
## raw mode there is no X = 1. `ref_x` is the healthy group mean of whatever
## is on the x-axis (≈ 0 for raw receptor values, exactly 0 for projection).
## The line is drawn before the points, so anything it carries ends up UNDER the
## data -- the rotated "Healthy mean" text it used to draw sat on top of the
## points in one panel and was itself obscured in others. The line is now
## unlabelled and the caption says what it is instead (see healthy_line_note),
## which keeps the panel free of text that can collide with data. Pass
## label = TRUE to get the old in-panel text back.
healthy_line <- function(ref_x, label = FALSE) {
  if (!SHOW_HEALTHY_LINE || is.na(ref_x)) return(NULL)
  out <- list(
    geom_vline(xintercept = ref_x, color = HEALTHY_LINE_COLOR,
               linetype = "dotted", linewidth = 0.9)
  )
  if (label) {
    out <- c(out, list(
      annotate("text", x = ref_x, y = Inf, label = "Healthy mean",
               color = HEALTHY_LINE_COLOR, angle = 90, hjust = 1.15, vjust = -0.4,
               size = FIG_ANNOT_SIZE * 0.8, fontface = "bold")
    ))
  }
  out
}

## Caption sentence replacing that in-panel label.
healthy_line_note <- function(ref_x) {
  if (!SHOW_HEALTHY_LINE || is.na(ref_x)) return("")
  sprintf("\nGreen dotted line = healthy group mean (x = %.4f).", ref_x)
}

## Per-group n, r, slope and slope p-value for one x/y pair. Shared by the
## faceted and overlay builders so both report the same numbers.
group_fit_summary <- function(plot_data, yvar, xvar = "X_disease") {
  plot_data %>%
    group_by(Plot_Group) %>%
    group_modify(~ {
      if (nrow(.x) < 3)
        return(tibble(n = nrow(.x), mean_x = NA_real_, mean_y = NA_real_,
                      slope = NA_real_, p = NA_real_, r = NA_real_))
      fit <- lm(.x[[yvar]] ~ .x[[xvar]])
      tibble(n      = nrow(.x),
             mean_x = mean(.x[[xvar]]),
             mean_y = mean(.x[[yvar]]),
             slope  = coef(fit)[2],
             p      = summary(fit)$coefficients[2, 4],
             r      = suppressWarnings(cor(.x[[xvar]], .x[[yvar]])))
    }) %>% ungroup()
}

## No plot title and no subtitle. The axis labels and the strip labels already
## say which receptor set, which outcome and which group a panel shows, and
## dropping them buys back ~1.5 in of panel height per figure -- which is what
## makes the panels legible once several are composed onto one page.
make_figure1 <- function(data, yvar, ylab, tag,
                         receptor_subset = RECEPTORS,
                         file_stub = NULL,
                         rs_name = "All_Receptors",
                         ref_x = NA_real_,
                         stats_mode = "corner") {

  stopifnot(stats_mode %in% names(STATS_MODES))
  
  plot_data <- data %>%
    filter(Plot_Group %in% PLOT_GROUPS, !is.na(.data[[yvar]]), !is.na(X_disease)) %>%
    mutate(Plot_Group = droplevels(factor(Plot_Group, levels = PLOT_GROUPS)))
  
  tt <- x_axis_ttest(plot_data, outcome = tag, rs_name = rs_name)
  
  summary_data <- group_fit_summary(plot_data, yvar) %>%
    mutate(label_text = stats_block(n, r, slope, p),
           stat_line  = stats_line(n, r, slope, p),
           stat_strip = stats_strip(n, r, slope, p))

  ## Strip labels are settled here, before the plot is built, so the header
  ## branch only has to add the facet. One line if the stats fit the strip's
  ## width, two if not -- ggplot clips strip text at the panel edge rather than
  ## wrapping it.
  plot_data <- strip_facet_data(
    plot_data, "Plot_Group", as.character(summary_data$Plot_Group),
    if (stats_fit_one_line(summary_data$stat_line, FIG1_W, 2))
      summary_data$stat_line else summary_data$stat_strip)

  p <- ggplot(plot_data, aes(x = X_disease, y = .data[[yvar]],
                             color = Plot_Group, shape = Plot_Group)) +
    healthy_line(ref_x) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "grey30",
                linetype = "solid", linewidth = FIG_LINE_SIZE,
                fill = "grey80", alpha = 0.4) +
    geom_point(size = FIG_POINT_SIZE, alpha = 0.85) +
    scale_color_manual(values = custom_colors) +
    scale_shape_manual(values = custom_shapes) +
    labs(x = wrap_lab(x_axis_label(rs_name), 48), y = wrap_lab(ylab, 22),
         caption = paste0(tt$label, healthy_line_note(ref_x))) +
    theme_fig(faceted = TRUE)

  if (stats_mode == "header") {
    ## Stats live on the strip, so the panel holds nothing but data. The label is
    ## plotmath and carries its own weights, so the strip element must be plain.
    p <- p + facet_wrap(~ .strip, scales = "fixed", labeller = label_parsed) +
      theme(strip.text = element_text(face = "plain", size = FIG_STRIP_SIZE,
                                      lineheight = 1.05,
                                      margin = margin(b = 6, t = 2)))
  } else {
    ## Measure the box against the faceted panel it will actually be drawn in,
    ## then place it. ref_x is passed so it does not land on the "Healthy mean"
    ## guide either.
    bf <- box_fracs(p + facet_wrap(~ Plot_Group, scales = "fixed"),
                    FIG1_W, FIG1_H, summary_data$label_text)
    lp <- place_corner_labels(plot_data, "Plot_Group", "X_disease", yvar,
                              ref_x = if (SHOW_HEALTHY_LINE) ref_x else NA_real_,
                              w_frac = bf[["w"]], h_frac = bf[["h"]])
    audit_corner(if (is.null(file_stub)) paste0("Figure1_", tag, "_", rs_name)
                 else file_stub, lp)
    lab_df <- left_join(lp$coords, summary_data, by = "Plot_Group")
    p <- p +
      facet_wrap(~ Plot_Group, scales = "fixed") +
      geom_label(data = lab_df, aes(x = lab_x, y = lab_y, label = label_text),
                 inherit.aes = FALSE,
                 hjust = lab_df$hjust, vjust = lab_df$vjust,
                 size = FIG_ANNOT_SIZE, lineheight = 1.05,
                 fontface = "bold", color = "black", fill = "white", alpha = 1,
                 linewidth = 0, label.padding = unit(0.35, "lines")) +
      coord_cartesian(ylim = lp$ylim)
  }

  if (is.null(file_stub)) file_stub <- paste0("Figure1_", tag, "_", rs_name)
  ggsave(fig_path(stats_mode, file_stub), p, width = FIG1_W, height = FIG1_H,
         dpi = 300)

  list(plot = p, summary = summary_data, data = plot_data, ttest = tt$row)
}

make_overlay <- function(data, yvar, ylab, tag, file_stub = NULL,
                         rs_name = "All_Receptors", ref_x = NA_real_,
                         stats_mode = "corner") {

  stopifnot(stats_mode %in% names(STATS_MODES))

  plot_data <- data %>%
    filter(Plot_Group %in% PLOT_GROUPS, !is.na(.data[[yvar]]), !is.na(X_disease)) %>%
    mutate(Plot_Group = droplevels(factor(Plot_Group, levels = PLOT_GROUPS)))

  tt <- x_axis_ttest(plot_data, outcome = tag, rs_name = rs_name)

  summary_data <- group_fit_summary(plot_data, yvar)

  p <- ggplot(plot_data, aes(X_disease, .data[[yvar]],
                             color = Plot_Group, shape = Plot_Group)) +
    healthy_line(ref_x) +
    geom_point(alpha = 0.6, size = FIG_POINT_SIZE) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linetype = "solid",
                linewidth = FIG_LINE_SIZE) +
    ## Both groups are circles now, so a shape legend would only repeat the
    ## colour legend.
    scale_shape_manual(values = custom_shapes, guide = "none") +
    labs(x = wrap_lab(x_axis_label(rs_name), 48), y = wrap_lab(ylab, 22),
         color = NULL, shape = NULL,
         caption = paste0(tt$label, healthy_line_note(ref_x))) +
    theme_fig(faceted = FALSE, legend = "bottom")

  if (stats_mode == "header") {
    ## An overlay has no strip, so the legend key text carries the stats. Still
    ## outside the panel, so still incapable of covering a point.
    key_labs <- setNames(
      paste0(as.character(summary_data$Plot_Group), ": ",
             stats_line(summary_data$n, summary_data$r,
                        summary_data$slope, summary_data$p)),
      as.character(summary_data$Plot_Group))
    p <- p +
      scale_color_manual(values = custom_colors, labels = key_labs) +
      guides(color = guide_legend(ncol = 1)) +
      ## One key per row, left-aligned and plain: a centred key this long
      ## overruns the device and ggplot clips it rather than wrapping.
      theme(legend.justification = "left",
            legend.text          = element_text(size = FIG_STRIP_STATS,
                                               face = "plain"),
            legend.margin        = margin(0, 0, 0, 0),
            legend.box.margin    = margin(0, 0, 0, 0),
            legend.key.spacing.y = unit(2, "pt"))
  } else {
    ## One shared box for both groups. The corner search runs on the pooled data
    ## because every point in the figure is drawn in this single panel, so
    ## "emptiest" has to be judged against all of them at once.
    box <- paste(sprintf("%s: %s", as.character(summary_data$Plot_Group),
                         stats_line(summary_data$n, summary_data$r,
                                    summary_data$slope, summary_data$p)),
                 collapse = "\n")
    bf <- box_fracs(p, OVL_W, OVL_H, box)
    lp <- place_corner_labels(
      plot_data %>% mutate(.all = factor("all")), ".all", "X_disease", yvar,
      ref_x = if (SHOW_HEALTHY_LINE) ref_x else NA_real_,
      extent_var = "Plot_Group", w_frac = bf[["w"]], h_frac = bf[["h"]])
    audit_corner(paste0(if (is.null(file_stub))
                          paste0("Figure1_", tag, "_overlay_", rs_name)
                        else file_stub), lp)
    p <- p +
      scale_color_manual(values = custom_colors) +
      annotate("label", x = lp$coords$lab_x[1], y = lp$coords$lab_y[1],
               hjust = lp$coords$hjust[1], vjust = lp$coords$vjust[1],
               label = box, size = FIG_ANNOT_SIZE, lineheight = 1.05,
               fontface = "bold", color = "black", fill = "white",
               linewidth = 0, label.padding = unit(0.35, "lines")) +
      coord_cartesian(ylim = lp$ylim)
  }

  if (is.null(file_stub)) file_stub <- paste0("Figure1_", tag, "_overlay_", rs_name)
  ggsave(fig_path(stats_mode, file_stub), p, width = OVL_W, height = OVL_H,
         dpi = 300)
  p
}

################################################################################
## 8. STATS FUNCTION
################################################################################

slope_stats <- function(data, yvar, xvar = "X_disease", gvar = "Plot_Group", min_n = MIN_N_SLOPE) {
  d <- data %>% filter(.data[[gvar]] %in% PLOT_GROUPS, !is.na(.data[[yvar]]), !is.na(.data[[xvar]])) %>%
    mutate(.g = droplevels(factor(.data[[gvar]], levels = PLOT_GROUPS)))
  keep <- names(which(table(d$.g) >= min_n))
  d <- d %>% filter(.g %in% keep) %>% mutate(.g = droplevels(.g))
  lev <- levels(d$.g)
  if (length(lev) < 2) return(NULL)
  
  f_add <- as.formula(paste(yvar, "~", xvar, "+ .g"))
  f_int <- as.formula(paste(yvar, "~", xvar, "* .g"))
  m_add <- lm(f_add, data = d)
  m_int <- lm(f_int, data = d)
  av <- anova(m_add, m_int)
  
  global <- tibble(outcome = yvar, test = "Slope homogeneity (X x Group interaction)",
                   df1 = av$Df[2], df2 = av$Res.Df[2], F_value = av$F[2],
                   p_value = av$`Pr(>F)`[2], adj_R2_full = summary(m_int)$adj.r.squared)
  
  cf <- coef(m_int); V <- vcov(m_int); nm <- names(cf); dfres <- m_int$df.residual
  svec <- function(g) {
    v <- setNames(numeric(length(cf)), nm)
    v[xvar] <- 1
    if (g != lev[1]) {
      trm <- paste0(xvar, ":.g", g)
      if (!trm %in% nm) stop("missing interaction term: ", trm)
      v[trm] <- 1
    }
    v
  }
  est_se <- function(v) {
    e  <- sum(v * cf)
    se <- sqrt(as.numeric(t(v) %*% V %*% v))
    c(est = e, se = se, t = e / se, p = 2 * pt(-abs(e / se), dfres))
  }
  
  per_group <- map_dfr(lev, function(g) {
    r <- est_se(svec(g))
    n <- sum(d$.g == g)
    tibble(outcome = yvar, group = g, n = n, slope = r["est"], se = r["se"],
           ci_lo = r["est"] - qt(.975, dfres) * r["se"],
           ci_hi = r["est"] + qt(.975, dfres) * r["se"],
           t = r["t"], p_value = r["p"],
           mean_X = mean(d[[xvar]][d$.g == g]), mean_Y = mean(d[[yvar]][d$.g == g]),
           pearson_r = cor(d[[xvar]][d$.g == g], d[[yvar]][d$.g == g]))
  })
  
  prs <- combn(lev, 2, simplify = FALSE)
  pairwise <- map_dfr(prs, function(pr) {
    r <- est_se(svec(pr[1]) - svec(pr[2]))
    tibble(outcome = yvar, group_1 = pr[1], group_2 = pr[2],
           slope_diff = r["est"], se = r["se"], t = r["t"], p_raw = r["p"])
  }) %>% mutate(p_FDR = p.adjust(p_raw, "fdr"),
                p_bonferroni = p.adjust(p_raw, "bonferroni"),
                sig_FDR = p_FDR < 0.05)
  
  means_test <- tibble(
    outcome = yvar, test = c("Group difference in X (ANOVA)", "Group difference in Y (ANOVA)"),
    F_value = c(summary(aov(d[[xvar]] ~ d$.g))[[1]][1, "F value"],
                summary(aov(d[[yvar]] ~ d$.g))[[1]][1, "F value"]),
    p_value = c(summary(aov(d[[xvar]] ~ d$.g))[[1]][1, "Pr(>F)"],
                summary(aov(d[[yvar]] ~ d$.g))[[1]][1, "Pr(>F)"]))
  
  list(global = global, per_group = per_group, pairwise = pairwise,
       means_test = means_test, model = m_int)
}

################################################################################
## 9. RUN EVERYTHING FOR EACH RECEPTOR VARIANT
################################################################################

all_xttests <- list()   # collects the X-axis t-tests across every receptor set

for (rs_name in names(receptor_sets)) {
  rs_cols <- receptor_sets[[rs_name]]
  
  message("\n=======================================================")
  message(sprintf(" PROCESSING RECEPTOR SET: %s", rs_name))
  message("=======================================================")
  
  # Recalculate X specifically for this receptor subset
  master_current <- master
  master_current$X_disease <- compute_x_value(master_current, rs_cols)
  
  # Healthy group mean on this x-axis -> where the green line goes
  ref_x <- mean(master_current$X_disease[master_current$Group == "H"], na.rm = TRUE)
  message(sprintf("  healthy reference line at x = %.5f", ref_x))
  
  figs      <- list()
  res_slope <- list()
  xtt       <- list()
  
  for (i in seq_len(nrow(OUTCOMES))) {
    yv <- OUTCOMES$var[i]
    yl <- OUTCOMES$label[i]
    tg <- OUTCOMES$tag[i]
    
    file_stub <- paste0("Fig1_", tg, "_", rs_name)

    ## Same figure drawn twice, once per stats placement, into the two folders
    ## STATS_MODES names. Filenames are identical in both, so anything that
    ## consumes them (Compose_Figure_Panels.R) needs only the folder swapped.
    for (sm in names(STATS_MODES)) {
      figs[[tg]] <- make_figure1(master_current, yvar = yv, ylab = yl, tag = tg,
                                 receptor_subset = rs_cols,
                                 file_stub = file_stub, rs_name = rs_name,
                                 ref_x = ref_x, stats_mode = sm)

      make_overlay(master_current, yvar = yv, ylab = yl, tag = tg,
                   file_stub = paste0(file_stub, "_overlay"), rs_name = rs_name,
                   ref_x = ref_x, stats_mode = sm)
    }

    ## Placement-independent: the t-test, slope models and centroids below are
    ## identical for both modes, so they are collected once.
    xtt[[tg]]       <- figs[[tg]]$ttest
    res_slope[[tg]] <- slope_stats(master_current, yvar = yv)
  }
  
  xttest_tbl   <- bind_rows(xtt)
  global_tbl   <- bind_rows(map(res_slope, "global"))
  pergroup_tbl <- bind_rows(map(res_slope, "per_group"))
  pairwise_tbl <- bind_rows(map(res_slope, "pairwise"))
  means_tbl    <- bind_rows(map(res_slope, "means_test"))
  
  all_xttests[[rs_name]] <- xttest_tbl
  
  cat(sprintf("\n----- X-AXIS T-TEST (%s vs %s) | receptor set: %s -----\n",
              FOCUS_GROUPS[1], FOCUS_GROUPS[2], rs_name))
  print(as.data.frame(xttest_tbl %>%
                        select(outcome, n_1, n_2, mean_X_1, mean_X_2,
                               mean_diff, t, df, p_value, cohens_d, significant)),
        row.names = FALSE, digits = 4)
  
  write_xlsx(list(Xaxis_ttest        = xttest_tbl,
                  Global_interaction = global_tbl,
                  Per_group_slopes   = pergroup_tbl,
                  Pairwise_slopes    = pairwise_tbl,
                  Group_mean_tests   = means_tbl),
             file.path(OUT_DIR, paste0("Stats_", rs_name, ".xlsx")))
  
  centroids <- map_dfr(names(figs), ~ figs[[.x]]$summary %>%
                         mutate(outcome = .x) %>%
                         select(outcome, Plot_Group, n, mean_x, mean_y, slope, p, r))
  write_xlsx(centroids, file.path(OUT_DIR, paste0("Centroids_", rs_name, ".xlsx")))
}

## One combined workbook with every X-axis t-test, all receptor sets in one sheet
xttest_all <- bind_rows(all_xttests)
write_xlsx(list(Xaxis_ttests_all_sets = xttest_all),
           file.path(OUT_DIR, "Xaxis_Ttests_ALL_receptor_sets.xlsx"))
write.csv(xttest_all, file.path(OUT_DIR, "Xaxis_Ttests_ALL_receptor_sets.csv"),
          row.names = FALSE)

message("\nDone. Everything written to: ", normalizePath(OUT_DIR))

################################################################################
## 10. FORMAT PUBLICATION-READY DEMOGRAPHICS TABLE  ->  outputs/
################################################################################

table_data <- read_excel(file.path(OUT_DIR, "Table1_Hydrocodone_vs_Tramadol.xlsx"))

pub_table <- table_data %>%
  gt() %>%
  cols_align(align = "left", columns = Characteristic) %>%
  cols_align(align = "center", columns = -Characteristic) %>%
  ## Column order is Characteristic | Hydrocodone | Tramadol | CBP+O | p-value,
  ## so these indices track the table built in section 6. Adding or removing a
  ## reference column means updating them.
  tab_spanner(label = "Opioid subgroup comparison", columns = 2:3) %>%
  tab_spanner(label = "Reference group",            columns = 4) %>%
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
                      " only. The CBP+O column is a descriptive reference group ",
                      "and is not included in the test."),
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
## 11. STANDALONE T-TEST: HT6 X VALUES (Hydrocodone vs Tramadol)
################################################################################

master_ht6 <- master
master_ht6$X_HT6 <- compute_x_value(master_ht6, c("R_5HT6"))

ht6_subset <- master_ht6 %>%
  filter(Plot_Group %in% FOCUS_GROUPS) %>%
  mutate(Plot_Group = droplevels(factor(Plot_Group, levels = FOCUS_GROUPS)))

ht6_ttest <- t.test(X_HT6 ~ Plot_Group, data = ht6_subset)

cat("\n===== T-TEST: HT6 values - Hydrocodone vs Tramadol =====\n")
print(ht6_ttest)


################################################################################
## 12. RECEPTOR GROUP MEANS — Healthy vs CBP-O
##
## Descriptive reference table: one row per receptor, giving the healthy and
## CBP-O means, their difference, and where CBP+O sits.
################################################################################

receptor_summary <- map_dfr(RECEPTORS, function(rc) {
  
  h <- master[[rc]][master$Group == "H"];     h <- h[!is.na(h)]
  o <- master[[rc]][master$Group == "CBP-O"]; o <- o[!is.na(o)]
  p <- master[[rc]][master$Group == "CBP+O"]; p <- p[!is.na(p)]
  
  gap     <- mean(o) - mean(h)
  sd_pool <- sd(c(h, o))
  tt      <- if (length(h) > 1 && length(o) > 1) t.test(o, h) else NULL
  
  tibble(
    Receptor      = sub("^R_", "", rc),
    n_H           = length(h),
    mean_H        = mean(h),
    sd_H          = if (length(h) > 1) sd(h) else NA_real_,
    n_CBPO        = length(o),
    mean_CBPO     = mean(o),
    sd_CBPO       = if (length(o) > 1) sd(o) else NA_real_,
    difference    = gap,
    t_value       = if (is.null(tt)) NA_real_ else unname(tt$statistic),
    p_value       = if (is.null(tt)) NA_real_ else tt$p.value,
    subject_SD    = sd_pool,
    gap_over_SD   = gap / sd_pool,
    n_CBPplusO    = length(p),
    mean_CBPplusO = if (length(p)) mean(p) else NA_real_,
    sd_CBPplusO   = if (length(p) > 1) sd(p) else NA_real_
  )
}) %>%
  mutate(p_FDR = p.adjust(p_value, "fdr"))

cat("\n===== RECEPTOR MEANS: Healthy vs CBP-O =====\n")
print(as.data.frame(receptor_summary), row.names = FALSE, digits = 4)

write_xlsx(list(Receptor_group_means = receptor_summary),
           file.path(OUT_DIR, "Receptor_Group_Means.xlsx"))
write.csv(receptor_summary, file.path(OUT_DIR, "Receptor_Group_Means.csv"),
          row.names = FALSE)

receptor_gt <- receptor_summary %>%
  transmute(
    Receptor,
    `Healthy (mean ± SD)`    = sprintf("%.4f \u00b1 %.4f", mean_H, sd_H),
    `CBP-O (mean ± SD)`      = sprintf("%.4f \u00b1 %.4f", mean_CBPO, sd_CBPO),
    `Difference (CBP-O − H)` = sprintf("%.4f", difference),
    `p (FDR)`                = map_chr(p_FDR, fmt_p),
    `CBP+O (mean ± SD)`      = sprintf("%.4f \u00b1 %.4f", mean_CBPplusO, sd_CBPplusO)
  ) %>%
  gt() %>%
  cols_align(align = "left",   columns = Receptor) %>%
  cols_align(align = "center", columns = -Receptor) %>%
  tab_style(style = cell_text(weight = "bold"),
            locations = cells_column_labels(everything())) %>%
  tab_style(style = cell_text(weight = "bold"),
            locations = cells_body(columns = Receptor)) %>%
  tab_footnote(
    footnote = paste0("Welch two-sample t-test, healthy vs CBP-O, ",
                      "Benjamini-Hochberg corrected across the 19 receptors."),
    locations = cells_column_labels(columns = `p (FDR)`)
  ) %>%
  tab_options(
    table.border.top.color = "black",    table.border.top.width = px(2),
    table.border.bottom.color = "black", table.border.bottom.width = px(2),
    column_labels.border.top.color = "black",    column_labels.border.top.width = px(2),
    column_labels.border.bottom.color = "black", column_labels.border.bottom.width = px(2),
    table_body.border.bottom.color = "black",    table_body.border.bottom.width = px(2),
    table_body.hlines.color = "transparent",
    table_body.vlines.color = "transparent"
  )

gtsave(receptor_gt, file.path(OUT_DIR, "Receptor_Group_Means.html"))

tryCatch(
  gtsave(receptor_gt, file.path(OUT_DIR, "Receptor_Group_Means.docx")),
  error = function(e) {
    message("gtsave(.docx) failed (", conditionMessage(e), ") — writing .rtf instead.")
    gtsave(receptor_gt, file.path(OUT_DIR, "Receptor_Group_Means.rtf"))
  }
)

message("\nReceptor group-means table written to: ", normalizePath(OUT_DIR))


################################################################################
## 14. OPIOID EXPOSURE (x) vs CLINICAL OUTCOMES (y)
##     — Hydrocodone group vs Tramadol group
##
## Predictors : log MME, log ROE
## Outcomes   : NRS, PC1, PC2, PC3
## => 8 combinations, each producing a faceted figure and an overlay figure,
##    plus one combined stats workbook.
##
## Self-contained: needs `master`, `FOCUS_GROUPS`, `OUT_DIR`, `custom_colors`,
## `custom_shapes`, `fmt_p`, and the packages already loaded.
################################################################################

## Predictors. Swap `var` to "MME" / "ROE" (and edit xlab) for untransformed axes.
PREDICTORS <- list(
  list(var = "log_MME", tag = "logMME",
       xlab = expression(log[10]~"MME (mg/day)")),
  list(var = "log_ROE", tag = "logROE",
       xlab = expression(log[10]~"ROE (mg/L)"))
)

XY_OUTCOMES <- tribble(
  ~var,  ~label,                                ~tag,
  "nrs", "Clinical Pain Score (NRS 0-10)",      "NRS",
  "PC1", "PC1 score — Functional disability",   "PC1",
  "PC2", "PC2 score — Pain quality/severity",   "PC2",
  "PC3", "PC3 score — Negative affect",         "PC3"
)

xy_slopes_all <- list()
xy_tests_all  <- list()

for (pr in PREDICTORS) {
  
  xv    <- pr$var
  xtag  <- pr$tag
  xlab  <- pr$xlab
  
  if (!xv %in% names(master)) {
    message(sprintf("Skipping %s — column not found in master", xv)); next
  }
  
  for (i in seq_len(nrow(XY_OUTCOMES))) {
    
    yv <- XY_OUTCOMES$var[i]
    yl <- XY_OUTCOMES$label[i]
    ytag <- XY_OUTCOMES$tag[i]
    combo <- paste0(xtag, "_vs_", ytag)
    
    d <- master %>%
      filter(Plot_Group %in% FOCUS_GROUPS,
             !is.na(.data[[xv]]), is.finite(.data[[xv]]),
             !is.na(.data[[yv]])) %>%
      mutate(Plot_Group = droplevels(factor(Plot_Group, levels = FOCUS_GROUPS)))
    
    if (nlevels(d$Plot_Group) < 2 || any(table(d$Plot_Group) < 3)) {
      message(sprintf("Skipping %s — fewer than 3 subjects in a group", combo)); next
    }
    
    ## ---- per-group regressions ------------------------------------------------
    slopes <- d %>%
      group_by(Plot_Group) %>%
      group_modify(~ {
        fit <- lm(.x[[yv]] ~ .x[[xv]])
        tibble(n         = nrow(.x),
               slope     = coef(fit)[2],
               intercept = coef(fit)[1],
               r         = cor(.x[[xv]], .x[[yv]]),
               p         = summary(fit)$coefficients[2, 4])
      }) %>% ungroup() %>%
      mutate(predictor = xtag, outcome = ytag,
             mean_x = map_dbl(Plot_Group, ~ mean(d[[xv]][d$Plot_Group == .x])),
             mean_y = map_dbl(Plot_Group, ~ mean(d[[yv]][d$Plot_Group == .x])),
             label_text = stats_block(n, r, slope, p),
             stat_line  = stats_line(n, r, slope, p),
             stat_strip = stats_strip(n, r, slope, p))

    ## Strip labels settled before the plot is built (see make_figure1).
    d <- strip_facet_data(
      d, "Plot_Group", as.character(slopes$Plot_Group),
      if (stats_fit_one_line(slopes$stat_line, XY_W, 2))
        slopes$stat_line else slopes$stat_strip)

    ## ---- group comparisons ----------------------------------------------------
    tt_x  <- t.test(d[[xv]] ~ d$Plot_Group)      # exposure difference
    tt_y  <- t.test(d[[yv]] ~ d$Plot_Group)      # outcome difference
    m_add <- lm(d[[yv]] ~ d[[xv]] + d$Plot_Group)
    m_int <- lm(d[[yv]] ~ d[[xv]] * d$Plot_Group)
    av    <- anova(m_add, m_int)                 # slope difference
    
    tests <- tibble(
      predictor = xtag, outcome = ytag,
      test      = c(paste0(xtag, " difference (Welch t)"),
                    paste0(ytag, " difference (Welch t)"),
                    "Slope difference (x by Group interaction)"),
      statistic = c(unname(tt_x$statistic), unname(tt_y$statistic), av$F[2]),
      df        = c(unname(tt_x$parameter), unname(tt_y$parameter), av$Df[2]),
      p_value   = c(tt_x$p.value, tt_y$p.value, av$`Pr(>F)`[2])
    )
    
    xy_slopes_all[[combo]] <- slopes
    xy_tests_all[[combo]]  <- tests
    
    cat(sprintf("\n===== %s vs %s: Hydrocodone vs Tramadol =====\n", xtag, ytag))
    print(as.data.frame(slopes %>% select(Plot_Group, n, mean_x, mean_y, slope, r, p)),
          row.names = FALSE, digits = 4)
    cat("\n")
    print(as.data.frame(tests %>% select(test, statistic, df, p_value)),
          row.names = FALSE, digits = 4)
    
    ## Wrapped by hand, one group per line. ggplot does not wrap captions, and at
    ## TEXT_SCALE = 2 the old single-line version ran off the right edge of the
    ## device and was silently clipped.
    cap <- sprintf(
      "%s — %s: mean = %.3f (n=%d)\n%s — %s: mean = %.3f (n=%d)\nWelch t(%.1f) = %.2f, p = %s\nSlope difference (interaction): F(%d, %d) = %.2f, p = %s",
      xtag, levels(d$Plot_Group)[1], slopes$mean_x[1], slopes$n[1],
      xtag, levels(d$Plot_Group)[2], slopes$mean_x[2], slopes$n[2],
      unname(tt_x$parameter), unname(tt_x$statistic), fmt_p(tt_x$p.value),
      av$Df[2], av$Res.Df[2], av$F[2], fmt_p(av$`Pr(>F)`[2]))
    
    ## ---- figures, one pair per stats placement --------------------------------
    ## No title and no subtitle: the axis labels and the strip labels already say
    ## which predictor, which outcome and which group each panel shows.
    for (sm in names(STATS_MODES)) {

      ## ---- faceted figure ----------------------------------------------------
      p_facet <- ggplot(d, aes(.data[[xv]], .data[[yv]],
                               color = Plot_Group, shape = Plot_Group)) +
        geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "grey30",
                    linetype = "solid", linewidth = FIG_LINE_SIZE,
                    fill = "grey80", alpha = 0.4) +
        geom_point(size = FIG_POINT_SIZE, alpha = 0.85) +
        scale_color_manual(values = custom_colors) +
        scale_shape_manual(values = custom_shapes) +
        labs(x = xlab, y = wrap_lab(yl, 22), caption = cap) +
        theme_fig(faceted = TRUE)

      if (sm == "header") {
        ## Bold group name over a smaller, plain stats line; the weights come
        ## from the plotmath label, so the strip element itself is plain.
        p_facet <- p_facet + facet_wrap(~ .strip, labeller = label_parsed) +
          theme(strip.text = element_text(face = "plain", size = FIG_STRIP_SIZE,
                                          lineheight = 1.05,
                                          margin = margin(b = 6, t = 2)))
      } else {
        bf     <- box_fracs(p_facet + facet_wrap(~ Plot_Group), XY_W, XY_H,
                            slopes$label_text)
        lp     <- place_corner_labels(d, "Plot_Group", xv, yv,
                                      w_frac = bf[["w"]], h_frac = bf[["h"]])
        audit_corner(paste0("Fig_", combo, "_faceted"), lp)
        lab_df <- left_join(lp$coords, slopes, by = "Plot_Group")
        p_facet <- p_facet +
          facet_wrap(~ Plot_Group) +
          geom_label(data = lab_df, aes(x = lab_x, y = lab_y, label = label_text),
                     inherit.aes = FALSE,
                     hjust = lab_df$hjust, vjust = lab_df$vjust,
                     size = FIG_ANNOT_SIZE, lineheight = 1.05,
                     fontface = "bold", color = "black", fill = "white", alpha = 1,
                     linewidth = 0, label.padding = unit(0.35, "lines")) +
          coord_cartesian(ylim = lp$ylim)
      }

      ggsave(fig_path(sm, paste0("Fig_", combo, "_faceted")),
             p_facet, width = XY_W, height = XY_H, dpi = 300)

      ## ---- overlay figure ----------------------------------------------------
      p_overlay <- ggplot(d, aes(.data[[xv]], .data[[yv]],
                                 color = Plot_Group, shape = Plot_Group)) +
        geom_point(alpha = 0.6, size = FIG_POINT_SIZE) +
        geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linetype = "solid",
                    alpha = 0.15, linewidth = FIG_LINE_SIZE) +
        scale_shape_manual(values = custom_shapes, guide = "none") +
        labs(x = xlab, y = wrap_lab(yl, 22), color = NULL, shape = NULL,
             caption = cap) +
        theme_fig(faceted = FALSE, legend = "bottom")

      if (sm == "header") {
        key_labs <- setNames(paste0(as.character(slopes$Plot_Group), ": ",
                                    slopes$stat_line),
                             as.character(slopes$Plot_Group))
        p_overlay <- p_overlay +
          scale_color_manual(values = custom_colors, labels = key_labs) +
          guides(color = guide_legend(ncol = 1)) +
          ## See make_overlay: a centred key this long gets clipped, not wrapped.
          theme(legend.justification = "left",
                legend.text          = element_text(size = FIG_STRIP_STATS,
                                                   face = "plain"),
                legend.margin        = margin(0, 0, 0, 0),
                legend.box.margin    = margin(0, 0, 0, 0),
                legend.key.spacing.y = unit(2, "pt"))
      } else {
        box <- paste(sprintf("%s: %s", as.character(slopes$Plot_Group),
                             slopes$stat_line), collapse = "\n")
        bf  <- box_fracs(p_overlay, XYO_W, XYO_H, box)
        lp  <- place_corner_labels(d %>% mutate(.all = factor("all")), ".all",
                                   xv, yv, extent_var = "Plot_Group",
                                   w_frac = bf[["w"]], h_frac = bf[["h"]])
        audit_corner(paste0("Fig_", combo, "_overlay"), lp)
        p_overlay <- p_overlay +
          scale_color_manual(values = custom_colors) +
          annotate("label", x = lp$coords$lab_x[1], y = lp$coords$lab_y[1],
                   hjust = lp$coords$hjust[1], vjust = lp$coords$vjust[1],
                   label = box, size = FIG_ANNOT_SIZE, lineheight = 1.05,
                   fontface = "bold", color = "black", fill = "white",
                   linewidth = 0, label.padding = unit(0.35, "lines")) +
          coord_cartesian(ylim = lp$ylim)
      }

      ggsave(fig_path(sm, paste0("Fig_", combo, "_overlay")),
             p_overlay, width = XYO_W, height = XYO_H, dpi = 300)
    }
  }
}

## ---- combined exports --------------------------------------------------------
xy_slopes_tbl <- bind_rows(xy_slopes_all) %>%
  select(predictor, outcome, Plot_Group, n, mean_x, mean_y, slope, intercept, r, p)

## FDR applied within each predictor x test-type family, across the 4 outcomes
xy_tests_tbl <- bind_rows(xy_tests_all) %>%
  group_by(predictor, test) %>%
  mutate(p_FDR = p.adjust(p_value, "fdr")) %>%
  ungroup() %>%
  select(predictor, outcome, test, statistic, df, p_value, p_FDR)

cat("\n===== SLOPE-DIFFERENCE TESTS (FDR across the 4 outcomes) =====\n")
print(as.data.frame(xy_tests_tbl %>%
                      filter(grepl("Slope difference", test)) %>%
                      select(predictor, outcome, statistic, df, p_value, p_FDR)),
      row.names = FALSE, digits = 4)

write_xlsx(list(Per_group_regression = xy_slopes_tbl,
                Group_comparisons    = xy_tests_tbl,
                Data_used            = master %>%
                  filter(Plot_Group %in% FOCUS_GROUPS) %>%
                  select(PIN, Plot_Group, MME, log_MME, ROE, log_ROE,
                         nrs, PC1, PC2, PC3, DOU, MQS_total)),
           file.path(OUT_DIR, "Stats_Exposure_vs_Outcomes.xlsx"))

message("\nMME/ROE vs NRS/PC1/PC2/PC3 figures and stats written to: ",
        normalizePath(OUT_DIR))

## ---- 14a. CORNER-PLACEMENT AUDIT ---------------------------------------------
## Every figure in figs_stats_corner/ put its stats box somewhere; this reports
## whether any of them ended up over a data point or a confidence ribbon. The
## whole point of the corner placement is that this comes out zero.
corner_audit_tbl <- if (length(CORNER_AUDIT))
  do.call(rbind, CORNER_AUDIT) else NULL

if (!is.null(corner_audit_tbl)) {
  covered <- corner_audit_tbl$pts_covered > 0 | corner_audit_tbl$rib_covered > 0

  cat("\n===== STATS-BOX PLACEMENT AUDIT (figs_stats_corner) =====\n")
  cat(sprintf("panels placed: %d\n", nrow(corner_audit_tbl)))
  print(table(corner = corner_audit_tbl$corner))
  cat(sprintf("panels needing an added blank band: %d\n",
              sum(corner_audit_tbl$pad_top > 0 | corner_audit_tbl$pad_bot > 0)))
  cat(sprintf("panels where the box covers data: %d\n", sum(covered)))

  if (any(covered)) {
    warning("Stats box overlaps data in ", sum(covered),
            " panel(s) -- see the Corner_placement_audit sheet.", call. = FALSE)
    print(corner_audit_tbl[covered, c("figure", "panel", "corner",
                                      "pts_covered", "rib_covered")],
          row.names = FALSE)
  } else {
    cat("No stats box covers a data point or a confidence ribbon.\n")
  }

  write_xlsx(list(Corner_placement_audit = corner_audit_tbl),
             file.path(OUT_DIR, "Figure_StatsBox_Placement_Audit.xlsx"))
}

################################################################################
## 15. MEDIATION ANALYSIS
##     X = log ROE   ->   M = 5-HT4 / 5-HT6 related activity   ->   Y = NRS
##
## MODEL
##   a  : X -> M            (exposure changes receptor-related activity)
##   b  : M -> Y | X        (receptor-related activity changes pain,
##                           holding exposure constant)
##   c' : X -> Y | M        (direct effect of exposure not routed through M)
##   c  : X -> Y            (total effect = c' + ab)
##   ab : indirect / mediated effect  <- the quantity of interest
##
## WHY EVERYTHING IS ESTIMATED SEPARATELY BY DRUG AND BY RECEPTOR
##   The c path reverses sign between subgroups (hydrocodone +, tramadol -,
##   interaction p = 0.008). A product ab that is positive in one group and
##   negative in the other averages toward zero when the groups are pooled, so
##   a pooled model would report "no mediation" even if both groups mediated
##   strongly in opposite directions. Formally: drug moderates the mediation.
##   Fitting one model per group is the stratified estimate of that moderated
##   mediation; section 15d tests the moderation itself.
##
## SCALING
##   ROE differs ~40x in scale between the two drugs, so raw-unit a, b and ab
##   are not comparable across groups. All primary estimates are standardized
##   WITHIN GROUP (z-scored on the estimation sample), which puts every path in
##   SD units and makes the hydrocodone/tramadol contrast meaningful.
##   Unstandardized point estimates are reported alongside for reference.
##
## INFERENCE
##   ab is a product of two normal-ish estimates and is not itself normal, so
##   CIs come from a nonparametric bootstrap (percentile method), not from a
##   Sobel test. Each replicate resamples subjects with replacement and refits
##   all three regressions. With n = 17-19 per cell, BCa is unreliable, so
##   percentile intervals are used deliberately.
##
## CAUSAL CAVEAT
##   These are cross-sectional data. Mediation notation implies a causal
##   ordering that the design cannot establish; treat ab as "variance in pain
##   statistically routed through receptor activity", not as a demonstrated
##   mechanism. Sequential-ignorability (no unmeasured M-Y confounding) is
##   untestable here.
##
## Self-contained: needs `master`, `FOCUS_GROUPS`, `OUT_DIR` and the packages
## already loaded at the top of the master script (dplyr, purrr, tibble,
## ggplot2, writexl). No new packages required.
################################################################################

## ---- 15a. CONFIG -------------------------------------------------------------

MED_X       <- "log_ROE"                      # exposure
MED_M       <- c("R_5HT4", "R_5HT6", "R_MOR")          # mediators, fitted one at a time
MED_Y       <- c("nrs")                       # outcome(s); add "PC2" if wanted
MED_GROUPS  <- FOCUS_GROUPS                   # Hydrocodone group / Tramadol group

## Covariates entered into all three regressions. Leave empty for the primary
## model: with n = 17-19 per cell, each covariate costs a df you cannot spare.
## Sex is the obvious candidate (87% vs 52% female) but is nearly collinear
## with group -- run it as a sensitivity check, not as the headline model.
MED_COVARS  <- character(0)                   # e.g. c("age", "female")

MED_BOOT    <- 5000                           # bootstrap replicates
MED_SEED    <- 42
MED_MIN_N   <- 12                             # skip cells smaller than this

## Pretty labels for figures
MED_LABELS <- c(log_ROE = "log10 ROE",
                log_MME = "log10 MME",
                R_MOR   = "MOR activity",
                R_5HT1A = "5-HT1A activity",
                R_5HT4  = "5-HT4 activity",
                R_5HT6  = "5-HT6 activity",
                nrs     = "NRS pain",
                PC2     = "PC2 pain quality")
med_lab <- function(v) if (!is.na(MED_LABELS[v])) unname(MED_LABELS[v]) else v

MED_DIR <- file.path(OUT_DIR, "mediation")
dir.create(MED_DIR, showWarnings = FALSE)

set.seed(MED_SEED)

## fmt_p may already exist from earlier in the script; define a fallback.
if (!exists("fmt_p")) {
  fmt_p <- function(p) ifelse(is.na(p), "NA",
                              ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}


## ---- 15b. CORE MACHINERY -----------------------------------------------------

## z-score the listed columns; returns NULL if any column is constant, which
## can happen in a degenerate bootstrap replicate.
med_scale <- function(df, vars) {
  for (v in vars) {
    s <- stats::sd(df[[v]])
    if (!is.finite(s) || s == 0) return(NULL)
    df[[v]] <- (df[[v]] - mean(df[[v]])) / s
  }
  df
}

## Look up a coefficient by variable name. R strips backticks from coefficient
## names when the variable name is syntactically valid, but keeps them when it
## is not, so a formula written with backticks can yield either form. This
## tries both instead of assuming.
med_pick <- function(nms, v) {
  hit <- match(v, nms)
  if (is.na(hit)) hit <- match(paste0("`", v, "`"), nms)
  if (is.na(hit)) stop("No coefficient found for variable: ", v,
                       "  (available: ", paste(nms, collapse = ", "), ")")
  hit
}

## Fit the three regressions and return every path coefficient.
med_fit_paths <- function(df, xv, mv, yv, covs = character(0)) {
  cov_rhs <- if (length(covs)) paste("+", paste(covs, collapse = " + ")) else ""

  m_a <- lm(as.formula(sprintf("`%s` ~ `%s` %s", mv, xv, cov_rhs)), data = df)
  m_b <- lm(as.formula(sprintf("`%s` ~ `%s` + `%s` %s", yv, mv, xv, cov_rhs)), data = df)
  m_c <- lm(as.formula(sprintf("`%s` ~ `%s` %s", yv, xv, cov_rhs)), data = df)
  
  ca <- summary(m_a)$coefficients
  cb <- summary(m_b)$coefficients
  cc <- summary(m_c)$coefficients
  
  ia <- med_pick(rownames(ca), xv)
  ib <- med_pick(rownames(cb), mv)
  ic <- med_pick(rownames(cb), xv)
  id <- med_pick(rownames(cc), xv)
  
  a  <- ca[ia, 1];  a_se  <- ca[ia, 2];  a_p  <- ca[ia, 4]
  b  <- cb[ib, 1];  b_se  <- cb[ib, 2];  b_p  <- cb[ib, 4]
  cp <- cb[ic, 1];  cp_se <- cb[ic, 2];  cp_p <- cb[ic, 4]
  ct <- cc[id, 1];  c_se  <- cc[id, 2];  c_p  <- cc[id, 4]
  
  ## SE of ab by the multivariate delta method (Sobel). Reported for reference
  ## only -- it assumes ab is normal, which it is not; the bootstrap CI is the
  ## interval to trust and to report.
  ab_se <- sqrt(b^2 * a_se^2 + a^2 * b_se^2)
  
  list(a = a, b = b, cp = cp, c = ct, ab = a * b,
       a_se = a_se, b_se = b_se, cp_se = cp_se, c_se = c_se, ab_se = ab_se,
       a_p = a_p, b_p = b_p, cp_p = cp_p, c_p = c_p,
       r2_m = summary(m_a)$r.squared, r2_y = summary(m_b)$r.squared)
}

## Nonparametric bootstrap of all paths. Re-standardizes inside every replicate
## so the interval reflects the standardized quantity actually being reported.
med_boot <- function(df, xv, mv, yv, covs, R, standardize = TRUE) {
  n   <- nrow(df)
  out <- matrix(NA_real_, nrow = R, ncol = 5,
                dimnames = list(NULL, c("a", "b", "cp", "c", "ab")))
  for (i in seq_len(R)) {
    dd <- df[sample.int(n, n, replace = TRUE), , drop = FALSE]
    if (standardize) {
      dd <- med_scale(dd, c(xv, mv, yv))
      if (is.null(dd)) next
    }
    res <- tryCatch(med_fit_paths(dd, xv, mv, yv, covs), error = function(e) NULL)
    if (!is.null(res)) out[i, ] <- c(res$a, res$b, res$cp, res$c, res$ab)
  }
  out[stats::complete.cases(out), , drop = FALSE]
}

## Percentile CI plus a two-sided bootstrap p (proportion of replicates on the
## far side of zero, doubled).
med_ci <- function(v, conf = 0.95) {
  v <- v[is.finite(v)]
  if (!length(v)) return(c(lo = NA, hi = NA, p = NA))
  q <- stats::quantile(v, c((1 - conf) / 2, 1 - (1 - conf) / 2), names = FALSE)
  p <- 2 * min(mean(v <= 0), mean(v >= 0))
  c(lo = q[1], hi = q[2], p = min(p, 1))
}


## ---- 15c. SINGLE-MEDIATOR MODELS, ONE PER GROUP x RECEPTOR -------------------

med_rows   <- list()
med_bootd  <- list()
med_data   <- list()

for (yv in MED_Y) {
  for (g in MED_GROUPS) {
    for (mv in MED_M) {
      
      cell <- sprintf("%s | %s -> %s -> %s", g, MED_X, mv, yv)
      keep <- c(MED_X, mv, yv, MED_COVARS)
      
      d <- master[master$Plot_Group == g, , drop = FALSE]
      d <- d[stats::complete.cases(d[, keep, drop = FALSE]), keep, drop = FALSE]
      d <- d[apply(d, 1, function(r) all(is.finite(r))), , drop = FALSE]
      
      if (nrow(d) < MED_MIN_N) {
        message(sprintf("Skipping %s -- n = %d < MED_MIN_N", cell, nrow(d)))
        next
      }
      
      ## point estimates, unstandardized and standardized
      raw <- med_fit_paths(d, MED_X, mv, yv, MED_COVARS)
      ds  <- med_scale(d, c(MED_X, mv, yv))
      std <- med_fit_paths(ds, MED_X, mv, yv, MED_COVARS)
      
      ## bootstrap the standardized solution
      bd <- med_boot(d, MED_X, mv, yv, MED_COVARS, MED_BOOT, standardize = TRUE)
      ci_ab <- med_ci(bd[, "ab"]);  ci_a  <- med_ci(bd[, "a"])
      ci_b  <- med_ci(bd[, "b"]);   ci_cp <- med_ci(bd[, "cp"])
      ci_c  <- med_ci(bd[, "c"])
      
      med_rows[[cell]] <- tibble(
        outcome   = yv,
        group     = g,
        mediator  = mv,
        n         = nrow(d),
        ## standardized (primary). SE = model standard error;
        ## lo/hi = bootstrap percentile 95% CI.
        a_std     = std$a,  a_se  = std$a_se,  a_p   = std$a_p,
        a_lo      = ci_a[["lo"]],  a_hi  = ci_a[["hi"]],
        b_std     = std$b,  b_se  = std$b_se,  b_p   = std$b_p,
        b_lo      = ci_b[["lo"]],  b_hi  = ci_b[["hi"]],
        cprime_std= std$cp, cp_se = std$cp_se, cp_p  = std$cp_p,
        cp_lo     = ci_cp[["lo"]], cp_hi = ci_cp[["hi"]],
        c_total_std = std$c, c_se = std$c_se,  c_p   = std$c_p,
        c_lo      = ci_c[["lo"]],  c_hi  = ci_c[["hi"]],
        ab_std    = std$ab,
        ab_se_delta = std$ab_se,
        ab_lo     = ci_ab[["lo"]],
        ab_hi     = ci_ab[["hi"]],
        ab_boot_p = ci_ab[["p"]],
        ab_sig    = !is.na(ci_ab[["lo"]]) & (ci_ab[["lo"]] * ci_ab[["hi"]] > 0),
        prop_med  = ifelse(abs(std$c) > 0.05, std$ab / std$c, NA_real_),
        ## unstandardized (reference)
        a_raw = raw$a, b_raw = raw$b, cprime_raw = raw$cp,
        c_total_raw = raw$c, ab_raw = raw$ab,
        R2_mediator = std$r2_m, R2_outcome = std$r2_y,
        n_boot_ok = nrow(bd)
      )
      
      med_bootd[[cell]] <- bd
      med_data[[cell]]  <- cbind(d, group = g, mediator = mv, outcome = yv)
      
      cat(sprintf("\n===== MEDIATION: %s =====\n", cell))
      cat(sprintf("n = %d\n", nrow(d)))
      cat(sprintf("  a  (X->M)      = %7.3f  SE %.3f  p = %s\n",
                  std$a,  std$a_se,  fmt_p(std$a_p)))
      cat(sprintf("  b  (M->Y|X)    = %7.3f  SE %.3f  p = %s\n",
                  std$b,  std$b_se,  fmt_p(std$b_p)))
      cat(sprintf("  c' (X->Y|M)    = %7.3f  SE %.3f  p = %s   <- DIRECT effect\n",
                  std$cp, std$cp_se, fmt_p(std$cp_p)))
      cat(sprintf("  c  (X->Y)      = %7.3f  SE %.3f  p = %s   <- TOTAL effect\n",
                  std$c,  std$c_se,  fmt_p(std$c_p)))
      cat(sprintf("  ab (indirect)  = %7.3f  95%% CI [%.3f, %.3f]  %s\n",
                  std$ab, ci_ab[["lo"]], ci_ab[["hi"]],
                  ifelse(ci_ab[["lo"]] * ci_ab[["hi"]] > 0,
                         "*** CI excludes zero", "CI includes zero")))
      cat(sprintf("  check: c' + ab = %.3f  (should equal c = %.3f)\n",
                  std$cp + std$ab, std$c))
    }
  }
}

med_tbl <- bind_rows(med_rows)


## ---- 15d. MODERATED MEDIATION: does the mediation differ by drug? ------------
##
## The stratified models in 15c show two separate pictures but never test
## whether they differ. This does. The index of moderated mediation is
##      ab(hydrocodone) - ab(tramadol)
## bootstrapped with resampling STRATIFIED WITHIN GROUP so both cell sizes are
## preserved in every replicate. Path differences (a, b, c') are bootstrapped
## the same way. Standardization is done within group inside each replicate,
## for the scale reason given in the header.

modmed_rows <- list()

for (yv in MED_Y) {
  for (mv in MED_M) {
    
    keep <- c(MED_X, mv, yv, MED_COVARS)
    dg <- lapply(MED_GROUPS, function(g) {
      d <- master[master$Plot_Group == g, , drop = FALSE]
      d <- d[stats::complete.cases(d[, keep, drop = FALSE]), keep, drop = FALSE]
      d[apply(d, 1, function(r) all(is.finite(r))), , drop = FALSE]
    })
    names(dg) <- MED_GROUPS
    if (any(vapply(dg, nrow, 1L) < MED_MIN_N)) next
    
    one_rep <- function(resample) {
      est <- lapply(dg, function(d) {
        dd <- if (resample) d[sample.int(nrow(d), nrow(d), replace = TRUE), , drop = FALSE] else d
        dd <- med_scale(dd, c(MED_X, mv, yv))
        if (is.null(dd)) return(NULL)
        tryCatch(med_fit_paths(dd, MED_X, mv, yv, MED_COVARS), error = function(e) NULL)
      })
      if (any(vapply(est, is.null, TRUE))) return(rep(NA_real_, 4))
      c(ab = est[[1]]$ab - est[[2]]$ab,
        a  = est[[1]]$a  - est[[2]]$a,
        b  = est[[1]]$b  - est[[2]]$b,
        cp = est[[1]]$cp - est[[2]]$cp)
    }
    
    pt <- one_rep(FALSE)
    names(pt) <- c("ab", "a", "b", "cp")
    bs <- t(vapply(seq_len(MED_BOOT), function(i) one_rep(TRUE), numeric(4)))
    colnames(bs) <- c("ab", "a", "b", "cp")
    bs <- bs[stats::complete.cases(bs), , drop = FALSE]
    
    ci <- lapply(c("ab", "a", "b", "cp"), function(k) med_ci(bs[, k]))
    names(ci) <- c("ab", "a", "b", "cp")
    
    modmed_rows[[paste(yv, mv)]] <- tibble(
      outcome = yv, mediator = mv,
      contrast = paste(MED_GROUPS[1], "-", MED_GROUPS[2]),
      n_1 = nrow(dg[[1]]), n_2 = nrow(dg[[2]]),
      index_modmed = pt[["ab"]],
      index_lo = ci$ab[["lo"]], index_hi = ci$ab[["hi"]], index_p = ci$ab[["p"]],
      a_diff = pt[["a"]],  a_lo  = ci$a[["lo"]],  a_hi  = ci$a[["hi"]],  a_p  = ci$a[["p"]],
      b_diff = pt[["b"]],  b_lo  = ci$b[["lo"]],  b_hi  = ci$b[["hi"]],  b_p  = ci$b[["p"]],
      cp_diff = pt[["cp"]], cp_lo = ci$cp[["lo"]], cp_hi = ci$cp[["hi"]], cp_p = ci$cp[["p"]],
      n_boot_ok = nrow(bs)
    )
  }
}

modmed_tbl <- bind_rows(modmed_rows)

if (nrow(modmed_tbl)) {
  cat("\n===== INDEX OF MODERATED MEDIATION (bootstrap) =====\n")
  print(as.data.frame(modmed_tbl %>%
                        select(outcome, mediator, index_modmed, index_lo, index_hi, index_p,
                               b_diff, b_p, cp_diff, cp_p)),
        row.names = FALSE, digits = 3)
}


## ---- 15e. PARALLEL MEDIATORS: 5-HT4 AND 5-HT6 IN ONE MODEL -------------------
##
## The PI's instruction is to keep the receptors separate, and 15c does that.
## This block is the supplement that answers the obvious follow-up question:
## the two receptors are correlated (r ~ 0.30 hydrocodone, ~0.44 tramadol), so
## do they carry the SAME variance in pain or different variance? Both mediators
## enter the b-model together; each b is then adjusted for the other receptor.
## Report this as a robustness check, not as the primary analysis.

par_rows <- list()

for (yv in MED_Y) {
  for (g in MED_GROUPS) {
    
    keep <- c(MED_X, MED_M, yv, MED_COVARS)
    d <- master[master$Plot_Group == g, , drop = FALSE]
    d <- d[stats::complete.cases(d[, keep, drop = FALSE]), keep, drop = FALSE]
    d <- d[apply(d, 1, function(r) all(is.finite(r))), , drop = FALSE]
    if (nrow(d) < MED_MIN_N) next
    
    par_paths <- function(dd) {
      dd <- med_scale(dd, c(MED_X, MED_M, yv))
      if (is.null(dd)) return(NULL)
      cov_rhs <- if (length(MED_COVARS)) paste("+", paste(MED_COVARS, collapse = " + ")) else ""
      as_ <- vapply(MED_M, function(mv) {
        ca <- coef(lm(as.formula(sprintf("`%s` ~ `%s` %s", mv, MED_X, cov_rhs)),
                      data = dd))
        unname(ca[med_pick(names(ca), MED_X)])
      }, numeric(1))
      
      mb <- lm(as.formula(sprintf("`%s` ~ %s + `%s` %s", yv,
                                  paste0("`", MED_M, "`", collapse = " + "),
                                  MED_X, cov_rhs)), data = dd)
      cb  <- coef(mb)
      bs_ <- vapply(MED_M, function(mv) unname(cb[med_pick(names(cb), mv)]),
                    numeric(1))
      
      c(setNames(as_ * bs_, paste0("ab_", MED_M)),
        setNames(as_,       paste0("a_",  MED_M)),
        setNames(bs_,       paste0("b_",  MED_M)),
        cp       = unname(cb[med_pick(names(cb), MED_X)]),
        ab_total = sum(as_ * bs_))
    }
    
    pt <- par_paths(d)
    if (is.null(pt)) next
    bs <- t(vapply(seq_len(MED_BOOT), function(i) {
      r <- tryCatch(par_paths(d[sample.int(nrow(d), nrow(d), replace = TRUE), , drop = FALSE]),
                    error = function(e) NULL)
      if (is.null(r)) rep(NA_real_, length(pt)) else r
    }, numeric(length(pt))))
    colnames(bs) <- names(pt)
    bs <- bs[stats::complete.cases(bs), , drop = FALSE]
    
    par_rows[[paste(yv, g)]] <- bind_rows(lapply(names(pt), function(k) {
      ci <- med_ci(bs[, k])
      tibble(outcome = yv, group = g, n = nrow(d), quantity = k,
             estimate = unname(pt[[k]]),
             lo = ci[["lo"]], hi = ci[["hi"]], boot_p = ci[["p"]])
    }))
  }
}

parallel_tbl <- bind_rows(par_rows)


## ---- 15f. DIAGNOSTICS --------------------------------------------------------

diag_tbl <- bind_rows(lapply(MED_GROUPS, function(g) {
  keep <- c(MED_X, MED_M, MED_Y)
  d <- master[master$Plot_Group == g, , drop = FALSE]
  d <- d[stats::complete.cases(d[, keep, drop = FALSE]), keep, drop = FALSE]
  if (!nrow(d)) return(NULL)
  cmb <- utils::combn(keep, 2, simplify = FALSE)
  bind_rows(lapply(cmb, function(p) tibble(
    group = g, n = nrow(d), var_1 = p[1], var_2 = p[2],
    pearson_r = stats::cor(d[[p[1]]], d[[p[2]]]),
    p_value   = stats::cor.test(d[[p[1]]], d[[p[2]]])$p.value)))
}))


## ---- 15g. PATH DIAGRAMS ------------------------------------------------------
## One diagram per cell. Solid arrow = p < .05, dashed = not significant.
## The indirect effect and its bootstrap CI sit under the title.

## Every path is annotated as  coefficient (SE), p.
## c' is the DIRECT effect; the TOTAL effect c is printed under the X->Y arrow
## so the decomposition c = c' + ab is readable off the figure.
med_path_diagram <- function(row) {
  sig <- function(p) if (!is.na(p) && p < 0.05) "solid" else "22"
  
  nodes <- tibble(
    x = c(0.10, 0.50, 0.90), y = c(0.28, 0.82, 0.28),
    lab = c(med_lab(MED_X), med_lab(row$mediator), med_lab(row$outcome)))
  
  edges <- tibble(
    x    = c(0.16, 0.56, 0.20), y    = c(0.35, 0.76, 0.28),
    xend = c(0.44, 0.84, 0.80), yend = c(0.76, 0.35, 0.28),
    lty  = c(sig(row$a_p), sig(row$b_p), sig(row$cp_p)),
    lab  = c(sprintf("a = %.2f (SE %.2f)\np = %s",
                     row$a_std,  row$a_se,  fmt_p(row$a_p)),
             sprintf("b = %.2f (SE %.2f)\np = %s",
                     row$b_std,  row$b_se,  fmt_p(row$b_p)),
             sprintf("c' = %.2f (SE %.2f), p = %s   [DIRECT]",
                     row$cprime_std, row$cp_se, fmt_p(row$cp_p))),
    lx = c(0.25, 0.75, 0.50), ly = c(0.60, 0.60, 0.205))
  
  ## total effect, printed below the direct effect
  tot <- sprintf("c  = %.2f (SE %.2f), p = %s   [TOTAL = c' + ab]",
                 row$c_total_std, row$c_se, fmt_p(row$c_p))
  
  ggplot() +
    geom_segment(data = edges,
                 aes(x = x, y = y, xend = xend, yend = yend, linetype = I(lty)),
                 linewidth = 0.8, colour = "grey25",
                 arrow = arrow(length = unit(3.2, "mm"), type = "closed")) +
    geom_label(data = nodes, aes(x, y, label = lab),
               size = 4.6, fontface = "bold", label.padding = unit(4, "mm"),
               fill = "white", colour = "black") +
    geom_text(data = edges, aes(lx, ly, label = lab),
              size = 3.8, lineheight = 0.95, colour = "grey15") +
    annotate("text", x = 0.50, y = 0.135, label = tot,
             size = 3.8, colour = "grey15") +
    labs(title = sprintf("%s  --  %s", row$group, med_lab(row$mediator)),
         subtitle = sprintf(
           "n = %d   |   INDIRECT ab = %.3f, bootstrap 95%% CI [%.3f, %.3f]%s",
           row$n, row$ab_std, row$ab_lo, row$ab_hi,
           ifelse(row$ab_sig, "  *", "  (contains zero)")),
         caption = paste0(
           "Standardized paths (within-group z-scores). Solid arrow = p < .05.\n",
           "SE = model standard error. The ab interval is a bootstrap percentile CI, ",
           "not SE-based:\nab is a product of two estimates and is not normally distributed.")) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0.08, 1), expand = FALSE) +
    theme_void(base_size = 13) +
    theme(plot.title    = element_text(face = "bold", hjust = 0, margin = margin(b = 2)),
          plot.subtitle = element_text(hjust = 0, margin = margin(b = 10)),
          plot.caption  = element_text(hjust = 0, colour = "grey40", size = 8),
          plot.margin   = margin(12, 12, 10, 12))
}

for (i in seq_len(nrow(med_tbl))) {
  row <- med_tbl[i, ]
  p <- med_path_diagram(row)
  ggsave(file.path(MED_DIR, sprintf("Mediation_path_%s_%s_%s.png",
                                    row$outcome, gsub(" ", "", row$group), row$mediator)),
         p, width = 7.5, height = 5.2, dpi = 300)
}

## Effect decomposition: total, direct and indirect side by side, each with its
## bootstrap CI. This is the figure that answers "where did the effect go" --
## when c' sits on top of c and ab sits on zero, nothing is routed through M.
if (nrow(med_tbl)) {
  
  decomp <- bind_rows(
    med_tbl %>% transmute(group, mediator, outcome, n,
                          effect = "Total (c)",     est = c_total_std,
                          lo = c_lo,  hi = c_hi),
    med_tbl %>% transmute(group, mediator, outcome, n,
                          effect = "Direct (c')",   est = cprime_std,
                          lo = cp_lo, hi = cp_hi),
    med_tbl %>% transmute(group, mediator, outcome, n,
                          effect = "Indirect (ab)", est = ab_std,
                          lo = ab_lo, hi = ab_hi)
  ) %>%
    mutate(
      effect = factor(effect, levels = c("Indirect (ab)", "Direct (c')", "Total (c)")),
      cell   = paste0(group, "\n", vapply(mediator, med_lab, ""),
                      "  ->  ", vapply(outcome, med_lab, ""), "   (n = ", n, ")")
    )
  
  fp <- ggplot(decomp, aes(x = est, y = effect, colour = group)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.16, linewidth = 0.8) +
    geom_point(size = 3) +
    facet_wrap(~ cell, ncol = 2) +
    scale_colour_manual(values = if (exists("custom_colors")) custom_colors else NULL,
                        guide = "none") +
    labs(title = "Effect decomposition: total = direct + indirect",
         subtitle = sprintf(
           "%s -> receptor -> outcome. Standardized within group; bootstrap percentile 95%% CIs.",
           med_lab(MED_X)),
         x = "Effect (SD units)", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold", size = 9, lineheight = 1.1),
          panel.spacing = unit(1.1, "lines"))
  
  ggsave(file.path(MED_DIR, "Mediation_effect_decomposition.png"),
         fp, width = 10, height = 7, dpi = 300)
  
  ## Indirect effects alone, all cells on one axis (the original forest plot).
  fp2 <- med_tbl %>%
    mutate(cell = paste0(group, "\n", vapply(mediator, med_lab, ""),
                         "  ->  ", vapply(outcome, med_lab, ""))) %>%
    ggplot(aes(x = ab_std, y = cell, colour = group)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmin = ab_lo, xmax = ab_hi), height = 0.18, linewidth = 0.8) +
    geom_point(size = 3.2) +
    scale_colour_manual(values = if (exists("custom_colors")) custom_colors else NULL) +
    labs(title = "Indirect effects (ab) with bootstrap 95% CIs",
         subtitle = sprintf("%s -> receptor -> outcome, standardized within group",
                            med_lab(MED_X)),
         x = "Indirect effect (SD units)", y = NULL, colour = NULL) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
  
  ggsave(file.path(MED_DIR, "Mediation_indirect_effects_forest.png"),
         fp2, width = 9, height = 5.5, dpi = 300)
}


## ---- 15h. EXPORT -------------------------------------------------------------

write_xlsx(
  list(
    Single_mediator   = med_tbl,
    Moderated_med     = modmed_tbl,
    Parallel_mediator = parallel_tbl,
    Diagnostics       = diag_tbl,
    Data_used         = bind_rows(med_data)
  ),
  file.path(MED_DIR, "Stats_Mediation_ROE_receptor_NRS.xlsx")
)

message("\nMediation results, path diagrams and forest plot written to: ",
        normalizePath(MED_DIR))


## ---- 15i. OPTIONAL SENSITIVITY CHECKS ----------------------------------------
## Uncomment individually. Each one re-runs section 15 with a different setting;
## easiest is to change the config at 15a and re-source from there.
##
## 1) Adjust for sex (the main confound: 87% vs 52% female):
##      MED_COVARS <- c("female")
##
## 2) Use prescribed rather than measured exposure:
##      MED_X <- "log_MME"        # larger n, but MME is what was prescribed
##
## 3) Add pain quality as a second outcome:
##      MED_Y <- c("nrs", "PC2")
##
## 4) Restrict to the overlapping ROE range (addresses the limited-overlap
##    caveat raised in the lab meeting):
##      rng <- range(intersect(
##        range(master$log_ROE[master$Plot_Group == MED_GROUPS[1]], na.rm = TRUE),
##        range(master$log_ROE[master$Plot_Group == MED_GROUPS[2]], na.rm = TRUE)))
##      master_overlap <- master[is.na(master$log_ROE) |
##                               (master$log_ROE >= rng[1] & master$log_ROE <= rng[2]), ]
##    then swap `master` for `master_overlap` in 15c-15f.
################################################################################

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
## the packages already loaded (dplyr, purrr, tibble, ggplot2, writexl). No
## new packages required. Reuses `fmt_p` if it already exists in the session.
################################################################################

## ---- 16a. CONFIG -------------------------------------------------------------

MVR_X         <- "log_ROE"        # exposure predictor
MVR_RECEPTORS <- RECEPTORS        # full 19-receptor panel, defined in section 2
MVR_Y         <- "nrs"            # outcome; change to "PC2" etc. to re-run
MVR_GROUPS    <- FOCUS_GROUPS     # Hydrocodone group / Tramadol group
MVR_MIN_N     <- 12               # skip a group x receptor cell smaller than this

MVR_DIR <- file.path(OUT_DIR, "multivariable_regression")
dir.create(MVR_DIR, showWarnings = FALSE)

if (!exists("fmt_p")) {
  fmt_p <- function(p) ifelse(is.na(p), "NA",
                              ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

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
  mutate(p_Receptor_FDR = p.adjust(p_Receptor, method = "fdr"),
         p_ROE_FDR      = p.adjust(p_ROE,      method = "fdr")) %>%
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
  
  ok <- requireNamespace("ggrepel", quietly = TRUE)
  if (!ok) {
    message("ggrepel not installed -- installing for label placement in the dissociation scatter")
    install.packages("ggrepel")
  }
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
## Self-contained: needs `master`, `FOCUS_GROUPS`, `OUT_DIR` from earlier
## sections, plus lavaan and writexl.
################################################################################

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

if (!exists("MED_LABELS")) MED_LABELS <- character(0)
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

################################################################################
## 18. MULTIPLE REGRESSION: NRS ~ log ROE + 5-HT4 / 5-HT6
##     Three pre-specified models, fit independently in each drug group.
##
## MODELS (each fit separately within Hydrocodone group and Tramadol group)
##   M1  "HT4"      NRS = b0 + b_ROE*log_ROE + b_HT4*R_5HT4                + e
##   M2  "HT6"      NRS = b0 + b_ROE*log_ROE                + b_HT6*R_5HT6 + e
##   M3  "HT4_HT6"  NRS = b0 + b_ROE*log_ROE + b_HT4*R_5HT4 + b_HT6*R_5HT6 + e
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
##   M1 and M2 are numerically identical to section 16's R_5HT4 and R_5HT6
##   rows (same formula, same data); that is a deliberate cross-check.
##
## WHY M3 IS THE INTERESTING ONE
##   M1 and M2 each ask "does this receptor track pain once exposure is netted
##   out". M3 asks the harder question: do 5-HT4 and 5-HT6 carry SEPARATE
##   information about pain, or is one of them just standing in for the other?
##   If both coefficients survive in M3 at roughly their M1/M2 values, they are
##   independent signals. If one collapses toward zero, the two receptors were
##   sharing variance and only one of them is doing work. 18d tests this
##   formally with nested F tests (M1 vs M3, M2 vs M3).
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
##
## THE -Inf TRAP
##   log_ROE = log10(ROE) is -Inf for the 4 subjects per group with ROE == 0,
##   and complete.cases() does NOT drop -Inf. Every model frame below is
##   filtered with is.finite(), matching sections 16 and 17. Expect n = 19
##   (hydrocodone) and n = 17 (tramadol).
##
## Self-contained: needs `master`, `FOCUS_GROUPS`, `OUT_DIR`, `fmt_p`,
## optionally `custom_colors` and `MED_LABELS`, and the packages already
## loaded (dplyr, tidyr, tibble, purrr, ggplot2, writexl). No new packages --
## the 2x2 diagnostic panels use base graphics on a png() device precisely so
## that patchwork/cowplot are not needed.
################################################################################

## ---- 18a. CONFIG -------------------------------------------------------------

MREG_Y      <- "nrs"          # outcome
MREG_X      <- "log_ROE"      # exposure, present in every model
MREG_GROUPS <- FOCUS_GROUPS   # Hydrocodone group / Tramadol group
MREG_MIN_N  <- 12             # skip a group x model cell smaller than this

## The three models. Each entry is the full predictor vector, in the order the
## terms should appear. Change this list (and MREG_Y) and re-source from 18a to
## re-run the whole section on different predictors -- nothing downstream is
## hard-coded to 5-HT4 / 5-HT6 or to three models.
MREG_MODELS <- list(
  HT4     = c("log_ROE", "R_5HT4"),
  HT6     = c("log_ROE", "R_5HT6"),
  HT4_HT6 = c("log_ROE", "R_5HT4", "R_5HT6")
)

## Pretty names for the model columns/facets, in display order.
MREG_MODEL_LABELS <- c(
  HT4     = "M1: ROE + 5-HT4",
  HT6     = "M2: ROE + 5-HT6",
  HT4_HT6 = "M3: ROE + 5-HT4 + 5-HT6"
)

## Nested comparisons to run in 18d: c(reduced, full).
MREG_NESTED <- list(c("HT4", "HT4_HT6"), c("HT6", "HT4_HT6"))

MREG_DIR <- file.path(OUT_DIR, "multiple_regression_HT4_HT6")
dir.create(MREG_DIR, showWarnings = FALSE, recursive = TRUE)

MREG_COLS <- if (exists("custom_colors")) custom_colors else
  setNames(c("#1b9e77", "#d95f02")[seq_along(MREG_GROUPS)], MREG_GROUPS)

## Variable label lookup. Robust to keys that are not in MED_LABELS (unlike
## med_lab at section 15), same shape as section 17's pmed_lab.
mreg_lab <- function(v) {
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
  out[out == "log_ROE"] <- "log10 ROE"
  out
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

cat("\n================================================================\n")
cat("18. MULTIPLE REGRESSION: ", MREG_Y, " ~ log ROE + 5-HT4 / 5-HT6\n", sep = "")
cat("================================================================\n")

mreg_fits <- list()
mreg_coef_rows <- list()
mreg_fit_rows  <- list()

for (g in MREG_GROUPS) {
  d_g <- master[!is.na(master$Plot_Group) & master$Plot_Group == g, , drop = FALSE]
  message("Group: ", g, "  (", nrow(d_g), " subjects before filtering)")
  mreg_fits[[g]] <- list()

  for (mk in names(MREG_MODELS)) {
    res <- mreg_fit_one(d_g, MREG_MODELS[[mk]], MREG_Y, mk, g)
    if (is.null(res)) next
    mreg_fits[[g]][[mk]]            <- res
    mreg_coef_rows[[paste(g, mk)]]  <- res$coef
    mreg_fit_rows[[paste(g, mk)]]   <- res$fit
  }
}

mreg_coef_tbl <- dplyr::bind_rows(mreg_coef_rows)
mreg_fit_tbl  <- dplyr::bind_rows(mreg_fit_rows)

if (!nrow(mreg_coef_tbl)) stop("Section 18: no models could be fit -- check MREG_MODELS and `master`.")

mreg_coef_tbl <- mreg_coef_tbl %>%
  dplyr::mutate(
    model_lab = unname(MREG_MODEL_LABELS[model]),
    term_lab  = mreg_short(term),
    sig       = !is.na(p) & p < 0.05,
    p_fmt     = mreg_fmt_p(p)
  )

mreg_fit_tbl <- mreg_fit_tbl %>%
  dplyr::mutate(model_lab = unname(MREG_MODEL_LABELS[model]),
                F_p_fmt   = mreg_fmt_p(F_p))

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
## Does the second receptor add anything over the first? M1 and M2 are each
## nested inside M3, so an F test on the residual sums of squares is the right
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

    mreg_nested_rows[[paste(g, pair[1])]] <- tibble::tibble(
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

mreg_save <- function(p, file, w, h) {
  ggplot2::ggsave(file.path(MREG_DIR, file), p, width = w, height = h,
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
## whether 5-HT4 / 5-HT6 hold their M1/M2 values once they are in M3 together.
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
    ggplot2::expand_limits(y = c(0.4, length(levels(cf$term_lab)) + 0.7)) +
    ggplot2::labs(
      title    = sprintf("Regression coefficients -- %s", g),
      subtitle = sprintf("%s predicted from log10 ROE and receptor activity; estimate with 95%% CI",
                         mreg_lab(MREG_Y)),
      x = sprintf("Change in %s per unit of predictor (raw units)", mreg_lab(MREG_Y)),
      y = NULL, caption = CAP_RAW
    ) +
    mreg_theme()

  mreg_save(p_cf, sprintf("MREG_coefficients_%s.png", mreg_slug(g)), w = 8.6, h = 8.2)
}

## (4) The regression table as a figure -- the "Table 2" of a regression paper,
## three models side by side. Drawn as a ggplot text grid rather than gt::gtsave,
## which would need chromote/webshot2 (not installed).
for (g in names(mreg_fits)) {
  mods <- names(mreg_fits[[g]])
  if (!length(mods)) next

  terms_order <- c("(Intercept)", unique(unlist(MREG_MODELS)))
  cf <- mreg_coef_tbl %>% dplyr::filter(group == g)

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

  ft <- mreg_fit_tbl %>% dplyr::filter(group == g)
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

  hdr <- data.frame(x = seq_len(ncol_tot), y = 0,
                    txt = c("", unname(MREG_MODEL_LABELS[mods])),
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
                         mreg_lab(MREG_Y), g),
      subtitle = "Raw-unit coefficient, standard error in parentheses",
      caption  = paste0("*** p < .001, ** p < .01, * p < .05.  ",
                        "Raw coefficients are not comparable across predictors (different scales); ",
                        "see the .xlsx for standardized betas.\n",
                        "VIF > 5 would indicate the two receptors are competing for the same variance.")
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

  mreg_save(p_tab, sprintf("MREG_table_%s.png", mreg_slug(g)),
            w = 2.6 * ncol_tot + 1.2, h = 0.38 * length(all_rows) + 2.4)
}


## ---- 18i. FIGURES ACROSS GROUPS ----------------------------------------------

## (5) Observed vs predicted. The most honest single picture of how much of NRS
## these models actually explain: points on the diagonal = good prediction,
## a flat cloud = the model knows nothing.
op_rows <- list()
for (g in names(mreg_fits)) for (mk in names(mreg_fits[[g]])) {
  res <- mreg_fits[[g]][[mk]]
  ft  <- mreg_fit_tbl[mreg_fit_tbl$group == g & mreg_fit_tbl$model == mk, ]
  op_rows[[paste(g, mk)]] <- tibble::tibble(
    group = g, model = mk,
    model_lab = unname(MREG_MODEL_LABELS[mk]),
    observed  = res$data[[res$yv]],
    predicted = unname(stats::fitted(res$model)),
    lab = sprintf("R2 = %.3f\nadj R2 = %.3f\nRMSE = %.2f\nn = %d",
                  ft$R2[1], ft$adj_R2[1], ft$RMSE[1], ft$n[1])
  )
}
op_df <- dplyr::bind_rows(op_rows)

if (nrow(op_df)) {
  op_df$model_lab <- factor(op_df$model_lab, levels = unname(MREG_MODEL_LABELS))
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

  mreg_save(p_op, "MREG_obs_vs_pred.png", w = 12.5, h = 8.6)
}

## (6) Model comparison: adjusted R2 and AIC side by side, with the nested F
## tests annotated. adj R2 and AIC both penalise the extra parameter, so if M3
## does not beat M1/M2 on these it is not earning its keep.
if (nrow(mreg_fit_tbl)) {
  comp <- mreg_fit_tbl %>%
    dplyr::select(group, model, model_lab, R2, adj_R2, AIC) %>%
    tidyr::pivot_longer(c(R2, adj_R2, AIC), names_to = "metric", values_to = "value") %>%
    dplyr::mutate(
      metric = factor(metric, levels = c("R2", "adj_R2", "AIC"),
                      labels = c("R2 (in-sample)", "Adjusted R2", "AIC (lower is better)")),
      model_lab = factor(model_lab, levels = unname(MREG_MODEL_LABELS)),
      group = factor(group, levels = MREG_GROUPS)
    )

  ## Nested-F results as the subtitle, ONE LINE PER GROUP -- all four on a
  ## single line runs off the right edge of the device and gets clipped.
  sub_txt <- NULL
  if (nrow(mreg_nested_tbl)) {
    per_group <- vapply(unique(mreg_nested_tbl$group), function(gg) {
      r <- mreg_nested_tbl[mreg_nested_tbl$group == gg, ]
      paste0(gg, ":  ",
             paste(sprintf("+%s over %s: F(%d,%d) = %.2f, p = %s",
                           mreg_short(r$added_term), r$reduced,
                           r$df_num, r$df_den, r$F_stat, mreg_fmt_p(r$p)),
                   collapse = "   |   "))
    }, character(1))
    sub_txt <- paste0("Nested F tests (does the second receptor add anything?)\n",
                      paste(per_group, collapse = "\n"))
  }

  p_comp <- ggplot2::ggplot(comp, ggplot2::aes(model_lab, value, fill = group)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.78), width = 0.72,
                      alpha = 0.9) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", value)),
                       position = ggplot2::position_dodge(width = 0.78),
                       vjust = -0.45, size = 3.1, colour = "grey25") +
    ggplot2::scale_fill_manual(values = MREG_COLS, name = NULL) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.14))) +
    ggplot2::scale_x_discrete(labels = function(x) sub(": ", ":\n", x)) +
    ggplot2::facet_wrap(~ metric, scales = "free_y", nrow = 1) +
    ggplot2::labs(
      title    = "Model comparison",
      subtitle = sub_txt,
      x = NULL, y = NULL,
      caption = paste0("R2 can only rise when a predictor is added; adjusted R2 and AIC ",
                       "penalise the extra parameter, so they are the ones to compare.")
    ) +
    mreg_theme() +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 8.5, colour = "grey30"),
                   axis.text.x = ggplot2::element_text(size = 9))

  mreg_save(p_comp, "MREG_model_comparison.png", w = 13.5, h = 7)
}

## (7) Marginal (adjusted) effect of each receptor from the full model M3:
## predicted NRS across the observed range of that receptor with the other
## predictors held at their group means, 95% CI ribbon, raw data overlaid as
## partial residuals so the points sit on the same scale as the line.
MREG_FULL <- names(MREG_MODELS)[length(MREG_MODELS)]   # M3 by construction

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
    ggplot2::expand_limits(y = c(0.4, length(levels(idf$term_lab)) + 0.7)) +
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

  mreg_save(p_gd, "MREG_group_difference_forest.png", w = 9, h = 8.2)
}


## ---- 18j. EXPORT -------------------------------------------------------------

mreg_config <- tibble::tibble(
  setting = c("outcome", "exposure", "groups", "models", "min_n", "scaling",
              "n_rule", "output_dir"),
  value = c(
    MREG_Y, MREG_X, paste(MREG_GROUPS, collapse = " | "),
    paste(sprintf("%s: %s ~ %s", names(MREG_MODELS), MREG_Y,
                  vapply(MREG_MODELS, paste, character(1), collapse = " + ")),
          collapse = "  |  "),
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


## ---- 18k. NOTES --------------------------------------------------------------
## - Read the added-variable plots before the coefficient table. A coefficient
##   whose AV panel is a shapeless cloud with one far-out point is a coefficient
##   driven by that point, and the Cook's distance panel in the diagnostics
##   figure will name the subject.
## - n = 19 / 17 with 3-4 parameters in M3. These are pre-specified receptors,
##   so no multiplicity correction is applied -- but the CIs are wide and any
##   result here needs replication before it is more than a lead. Section 16's
##   19-receptor screen is where FDR belongs; this is not that.
## - Cross-check: M1's R_5HT4 coefficient and M2's R_5HT6 coefficient are the
##   same quantity as section 16's raw b_Receptor_raw for those receptors, and
##   the same as section 15's b path. Any disagreement means one of the three
##   sections is filtering rows differently.
## - The 4 subjects per group with ROE == 0 are dropped, not imputed, exactly as
##   in sections 15-17. If you want them in, change the log transform at section
##   2 to log10(ROE + 0.5*min(ROE[ROE > 0])) and re-source; do it as a labelled
##   sensitivity analysis, because the offset is arbitrary.
## - To re-run on a different outcome or different receptors, edit MREG_Y and
##   MREG_MODELS at 18a (and MREG_MODEL_LABELS / MREG_NESTED to match) and
##   re-source from 18a. Nothing downstream is hard-coded to two receptors or
##   three models.
################################################################################
