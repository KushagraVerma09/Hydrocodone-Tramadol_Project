###############################################################################
## 01_load_data.R
##
## Loads all five spreadsheets, harmonises them, and merges them into ONE
## master dataset (`master`). Also defines the receptor vocabulary
## (RECEPTORS_RAW/RECEPTORS) and the RECEPTOR SETS used by every figure and
## test downstream (`receptor_sets`, just below) -- add, remove or rename a
## receptor set there and every script that loops over receptor sets
## (03_figures.R) picks it up automatically.
##
## Run standalone: Rscript R/01_load_data.R  (writes the merged dataset to
## New_Outputs/ and stops -- makes no figures or tests).
## Sourced by: every other 0N_*.R script and Run_All.R.
###############################################################################

.HTK_PROJECT_DIR <- "/Users/kushagraverma/Work/Projects/Hydrocodone+Tramadol_Project"
## Always re-sourced (not gated on .HTK_CONFIG_LOADED): 00_config.R now
## re-evaluates its own values on every source() and only guards its one
## expensive step (packages) internally, so this stays cheap even when several
## 0N_*.R scripts each source it in the same R session -- see 00_config.R.
source(file.path(.HTK_PROJECT_DIR, "R", "00_config.R"))

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

## Receptor-set key -> the receptor name(s) as they should READ on a figure.
##
## Two things happen here. Every serotonin receptor gets its "5" back
## ("HT4" -> "5HT4"), because "HT4" on its own is not the receptor's name. And
## the structural words in the key become symbols ("HT6_plus_HT4" ->
## "5HT6 + 5HT4"), so the label is the receptor list and nothing else.
##
## This is a RENDER-time map only, exactly like GROUP_LABELS in 00_config.R:
## receptor_sets below, every output filename, and every sheet/column name stay
## keyed on the raw names, so nothing downstream can break on it.
receptor_set_label <- function(rs_name) {
  ## Sets that are a count rather than a list of receptors -- naming all 15, or
  ## all six, on an axis would not fit and would not be read.
  special <- c(All_Receptors = "All Receptors",
               Six_Specified = "Six Receptors")
  if (rs_name %in% names(special)) return(unname(special[rs_name]))

  toks <- strsplit(rs_name, "_", fixed = TRUE)[[1]]
  toks <- toks[!toks %in% c("Only", "plus", "and")]  # joiners, not receptors
  ## "HT1A" -> "5HT1A"; an already-prefixed "5HT1A" is left alone.
  toks <- sub("^5?HT", "5HT", toks)
  paste(toks, collapse = " + ")
}

## Human-readable x-axis title for a given receptor set: "<receptors> Activity",
## e.g. "5HT6 + 5HT4 Activity", "MOR Activity", "All Receptors Activity".
##
## Deliberately SHORT: 03_figures.R no longer wraps axis titles, so this has to
## fit on one line at FIG_AXIS_TITLE.
x_axis_label <- function(rs_name) {
  pretty <- receptor_set_label(rs_name)
  if (X_MODE == "raw")
    paste0(pretty, " Activity")
  else
    paste0("Disease Axis (0 = Healthy, 1 = CBP-O): ", pretty)
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
