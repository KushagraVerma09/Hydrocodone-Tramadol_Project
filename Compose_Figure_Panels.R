###############################################################################
## COMPOSE MULTI-PANEL FIGURES  ->  outputs/composites/
##
## Takes PNGs that Hydrocodone_Tramadol_Kush_All_Graphs.R already wrote into
## ./New_Outputs/figs_stats_{header,corner}/ and glues them into single
## multi-panel figures laid out like the "Brain Opioid ALFF Figure N.docx"
## files: one composite image, panels tagged a, b, c, d ... in reading order,
## then an empty caption line underneath for you to type the figure legend into
## by hand.
##
## Every figure is built twice, once from each of the two stats-placement
## folders, into New_Outputs/composites/stats_header/ and .../stats_corner/.
##
## The original single-panel PNGs are NOT touched or overwritten.
##
## To change a figure: edit FIGURE_PANELS below. Add or remove a filename,
## change `ncol`, or add a whole new list(...) block. Nothing else needs to
## change. Any file listed that does not exist is skipped with a warning.
##
## Run AFTER the main script:
##   source("Compose_Figure_Panels.R")
###############################################################################


## ---- 0. PACKAGES -------------------------------------------------------------
pkgs <- c("magick", "officer")
missing <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(missing)) install.packages(missing)
invisible(lapply(pkgs, library, character.only = TRUE))

## ---- 1. PATHS ----------------------------------------------------------------
## DATA_DIR/OUT_DIR/STATS_MODES come from R/00_config.R -- the same shared
## config the pipeline scripts use -- so a change to the output folder name or
## the stats-placement folder names only has to be made in one place.
.HTK_PROJECT_DIR <- "/Users/kushagraverma/Work/Projects/Hydrocodone+Tramadol_Project"
source(file.path(.HTK_PROJECT_DIR, "R", "00_config.R"))

IN_DIR <- OUT_DIR                                         # where the PNGs live

## The main pipeline writes each scatter figure twice, once per placement of
## the n / r / slope / p text, into the two STATS_MODES folders. Filenames are
## identical in both, so one FIGURE_PANELS list below covers them and each
## gets its own composites folder. Local names here ("stats_header"/
## "stats_corner") only control this script's own output subfolder naming.
STATS_VARIANTS <- setNames(STATS_MODES, paste0("stats_", names(STATS_MODES)))

COMP_ROOT <- file.path(IN_DIR, "composites")             # where composites go

## ---- 2. LAYOUT SETTINGS ------------------------------------------------------
PANEL_WIDTH_PX <- 1600      # every panel is rescaled to this width before tiling
TAG_LABELS     <- letters   # panel tags: a, b, c, ...  (use LETTERS for A, B, C)
TAG_SIZE       <- 78        # point size of the panel tag
TAG_PAD_X_PX   <- 20        # white gutter on the left/right of each panel
## Height of the band added above each panel for its tag. Set smaller than
## TAG_SIZE on purpose: the figures carry a little white margin of their own
## above the strip label, so the glyph can drop into that and sit LEVEL with the
## top of the graph. A band taller than the glyph leaves the letter floating
## above its panel, which is what it used to do.
TAG_PAD_Y_PX   <- 46        # white band above each panel, holds the tag
TAG_OFFSET_PX  <- 4         # nudge the tag down from the very top edge
GAP_PX         <- 30        # white space between panels
DOC_WIDTH_IN   <- 6.5       # image width inside the .docx (US Letter, 1" margins)
OUT_DPI        <- 300       # dpi metadata stamped on the exported PNG/TIFF
WRITE_TIFF     <- TRUE      # the reference .docx files embed TIFF; set FALSE to skip

## ---- 3. THE FIGURES ----------------------------------------------------------
## name  : output file stem (composite PNG/TIFF/DOCX all take this name)
## title : the "Figure N." line placed above the caption slot in the .docx
## ncol  : how many panels per row
## files : panel images, in reading order, relative to outputs/
##         (swap "Fig1_NRS_HT4_Only.png" -> "Fig1_NRS_HT4_Only_overlay.png"
##          etc. if you want the overlay versions instead)
##
## This list is a curated selection, not derived from receptor_sets/OUTCOMES --
## which four panels belong on which manuscript figure is an editorial choice,
## not something that should auto-expand every time a receptor set is added.
## Filenames follow R/03_figures.R's own naming convention
## ("Fig1_<outcome tag>_<receptor set name>.png" for section 9,
## "Fig_<predictor>_vs_<outcome>_faceted.png" for section 14) -- if a receptor
## set or outcome tag is renamed there, the matching entries below need a
## matching (manual) rename, or they'll skip with a warning (see :188).

FIGURE_PANELS <- list(

  ## logMME no longer gets a composite of its own. Only its NRS panel carried
  ## the result worth showing, so that one panel rides along on the logROE
  ## figure below as panel e; the three logMME PC panels are dropped from the
  ## composites entirely (their single-panel PNGs are still written by the main
  ## script and untouched in figs_stats_{header,corner}/).
  list(
    name  = "Composite_Fig2_logROE_vs_outcomes",
    title = "Figure 2.",
    ncol  = 2,
    files = c("Fig_logROE_vs_NRS_faceted.png",
              "Fig_logROE_vs_PC1_faceted.png",
              "Fig_logROE_vs_PC2_faceted.png",
              "Fig_logROE_vs_PC3_faceted.png",
              ## panel e -- the ex-Figure 1 panel a
              "Fig_logMME_vs_NRS_faceted.png")
  ),

  list(
    name  = "Composite_Fig3_NRS_PC2_by_HT4_HT6",
    title = "Figure 3.",
    ncol  = 2,
    files = c("Fig1_NRS_HT4_Only.png",
              "Fig1_NRS_HT6_Only.png",
              "Fig1_PC2_HT4_Only.png",
              "Fig1_PC2_HT6_Only.png")
  ),

  ## CB1 across all four outcomes. Exploratory: CB1 is not in any of the
  ## mediation or multiple-regression models, so it gets its own figure rather
  ## than sharing one with the receptors that are.
  list(
    name  = "Composite_Fig4_CB1_vs_outcomes",
    title = "Figure 4.",
    ncol  = 2,
    files = c("Fig1_NRS_CB1_Only.png",
              "Fig1_PC1_CB1_Only.png",
              "Fig1_PC2_CB1_Only.png",
              "Fig1_PC3_CB1_Only.png")
  ),
  list(
    name  = "Composite_Fig5_Serotonine_vs_SOWS",
    title = "Figure 5.",
    ncol  = 2,
    files = c("Fig1_SOWS_HT1A_Only.png",
              "Fig1_SOWS_HT2A_Only.png",
              "Fig1_SOWS_HT4_Only.png",
              "Fig1_SOWS_HT6_Only.png")
  )
  
)

## ---- 4. COMPOSITION HELPERS --------------------------------------------------

## Load one panel, normalise its width, and stamp its tag letter in the
## top-left corner on a small white gutter so nothing in the plot is covered.
load_panel <- function(path, tag) {
  img <- image_read(path)
  img <- image_scale(img, paste0(PANEL_WIDTH_PX, "x"))

  ## Sides first. image_border's "WxH" geometry is symmetric, so it cannot give
  ## a band on the top alone -- and a matching band on the BOTTOM is wasted
  ## space that shows up as a wide gutter between rows of the composite.
  img <- image_border(img, color = "white", geometry = paste0(TAG_PAD_X_PX, "x0"))

  ## So grow the canvas upward instead: extent to a taller size anchored south,
  ## which puts all the new white space above the plot for the tag to sit in.
  inf <- image_info(img)
  img <- image_extent(img,
                      geometry = paste0(inf$width, "x", inf$height + TAG_PAD_Y_PX),
                      gravity = "south", color = "white")

  image_annotate(img, tag,
                 gravity  = "northwest",
                 location = paste0("+8+", TAG_OFFSET_PX),
                 size     = TAG_SIZE,
                 weight   = 700,
                 color    = "black")
}

## Pad every image in a set out to a common width/height on a white ground,
## anchored top-left, so rows and columns line up.
pad_to <- function(imgs, w, h) {
  image_join(lapply(imgs, function(im)
    image_extent(im, paste0(w, "x", h), gravity = "northwest", color = "white")))
}

## Lay panels out into `ncol` columns and return one stitched magick image.
tile_panels <- function(imgs, ncol) {
  n    <- length(imgs)
  rows <- split(seq_len(n), ceiling(seq_len(n) / ncol))

  row_imgs <- lapply(rows, function(idx) {
    set <- image_join(imgs[idx])
    inf <- image_info(set)
    r   <- image_append(pad_to(as.list(set), max(inf$width), max(inf$height)))
    ## trailing spacer keeps a short final row left-aligned rather than centred
    image_border(r, "white", paste0("0x", GAP_PX %/% 2))
  })

  stacked <- image_join(row_imgs)
  inf     <- image_info(stacked)
  image_append(pad_to(as.list(stacked), max(inf$width), max(inf$height)),
               stack = TRUE)
}

## ---- 5. BUILD EACH FIGURE, ONCE PER STATS VARIANT ----------------------------
for (variant in names(STATS_VARIANTS)) {

  SRC_DIR  <- file.path(IN_DIR, STATS_VARIANTS[[variant]])
  COMP_DIR <- file.path(COMP_ROOT, variant)

  if (!dir.exists(SRC_DIR)) {
    warning("No such panel folder: ", SRC_DIR,
            " - run the main script first. Skipping ", variant, ".", call. = FALSE)
    next
  }
  dir.create(COMP_DIR, showWarnings = FALSE, recursive = TRUE)

  message("\n=== ", variant, " -> ", COMP_DIR, " ===")

  for (fig in FIGURE_PANELS) {

    paths   <- file.path(SRC_DIR, fig$files)
    present <- file.exists(paths)

    if (any(!present)) {
      warning("Skipping missing panel(s) for ", fig$name, ": ",
              paste(fig$files[!present], collapse = ", "), call. = FALSE)
    }
    paths <- paths[present]
    if (!length(paths)) {
      warning("No panels found for ", fig$name, " - nothing written.", call. = FALSE)
      next
    }

    panels <- Map(load_panel, paths, TAG_LABELS[seq_along(paths)])
    sheet  <- tile_panels(image_join(panels), fig$ncol)
    sheet  <- image_border(sheet, "white", paste0(GAP_PX, "x", GAP_PX))

    ## -- raster exports --------------------------------------------------------
    png_path <- file.path(COMP_DIR, paste0(fig$name, ".png"))
    image_write(sheet, png_path, format = "png", density = OUT_DPI)

    if (WRITE_TIFF) {
      image_write(sheet,
                  file.path(COMP_DIR, paste0(fig$name, ".tiff")),
                  format = "tiff", density = OUT_DPI, compression = "LZW")
    }

    ## -- .docx: image sized to the text column, caption slot underneath --------
    info   <- image_info(sheet)
    doc_h  <- DOC_WIDTH_IN * info$height / info$width

    doc <- read_docx()
    doc <- body_add_img(doc, src = png_path,
                        width = DOC_WIDTH_IN, height = doc_h)
    doc <- body_add_par(doc, "")
    doc <- body_add_fpar(doc, fpar(ftext(fig$title, prop = fp_text(bold = TRUE))))
    doc <- body_add_par(doc, "")   # <- type the figure legend here
    print(doc, target = file.path(COMP_DIR, paste0(fig$name, ".docx")))

    message("Wrote ", fig$name, "  (", length(paths), " panels, ",
            info$width, "x", info$height, " px)")
  }
}

message("\nComposite figures written to: ", COMP_ROOT)