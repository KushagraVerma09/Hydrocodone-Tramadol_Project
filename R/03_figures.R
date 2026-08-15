###############################################################################
## 03_figures.R
##
## Figure 1 / overlay engine, run once per receptor set (receptor_sets, from
## 01_load_data.R) x clinical outcome, plus the standalone HT6 t-test,
## receptor-group-means table, exposure-vs-outcome figures, and the stats-box
## placement audit. Restricted to Hydrocodone vs Tramadol throughout
## (FOCUS_GROUPS, from 00_config.R).
##
## The publication-formatted demographics table (section 10 in the original
## monolithic script) now lives in 02_table1.R instead -- it's Table 1's own
## output, not a figure, and keeping it there means running 02_table1.R alone
## produces the finished .xlsx AND .html/.docx without needing this script.
## Section numbers below jump from 9 to 11 for that reason, not a mistake.
##
## Run standalone: Rscript R/03_figures.R
## Sourced by: Run_All.R. Nothing later in the pipeline (04-07) depends on
## anything this script defines -- they only need master/config/01_load_data.
###############################################################################

.HTK_PROJECT_DIR <- "/Users/kushagraverma/Work/Projects/Hydrocodone+Tramadol_Project"
## Always re-sourced -- see the note above this line in 01_load_data.R.
source(file.path(.HTK_PROJECT_DIR, "R", "00_config.R"))
if (!exists("master"))              source(file.path(.HTK_PROJECT_DIR, "R", "01_load_data.R"))

################################################################################
## 7. FIGURES AND OVERLAYS  (Hydrocodone group vs Tramadol group)
################################################################################

## Axis titles are NEVER wrapped (see the note above labs() in the builders), so
## every label here has to fit on ONE line. They are kept short and to a single
## "NAME (description)" shape on purpose: a rotated y-axis title is bounded by the
## panel HEIGHT and ggplot clips rather than wraps it, and holding every label to
## one line is what keeps the left-hand axis furniture the same height from figure
## to figure -- which is what makes the margins line up across a composite.
OUTCOMES <- tribble(
  ~var,   ~label,                          ~tag,
  "nrs",  "Clinical Pain (NRS 0-10)",      "NRS",
  "PC1",  "PC1 (Functional Disability)",   "PC1",
  "PC2",  "PC2 (Pain Quality/Severity)",   "PC2",
  "PC3",  "PC3 (Negative Affect)",         "PC3",
  "SOWS", "SOWS (Withdrawal Symptoms)",    "SOWS"
)

PLOT_GROUPS <- FOCUS_GROUPS

## custom_colors, custom_shapes and STATS_MODES come from 00_config.R -- they
## are shared with 04-07 (which color figures the same way defensively) and
## with Compose_Figure_Panels.R (which reads the STATS_MODES folder names).

## ---- 7a. WHERE THE n / r / slope / p TEXT GOES --------------------------------
## Two placements, both produced on every run into parallel folders that hold
## identical filenames, so the same figure can be compared side by side:
##
##   "header" -> the stats become a second line on the strip above each panel
##               (and the legend key text on the overlay figures). The text is
##               outside the panel, so it cannot touch the data at all.
##   "corner" -> the stats stay inside the panel, in a box placed in whichever
##               corner is emptiest (see pick_corner below).
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
##
## The defaults MUST match what the geom_label() calls actually draw with, or the
## measurement is of a box that doesn't exist and the placement puts the real one
## on top of the data. That is why these are the shared 00_config.R constants and
## not the bare literals they used to be -- the two sides now cannot drift apart.
label_size_in <- function(txt, size_mm = FIG_ANNOT_SIZE,
                          lineheight = FIG_ANNOT_LINEHEIGHT,
                          pad_lines  = FIG_ANNOT_PAD) {
  fs <- size_mm * ggplot2::.pt
  g  <- grid::textGrob(txt, gp = grid::gpar(fontsize = fs,
                                            fontface = FIG_ANNOT_FACE,
                                            fontfamily = FIG_FONT,
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

## Facet strip carrying a bold group name above a plain stats line, both at the
## SAME size (FIG_STRIP_SIZE) -- matching the group name is the point, not an
## independent size, so there is nothing left here for FIG_STRIP_STATS to
## drive; that constant only sizes the overlay legend text now (see theme_fig).
##
## ggplot styles a strip with ONE element_text, so two weights in one strip
## normally needs ggtext. plotmath does it without the dependency: atop() stacks
## the lines and bold() weights the first. (An earlier version wrapped the
## second line in scriptstyle(), which is a FIXED ~0.71x plotmath scale --
## constant regardless of FIG_STRIP_STATS -- and is why that config value used
## to look like it did nothing.) The catch is that a lookup vector passed to
## labeller() is never parsed, so the panels are faceted on a column whose
## LEVELS are the plotmath strings and label_parsed is applied to those.
## `stats_txt` may contain newlines; each becomes another stacked atop() line.
strip_math <- function(name, stats_txt) {
  ## Quotes would end the plotmath string early; nothing else needs escaping.
  clean <- function(s) gsub('"', "", s, fixed = TRUE)
  vapply(seq_along(name), function(i) {
    lines <- strsplit(clean(stats_txt[i]), "\n", fixed = TRUE)[[1]]
    inner <- sprintf('plain("%s")', lines)
    body  <- inner[length(inner)]
    for (k in rev(seq_len(length(inner) - 1)))
      body <- sprintf("atop(%s, %s)", inner[k], body)
    sprintf('atop(bold("%s"), %s)', clean(name[i]), body)
  }, character(1))
}

## Add the plotmath strip column to `plot_data`, ready for
## facet_wrap(~ .strip, labeller = label_parsed).
## glab() relabels the DISPLAYED group name ("Hydrocodone group" -> "CBP+Hydrocodone",
## per GROUP_LABELS in 00_config.R); the lookup vector stays NAMED by the real
## group values, so the facet still keys on the data and only the strip text
## changes.
strip_facet_data <- function(plot_data, gvar, names_in_order, stats_txt) {
  labs_v <- setNames(strip_math(glab(names_in_order), stats_txt), names_in_order)
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

## TRUE if every stats line fits on one strip line. The stats line now renders
## at the same size as the group name (FIG_STRIP_SIZE, no scriptstyle shrink
## -- see strip_math()), so it must be measured at that same size or this
## under-estimates its width and lets text that will actually clip pass as
## "fits on one line".
stats_fit_one_line <- function(txt, fig_w, n_panels) {
  all(vapply(txt, text_width_in, numeric(1), FIG_STRIP_SIZE) <
        0.97 * strip_budget_in(fig_w, n_panels))
}

## Shared theme for every scatter figure. Pulling it out of the four builders is
## what keeps the sizes consistent between them -- they had drifted apart.
## "Floating" axes: no gridlines and no panel box, so the only rules on the figure
## are the two that carry information -- a bold bottom x-axis and left y-axis with
## outward-facing ticks. Everything the data doesn't need is gone.
##
## For a FACETED figure the axis.line theme element is switched off here and
## drawn per panel instead, via panel_axis_lines() (added by the caller as a
## geom layer). axis.line.x.bottom/axis.line.y.left are genuinely outer-edge-only
## elements in ggplot2 -- in a facet_wrap they draw on the left column and bottom
## row of the WHOLE layout, not on every panel -- so in a 1x2 figure the right
## ("CBP+T") panel used to get only its bottom rule and no left-hand one: a full
## "L" axis on the left panel, a bare "_" on the right. panel_axis_lines() gives
## every panel both strokes.
theme_fig <- function(faceted = TRUE, legend = "none") {
  th <- theme_minimal(base_size = FIG_BASE_SIZE, base_family = FIG_FONT) +
    theme(
      panel.grid         = element_blank(),   # no major, no minor
      panel.border       = element_blank(),   # no top/right box
      axis.line.x.bottom = if (faceted) element_blank() else
        element_line(colour = "black", linewidth = FIG_AXIS_LINE, lineend = "square"),
      axis.line.y.left   = if (faceted) element_blank() else
        element_line(colour = "black", linewidth = FIG_AXIS_LINE, lineend = "square"),
      axis.ticks         = element_line(colour = "black", linewidth = FIG_TICK_LINE),
      axis.ticks.length  = FIG_TICK_LEN,      # positive length => points outward
      axis.title         = element_text(size = FIG_AXIS_TITLE, colour = "black"),
      axis.text          = element_text(size = FIG_AXIS_TEXT, colour = FIG_AXIS_TEXT_COL,
                                        face = "bold"),
      legend.position    = legend,
      legend.text        = element_text(size = FIG_STRIP_STATS),
      ## A bottom legend otherwise sits on a wide band of its own padding, which
      ## on the overlay figures was a visible strip of nothing between the x-axis
      ## title and the keys. Header mode overrides these again for its long keys.
      legend.box.spacing = unit(0.3, "lines"),
      legend.margin      = margin(0, 0, 0, 0),
      plot.margin        = margin(4, 8, 2, 4),
      ## Wide enough that the LAST x tick label of one panel and the FIRST of the
      ## next do not touch. At 1.1 lines they collided ("0.50" running straight
      ## into "-0.75") once the tick labels were scaled up -- the gap has to grow
      ## with the type, and this is the axis furniture most likely to be read
      ## wrong when it doesn't.
      panel.spacing      = unit(2.2, "lines")
    )
  if (faceted)
    th <- th + theme(
      ## Bold, because this strip carries the group name. Header-mode figures
      ## override it to "plain" -- there the label is plotmath and supplies its
      ## own bold() for the name and scriptstyle(plain()) for the stats line.
      strip.text = element_text(size = FIG_STRIP_SIZE, face = "bold",
                                lineheight = 0.95,
                                margin = margin(b = 3, t = 1))
    )
  th
}

## Draws a left + bottom axis rule (an "L") in EVERY facet panel, using -Inf/Inf
## endpoints so each one lands exactly on that panel's own rendered range --
## fixed-scale, free-scale, or expanded by coord_cartesian(ylim = ...) for the
## corner stats box, it doesn't matter which, the segment still lands on the
## edge. annotate() (rather than geom_segment + a dummy data frame) is what
## replicates the same segment into every panel without a spurious "aesthetics
## recycled" warning.
## Companion to theme_fig(faceted = TRUE), which blanks the normal (outer-edge-
## only) axis.line element so this doesn't draw a double-thickness line on the
## panel that used to get the real one.
panel_axis_lines <- function() {
  list(
    annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf,
             colour = "black", linewidth = FIG_AXIS_LINE, lineend = "square"),
    annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf,
             colour = "black", linewidth = FIG_AXIS_LINE, lineend = "square")
  )
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

## Axis titles used to be wrapped onto two lines by a wrap_lab() helper here.
## They are not any more: the labels in OUTCOMES and x_axis_label() were shortened
## instead, so every one fits on a single line and a title costs one line of
## furniture on every figure rather than one or two depending on its wording.
## That is what holds the panel size steady from figure to figure.

## p-value as it should read INSIDE a caption. Pasting fmt_p() after "p = " gives
## "p = <10^-5" for a very small p, which reads badly; this switches the operator.
fmt_p_cap <- function(p) {
  if (is.na(p)) return("p = --")
  if (p < 1e-5) return("p < 10^-5")
  paste0("p = ", fmt_p(p))
}

## The two group means in a caption span very different magnitudes depending on
## the figure -- log ROE sits around -2.8, raw receptor activity around 0.001 --
## so a fixed number of decimals either wastes width or prints "0.00". Precision
## is chosen from the larger of the pair and applied to BOTH, so the two numbers
## are always directly comparable.
fmt_mean_pair <- function(a, b) {
  m  <- suppressWarnings(max(abs(c(a, b)), na.rm = TRUE))
  dp <- if (is.finite(m) && m >= 1) 2 else 4
  sprintf(paste0("%.", dp, "f"), c(a, b))
}

## The caption is one horizontal line, in " | "-separated fields, to spend as
## little vertical space as possible. On the narrower overlay figures a long one
## can still overrun the width, so measure it (text_width_in, above) and fold it
## once at a field boundary only when it actually does not fit -- a fold at " | "
## keeps whole statistics together, which an ordinary text wrap would not.
fit_caption <- function(txt, fig_w, size_pt = FIG_CAPTION_SIZE) {
  if (is.na(txt) || !nzchar(txt)) return(txt)
  ## 0.35 in covers plot.margin's left+right (12pt = 0.17 in) with slack. Measured
  ## against the png device this is accurate to better than 0.02 in, so a caption
  ## that passes here is genuinely not clipped.
  if (text_width_in(txt, size_pt) <= fig_w - 0.35) return(txt)
  parts <- strsplit(txt, " | ", fixed = TRUE)[[1]]
  if (length(parts) < 3) return(txt)
  paste(paste(parts[1], parts[2], sep = " | "),
        paste(parts[-(1:2)], collapse = " | "), sep = "\n")
}

## Point size at which the longest overlay legend key fits the figure width.
##
## The overlay legend puts the whole stats line into the key label
## ("CBP+Hydrocodone: n = 25, r = 0.25, slope = 1.69 (p = 0.232)"), which on the
## 11 in overlay is close to the full width. ggplot CLIPS a key that overruns --
## it does not wrap or shrink it -- so the trailing "(p = ...)" silently
## disappears off the right edge, which is the one part of the label a reader is
## most likely to be looking for. Left-justifying (below) buys some room but not
## enough once the type is scaled up, so measure the label and step the size down
## only when it genuinely does not fit.
##
## The 1.0 in budget covers plot.margin left+right plus the key glyph and its
## padding; the 9 pt floor stops a pathological label from shrinking the legend
## into illegibility -- at that point the label itself is the thing to shorten.
fit_legend_size <- function(txt, fig_w, size_pt = FIG_STRIP_STATS) {
  if (!length(txt)) return(size_pt)
  widest <- max(vapply(txt, text_width_in, numeric(1), size_pt))
  budget <- fig_w - 1.0
  if (widest <= budget) return(size_pt)
  max(size_pt * budget / widest, 9)
}

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
  
  ## One line, " | "-separated. The SDs, the per-group n, the t degrees of
  ## freedom, the significance stars and Cohen's d have all come out of the TEXT:
  ## n is already on the facet strip / in the corner box, and every one of them --
  ## Cohen's d included -- is still carried in `row` and written to the
  ## Xaxis_ttest sheet. Nothing is lost, it just isn't spent on figure height.
  ## Same render-time naming as the axis title (receptor_set_label, 01) so the
  ## box and the axis under it never disagree: "5HT4", not "HT4_Only".
  xtag  <- if (is.na(rs_name)) "Receptor activity" else receptor_set_label(rs_name)
  mns   <- fmt_mean_pair(mean(x1), mean(x2))
  label <- sprintf("%s mean: %s = %s, %s = %s | Welch t = %.2f, %s",
                   xtag, glab(g1), mns[1], glab(g2), mns[2],
                   unname(tt$statistic), fmt_p_cap(tt$p.value))
  
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
    panel_axis_lines() +
    healthy_line(ref_x) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "grey30",
                linetype = "solid", linewidth = FIG_LINE_SIZE,
                fill = "grey80", alpha = 0.4) +
    geom_point(size = FIG_POINT_SIZE, alpha = 0.85) +
    scale_color_manual(values = custom_colors, labels = GROUP_LABELS) +
    scale_shape_manual(values = custom_shapes) +
    ## Titles are NOT wrapped: OUTCOMES and x_axis_label() are short enough to fit
    ## on one line, which is what keeps the axis furniture -- and so the panel --
    ## the same size on every figure. No caption: the n / r / slope / p and t-test
    ## stats stay in the in-panel box / strip (stats_mode), nothing below the plot.
    labs(x = x_axis_label(rs_name), y = ylab) +
    theme_fig(faceted = TRUE)

  if (stats_mode == "header") {
    ## Stats live on the strip, so the panel holds nothing but data. The label is
    ## plotmath and carries its own weights, so the strip element must be plain.
    p <- p + facet_wrap(~ .strip, scales = "fixed", labeller = label_parsed) +
      theme(strip.text = element_text(face = "plain", size = FIG_STRIP_SIZE,
                                      lineheight = 0.95,
                                      margin = margin(b = 6, t = 2)))
  } else {
    ## Measure the box against the faceted panel it will actually be drawn in,
    ## then place it. ref_x is passed so it does not land on the "Healthy mean"
    ## guide either.
    bf <- box_fracs(p + facet_wrap(~ Plot_Group, scales = "fixed",
                             labeller = as_labeller(GROUP_LABELS)),
                    FIG1_W, FIG1_H, summary_data$label_text)
    lp <- place_corner_labels(plot_data, "Plot_Group", "X_disease", yvar,
                              ref_x = if (SHOW_HEALTHY_LINE) ref_x else NA_real_,
                              w_frac = bf[["w"]], h_frac = bf[["h"]])
    audit_corner(if (is.null(file_stub)) paste0("Figure1_", tag, "_", rs_name)
                 else file_stub, lp)
    lab_df <- left_join(lp$coords, summary_data, by = "Plot_Group")
    p <- p +
      facet_wrap(~ Plot_Group, scales = "fixed",
                 labeller = as_labeller(GROUP_LABELS)) +
      geom_label(data = lab_df, aes(x = lab_x, y = lab_y, label = label_text),
                 inherit.aes = FALSE,
                 hjust = lab_df$hjust, vjust = lab_df$vjust,
                 size = FIG_ANNOT_SIZE, lineheight = FIG_ANNOT_LINEHEIGHT,
                   family = FIG_FONT,
                 fontface = FIG_ANNOT_FACE, color = FIG_ANNOT_COL, fill = "white", alpha = 1,
                 linewidth = 0, label.padding = unit(FIG_ANNOT_PAD, "lines")) +
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
    labs(x = x_axis_label(rs_name), y = ylab, color = NULL, shape = NULL) +
    theme_fig(faceted = FALSE, legend = "bottom")

  if (stats_mode == "header") {
    ## An overlay has no strip, so the legend key text carries the stats. Still
    ## outside the panel, so still incapable of covering a point.
    key_labs <- setNames(
      paste0(glab(summary_data$Plot_Group), ": ",
             stats_line(summary_data$n, summary_data$r,
                        summary_data$slope, summary_data$p)),
      as.character(summary_data$Plot_Group))
    p <- p +
      scale_color_manual(values = custom_colors, labels = key_labs) +
      guides(color = guide_legend(ncol = 1)) +
      ## One key per row, left-aligned and plain: a centred key this long
      ## overruns the device and ggplot clips it rather than wrapping.
      ## legend.location = "plot" justifies against the whole figure rather than
      ## against the panel, so the keys start at the far left and get the full
      ## width to run into -- without it the block stays centred and loses about
      ## 2 in of usable width to the leading offset before anything is measured.
      theme(legend.justification = "left",
            legend.location      = "plot",
            legend.text          = element_text(size = fit_legend_size(key_labs, OVL_W),
                                               face = "plain"),
            legend.margin        = margin(0, 0, 0, 0),
            legend.box.margin    = margin(0, 0, 0, 0),
            legend.key.spacing.y = unit(2, "pt"))
  } else {
    ## One shared box for both groups. The corner search runs on the pooled data
    ## because every point in the figure is drawn in this single panel, so
    ## "emptiest" has to be judged against all of them at once.
    box <- paste(sprintf("%s: %s", glab(summary_data$Plot_Group),
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
      scale_color_manual(values = custom_colors, labels = GROUP_LABELS) +
      annotate("label", x = lp$coords$lab_x[1], y = lp$coords$lab_y[1],
               hjust = lp$coords$hjust[1], vjust = lp$coords$vjust[1],
               label = box, size = FIG_ANNOT_SIZE, lineheight = FIG_ANNOT_LINEHEIGHT,
                   family = FIG_FONT,
               fontface = FIG_ANNOT_FACE, color = FIG_ANNOT_COL, fill = "white",
               linewidth = 0, label.padding = unit(FIG_ANNOT_PAD, "lines")) +
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
## 11. STANDALONE T-TEST: HT6 X VALUES (Hydrocodone vs Tramadol)
################################################################################

master_ht6 <- master
## Reads the HT6 receptor column(s) from receptor_sets (01_load_data.R)
## rather than hardcoding "R_5HT6" a second time, so a rename/remap of that
## receptor set is picked up here automatically instead of silently going
## stale.
master_ht6$X_HT6 <- compute_x_value(master_ht6, receptor_sets[["HT6_Only"]])

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

## Derived from OUTCOMES rather than restating the labels, so the y-axis title for
## a given outcome is character-for-character the same here as on the section 7/9
## figures. Restating them is how the two lists drifted before, and identical
## titles are what make the axis furniture -- and so the panel size and margins --
## line up when these figures are composed onto one page.
XY_OUTCOMES <- OUTCOMES %>% filter(var %in% c("nrs", "PC1", "PC2", "PC3"))

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
    
    ## ONE line, " | "-separated, so the caption costs a single line of figure
    ## height instead of four. fit_caption() folds it once at a " | " boundary
    ## only if it would actually overrun the device width.
    mns <- fmt_mean_pair(slopes$mean_x[1], slopes$mean_x[2])
    cap <- sprintf(
      "%s mean: %s = %s, %s = %s | Welch t = %.2f, %s | Slope diff: F = %.2f, %s",
      xtag, glab(levels(d$Plot_Group)[1]), mns[1],
            glab(levels(d$Plot_Group)[2]), mns[2],
      unname(tt_x$statistic), fmt_p_cap(tt_x$p.value),
      av$F[2], fmt_p_cap(av$`Pr(>F)`[2]))
    
    ## ---- figures, one pair per stats placement --------------------------------
    ## No title and no subtitle: the axis labels and the strip labels already say
    ## which predictor, which outcome and which group each panel shows.
    for (sm in names(STATS_MODES)) {

      ## ---- faceted figure ----------------------------------------------------
      p_facet <- ggplot(d, aes(.data[[xv]], .data[[yv]],
                               color = Plot_Group, shape = Plot_Group)) +
        panel_axis_lines() +
        geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "grey30",
                    linetype = "solid", linewidth = FIG_LINE_SIZE,
                    fill = "grey80", alpha = 0.4) +
        geom_point(size = FIG_POINT_SIZE, alpha = 0.85) +
        scale_color_manual(values = custom_colors, labels = GROUP_LABELS) +
        scale_shape_manual(values = custom_shapes) +
        ## No caption: stats stay in the in-panel box / strip, nothing below the plot.
        labs(x = xlab, y = yl) +
        theme_fig(faceted = TRUE)

      if (sm == "header") {
        ## Bold group name over a plain stats line, same size; the weights come
        ## from the plotmath label, so the strip element itself is plain.
        p_facet <- p_facet + facet_wrap(~ .strip, labeller = label_parsed) +
          theme(strip.text = element_text(face = "plain", size = FIG_STRIP_SIZE,
                                          lineheight = 0.95,
                                          margin = margin(b = 6, t = 2)))
      } else {
        bf     <- box_fracs(p_facet + facet_wrap(~ Plot_Group,
                                       labeller = as_labeller(GROUP_LABELS)), XY_W, XY_H,
                            slopes$label_text)
        lp     <- place_corner_labels(d, "Plot_Group", xv, yv,
                                      w_frac = bf[["w"]], h_frac = bf[["h"]])
        audit_corner(paste0("Fig_", combo, "_faceted"), lp)
        lab_df <- left_join(lp$coords, slopes, by = "Plot_Group")
        p_facet <- p_facet +
          facet_wrap(~ Plot_Group, labeller = as_labeller(GROUP_LABELS)) +
          geom_label(data = lab_df, aes(x = lab_x, y = lab_y, label = label_text),
                     inherit.aes = FALSE,
                     hjust = lab_df$hjust, vjust = lab_df$vjust,
                     size = FIG_ANNOT_SIZE, lineheight = FIG_ANNOT_LINEHEIGHT,
                   family = FIG_FONT,
                     fontface = FIG_ANNOT_FACE, color = FIG_ANNOT_COL, fill = "white", alpha = 1,
                     linewidth = 0, label.padding = unit(FIG_ANNOT_PAD, "lines")) +
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
        labs(x = xlab, y = yl, color = NULL, shape = NULL) +
        theme_fig(faceted = FALSE, legend = "bottom")

      if (sm == "header") {
        key_labs <- setNames(paste0(glab(slopes$Plot_Group), ": ",
                                    slopes$stat_line),
                             as.character(slopes$Plot_Group))
        p_overlay <- p_overlay +
          scale_color_manual(values = custom_colors, labels = key_labs) +
          guides(color = guide_legend(ncol = 1)) +
          ## See make_overlay: a centred key this long gets clipped, not wrapped.
          theme(legend.justification = "left",
                legend.location      = "plot",
                legend.text          = element_text(size = fit_legend_size(key_labs, XYO_W),
                                                   face = "plain"),
                legend.margin        = margin(0, 0, 0, 0),
                legend.box.margin    = margin(0, 0, 0, 0),
                legend.key.spacing.y = unit(2, "pt"))
      } else {
        box <- paste(sprintf("%s: %s", glab(slopes$Plot_Group),
                             slopes$stat_line), collapse = "\n")
        bf  <- box_fracs(p_overlay, XYO_W, XYO_H, box)
        lp  <- place_corner_labels(d %>% mutate(.all = factor("all")), ".all",
                                   xv, yv, extent_var = "Plot_Group",
                                   w_frac = bf[["w"]], h_frac = bf[["h"]])
        audit_corner(paste0("Fig_", combo, "_overlay"), lp)
        p_overlay <- p_overlay +
          scale_color_manual(values = custom_colors, labels = GROUP_LABELS) +
          annotate("label", x = lp$coords$lab_x[1], y = lp$coords$lab_y[1],
                   hjust = lp$coords$hjust[1], vjust = lp$coords$vjust[1],
                   label = box, size = FIG_ANNOT_SIZE, lineheight = FIG_ANNOT_LINEHEIGHT,
                   family = FIG_FONT,
                   fontface = FIG_ANNOT_FACE, color = FIG_ANNOT_COL, fill = "white",
                   linewidth = 0, label.padding = unit(FIG_ANNOT_PAD, "lines")) +
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


################################################################################
## 14a. THE TWO EXPOSURE MEASURES AGAINST EACH OTHER: MME (x) vs ROE (y)
##      — Hydrocodone group vs Tramadol group
##
## Section 14 puts MME and ROE on the x axis SEPARATELY, against the clinical
## outcomes. This asks the prior question: within a drug group, do the two
## exposure measures track each other at all? A high correlation means they are
## near-interchangeable and section 14's two rows of figures are telling one
## story twice; a low one means they carry different information.
##
## BOTH SCALES ARE DRAWN, and they are not interchangeable:
##   log-log : matches every other exposure figure in this file, but log10(0)
##             is -Inf, so the subjects with ROE = 0 (4 per group) drop out.
##   raw     : keeps those subjects, at the cost of a very skewed x and y.
## Reporting only one would hide either 8 subjects or the skew, so both are
## produced and the n is on every panel.
##
## Same conventions as section 14 throughout -- shared helpers, both stats
## placements, faceted + overlay -- so these figures compose with the rest.
################################################################################

## Both members of a pair are exposure measures, so unlike section 14 the x is
## not "the predictor" in any causal sense: this is a correlation, and the
## slope is reported only because the stats helpers print it.
##
## show_slope: on the RAW scale the slope is ROE per mg/day of MME, i.e. about
## 1.7e-04, which the shared stats helpers print as "slope = 0.00" -- a number
## that says nothing and reads as a flat line. Those helpers' %.2f is calibrated
## for every other figure in this file and is not worth changing for one panel,
## so the raw figures print n / r / p instead and leave the slope to the
## workbook, where it has room for its real magnitude. On the log-log scale the
## slope IS meaningful at 2dp (it is an elasticity: a slope near 1 means ROE
## rises proportionally with MME), so there it stays.
MME_ROE_SCALES <- list(
  list(xv = "MME",     yv = "ROE",     tag = "MME_vs_ROE",
       xlab = "MME (mg/day)",
       ylab = "ROE (mg/L)",
       show_slope = FALSE),
  list(xv = "log_MME", yv = "log_ROE", tag = "logMME_vs_logROE",
       xlab = expression(log[10]~"MME (mg/day)"),
       ylab = expression(log[10]~"ROE (mg/L)"),
       show_slope = TRUE)
)

## Slope-free counterparts of stats_block / stats_line / stats_strip, used when
## show_slope is FALSE. Same field order and same separators as the originals.
stats_block_nos <- function(n, r, p)
  sprintf("n = %d\nr = %.2f\np = %s", n, r, vapply(p, fmt_p, character(1)))
stats_line_nos <- function(n, r, p)
  sprintf("n = %d, r = %.2f, p = %s", n, r, vapply(p, fmt_p, character(1)))
stats_strip_nos <- function(n, r, p)
  sprintf("n = %d,  r = %.2f\np = %s", n, r, vapply(p, fmt_p, character(1)))

## Fisher r-to-z: is the correlation in one group different from the other?
## Two independent samples, so the z is the difference of the transformed r's
## over the SE of that difference. Needs n > 3 per group to be defined.
fisher_r_diff <- function(r1, n1, r2, n2) {
  if (any(!is.finite(c(r1, r2))) || n1 < 4 || n2 < 4)
    return(list(z = NA_real_, p = NA_real_))
  zf <- function(r) 0.5 * log((1 + r) / (1 - r))
  z  <- (zf(r1) - zf(r2)) / sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
  list(z = z, p = 2 * pnorm(-abs(z)))
}

mr_stats_all <- list()
mr_tests_all <- list()

for (sc in MME_ROE_SCALES) {

  xv <- sc$xv; yv <- sc$yv; combo <- sc$tag

  if (!all(c(xv, yv) %in% names(master))) {
    message(sprintf("Skipping %s — column not found in master", combo)); next
  }

  ## is.finite() on BOTH, not just !is.na(): under the log scale the ROE = 0
  ## subjects are -Inf rather than NA, and an -Inf would silently poison the
  ## correlation and the axis range alike.
  d <- master %>%
    filter(Plot_Group %in% FOCUS_GROUPS,
           is.finite(.data[[xv]]), is.finite(.data[[yv]])) %>%
    mutate(Plot_Group = droplevels(factor(Plot_Group, levels = FOCUS_GROUPS)))

  if (nlevels(d$Plot_Group) < 2 || any(table(d$Plot_Group) < 3)) {
    message(sprintf("Skipping %s — fewer than 3 subjects in a group", combo)); next
  }

  ## ---- per-group correlation + regression -------------------------------------
  mr_stats <- d %>%
    group_by(Plot_Group) %>%
    group_modify(~ {
      ct  <- cor.test(.x[[xv]], .x[[yv]])                  # Pearson, with CI
      sp  <- suppressWarnings(cor.test(.x[[xv]], .x[[yv]], method = "spearman"))
      fit <- lm(.x[[yv]] ~ .x[[xv]])
      tibble(n         = nrow(.x),
             r         = unname(ct$estimate),
             r_ci_lo   = ct$conf.int[1],
             r_ci_hi   = ct$conf.int[2],
             p         = ct$p.value,
             rho       = unname(sp$estimate),
             p_spearman = sp$p.value,
             slope     = coef(fit)[2],
             intercept = coef(fit)[1])
    }) %>% ungroup() %>%
    mutate(scale = combo,
           mean_x = map_dbl(Plot_Group, ~ mean(d[[xv]][d$Plot_Group == .x])),
           mean_y = map_dbl(Plot_Group, ~ mean(d[[yv]][d$Plot_Group == .x])),
           label_text = if (sc$show_slope) stats_block(n, r, slope, p)
                        else stats_block_nos(n, r, p),
           stat_line  = if (sc$show_slope) stats_line(n, r, slope, p)
                        else stats_line_nos(n, r, p),
           stat_strip = if (sc$show_slope) stats_strip(n, r, slope, p)
                        else stats_strip_nos(n, r, p))

  ## Strip labels settled before the plot is built (see make_figure1).
  d <- strip_facet_data(
    d, "Plot_Group", as.character(mr_stats$Plot_Group),
    if (stats_fit_one_line(mr_stats$stat_line, XY_W, 2))
      mr_stats$stat_line else mr_stats$stat_strip)

  ## ---- group comparisons --------------------------------------------------------
  ## Two different questions, both worth having: the interaction F asks whether
  ## the SLOPES differ (units of y per unit of x), Fisher's z asks whether the
  ## CORRELATIONS differ (how tightly the two measures track, unit-free). A pair
  ## of groups can easily differ on one and not the other.
  m_add <- lm(d[[yv]] ~ d[[xv]] + d$Plot_Group)
  m_int <- lm(d[[yv]] ~ d[[xv]] * d$Plot_Group)
  av    <- anova(m_add, m_int)
  fz    <- fisher_r_diff(mr_stats$r[1], mr_stats$n[1],
                         mr_stats$r[2], mr_stats$n[2])

  mr_tests <- tibble(
    scale     = combo,
    test      = c("Slope difference (x by Group interaction)",
                  "Correlation difference (Fisher r-to-z)"),
    statistic = c(av$F[2], fz$z),
    df        = c(av$Df[2], NA_real_),
    p_value   = c(av$`Pr(>F)`[2], fz$p)
  )

  mr_stats_all[[combo]] <- mr_stats
  mr_tests_all[[combo]] <- mr_tests

  cat(sprintf("\n===== %s: Hydrocodone vs Tramadol =====\n", combo))
  print(as.data.frame(mr_stats %>%
                        select(Plot_Group, n, r, r_ci_lo, r_ci_hi, p, rho, slope)),
        row.names = FALSE, digits = 4)
  cat("\n")
  print(as.data.frame(mr_tests %>% select(test, statistic, df, p_value)),
        row.names = FALSE, digits = 4)

  ## ---- figures, one pair per stats placement ------------------------------------
  for (sm in names(STATS_MODES)) {

    ## ---- faceted figure --------------------------------------------------------
    p_facet <- ggplot(d, aes(.data[[xv]], .data[[yv]],
                             color = Plot_Group, shape = Plot_Group)) +
      panel_axis_lines() +
      geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "grey30",
                  linetype = "solid", linewidth = FIG_LINE_SIZE,
                  fill = "grey80", alpha = 0.4) +
      geom_point(size = FIG_POINT_SIZE, alpha = 0.85) +
      scale_color_manual(values = custom_colors, labels = GROUP_LABELS) +
      scale_shape_manual(values = custom_shapes) +
      labs(x = sc$xlab, y = sc$ylab) +
      theme_fig(faceted = TRUE)

    if (sm == "header") {
      p_facet <- p_facet + facet_wrap(~ .strip, labeller = label_parsed) +
        theme(strip.text = element_text(face = "plain", size = FIG_STRIP_SIZE,
                                        lineheight = 0.95,
                                        margin = margin(b = 6, t = 2)))
    } else {
      bf     <- box_fracs(p_facet + facet_wrap(~ Plot_Group,
                                     labeller = as_labeller(GROUP_LABELS)),
                          XY_W, XY_H, mr_stats$label_text)
      lp     <- place_corner_labels(d, "Plot_Group", xv, yv,
                                    w_frac = bf[["w"]], h_frac = bf[["h"]])
      audit_corner(paste0("Fig_", combo, "_faceted"), lp)
      lab_df <- left_join(lp$coords, mr_stats, by = "Plot_Group")
      p_facet <- p_facet +
        facet_wrap(~ Plot_Group, labeller = as_labeller(GROUP_LABELS)) +
        geom_label(data = lab_df, aes(x = lab_x, y = lab_y, label = label_text),
                   inherit.aes = FALSE,
                   hjust = lab_df$hjust, vjust = lab_df$vjust,
                   size = FIG_ANNOT_SIZE, lineheight = FIG_ANNOT_LINEHEIGHT,
                   family = FIG_FONT,
                   fontface = FIG_ANNOT_FACE, color = FIG_ANNOT_COL, fill = "white",
                   alpha = 1, linewidth = 0,
                   label.padding = unit(FIG_ANNOT_PAD, "lines")) +
        coord_cartesian(ylim = lp$ylim)
    }

    ggsave(fig_path(sm, paste0("Fig_", combo, "_faceted")),
           p_facet, width = XY_W, height = XY_H, dpi = 300)

    ## ---- overlay figure --------------------------------------------------------
    p_overlay <- ggplot(d, aes(.data[[xv]], .data[[yv]],
                               color = Plot_Group, shape = Plot_Group)) +
      geom_point(alpha = 0.6, size = FIG_POINT_SIZE) +
      geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linetype = "solid",
                  alpha = 0.15, linewidth = FIG_LINE_SIZE) +
      scale_shape_manual(values = custom_shapes, guide = "none") +
      labs(x = sc$xlab, y = sc$ylab, color = NULL, shape = NULL) +
      theme_fig(faceted = FALSE, legend = "bottom")

    if (sm == "header") {
      key_labs <- setNames(paste0(glab(mr_stats$Plot_Group), ": ",
                                  mr_stats$stat_line),
                           as.character(mr_stats$Plot_Group))
      p_overlay <- p_overlay +
        scale_color_manual(values = custom_colors, labels = key_labs) +
        guides(color = guide_legend(ncol = 1)) +
        ## See make_overlay: a centred key this long gets clipped, not wrapped.
        theme(legend.justification = "left",
              legend.location      = "plot",
              legend.text          = element_text(size = fit_legend_size(key_labs, XYO_W),
                                                  face = "plain"),
              legend.margin        = margin(0, 0, 0, 0),
              legend.box.margin    = margin(0, 0, 0, 0),
              legend.key.spacing.y = unit(2, "pt"))
    } else {
      box <- paste(sprintf("%s: %s", glab(mr_stats$Plot_Group),
                           mr_stats$stat_line), collapse = "\n")
      bf  <- box_fracs(p_overlay, XYO_W, XYO_H, box)
      lp  <- place_corner_labels(d %>% mutate(.all = factor("all")), ".all",
                                 xv, yv, extent_var = "Plot_Group",
                                 w_frac = bf[["w"]], h_frac = bf[["h"]])
      audit_corner(paste0("Fig_", combo, "_overlay"), lp)
      p_overlay <- p_overlay +
        scale_color_manual(values = custom_colors, labels = GROUP_LABELS) +
        annotate("label", x = lp$coords$lab_x[1], y = lp$coords$lab_y[1],
                 hjust = lp$coords$hjust[1], vjust = lp$coords$vjust[1],
                 label = box, size = FIG_ANNOT_SIZE,
                 lineheight = FIG_ANNOT_LINEHEIGHT, family = FIG_FONT,
                 fontface = FIG_ANNOT_FACE, color = FIG_ANNOT_COL, fill = "white",
                 linewidth = 0, label.padding = unit(FIG_ANNOT_PAD, "lines")) +
        coord_cartesian(ylim = lp$ylim)
    }

    ggsave(fig_path(sm, paste0("Fig_", combo, "_overlay")),
           p_overlay, width = XYO_W, height = XYO_H, dpi = 300)
  }
}

## ---- combined export -----------------------------------------------------------
## No FDR here: two scales of the same two variables are not four independent
## questions, they are one question asked twice, so adjusting across them would
## be a penalty for reporting the sensitivity check rather than hiding it.
if (length(mr_stats_all)) {
  mr_stats_tbl <- bind_rows(mr_stats_all) %>%
    select(scale, Plot_Group, n, mean_x, mean_y, r, r_ci_lo, r_ci_hi, p,
           rho, p_spearman, slope, intercept)
  mr_tests_tbl <- bind_rows(mr_tests_all) %>%
    select(scale, test, statistic, df, p_value)

  write_xlsx(list(Per_group_correlation = mr_stats_tbl,
                  Group_comparisons     = mr_tests_tbl,
                  Data_used             = master %>%
                    filter(Plot_Group %in% FOCUS_GROUPS) %>%
                    select(PIN, Plot_Group, MME, log_MME, ROE, log_ROE)),
             file.path(OUT_DIR, "Stats_MME_vs_ROE.xlsx"))

  message("\nMME vs ROE figures and stats written to: ", normalizePath(OUT_DIR))
}

## ---- 14b. CORNER-PLACEMENT AUDIT ---------------------------------------------
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

