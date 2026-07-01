# Shape + type validators applied AFTER column resolution. Public-API
# boundary only; internal helpers trust the resolved frames.

#' @keywords internal
#' @noRd
validate_volc_df <- function(volc_df, cols, df_name = "volc_df") {
  if (!is.data.frame(volc_df)) {
    ev_abort(
      "{.arg {df_name}} must be a data.frame, not {.cls {class(volc_df)[1]}}.",
      class = "enrichVolcano_input_error"
    )
  }
  if (nrow(volc_df) == 0) {
    ev_abort(
      "{.arg {df_name}} has zero rows.",
      class = "enrichVolcano_input_error"
    )
  }
  if (!is.numeric(volc_df[[cols$logfc]])) {
    ev_abort(
      "Column {.val {cols$logfc}} of {.arg {df_name}} must be numeric.",
      class = "enrichVolcano_column_error"
    )
  }
  p <- volc_df[[cols$pval]]
  if (!is.numeric(p)) {
    ev_abort(
      "Column {.val {cols$pval}} of {.arg {df_name}} must be numeric.",
      class = "enrichVolcano_column_error"
    )
  }
  bad_p <- !is.na(p) & (p < 0 | p > 1)
  if (any(bad_p)) {
    ev_abort(
      c("Column {.val {cols$pval}} has values outside [0, 1].",
        "i" = "Found {sum(bad_p)} row{?s} with `p < 0` or `p > 1`."
      ),
      class = "enrichVolcano_data_error"
    )
  }
  invisible(volc_df)
}

#' @keywords internal
#' @noRd
validate_enrich_df <- function(enrich_df, cols, df_name = "enrich_df") {
  if (!is.data.frame(enrich_df)) {
    ev_abort(
      "{.arg {df_name}} must be a data.frame, not {.cls {class(enrich_df)[1]}}.",
      class = "enrichVolcano_input_error"
    )
  }
  if (nrow(enrich_df) == 0) {
    ev_abort(
      "{.arg {df_name}} has zero rows.",
      class = "enrichVolcano_input_error"
    )
  }
  nes <- enrich_df[[cols$nes]]
  if (!is.numeric(nes)) {
    ev_abort(
      "Column {.val {cols$nes}} of {.arg {df_name}} must be numeric.",
      class = "enrichVolcano_column_error"
    )
  }
  if (all(is.na(nes))) {
    ev_abort(
      c("Every {.val {cols$nes}} value is `NA`; cannot colour arcs.",
        "i" = "Check that your enrichment table has NES values."
      ),
      class = "enrichVolcano_data_error"
    )
  }
  padj <- enrich_df[[cols$padj]]
  if (!is.numeric(padj)) {
    ev_abort(
      "Column {.val {cols$padj}} of {.arg {df_name}} must be numeric.",
      class = "enrichVolcano_column_error"
    )
  }
  bad_padj <- !is.na(padj) & (padj < 0 | padj > 1)
  if (any(bad_padj)) {
    ev_abort(
      c("Column {.val {cols$padj}} has values outside [0, 1].",
        "i" = "Found {sum(bad_padj)} row{?s}."
      ),
      class = "enrichVolcano_data_error"
    )
  }
  dups <- duplicated(enrich_df[[cols$term]])
  if (any(dups)) {
    ev_warn(
      c("Duplicate values in {.val {cols$term}}; keeping the row with the lowest padj per term.",
        "i" = "{sum(dups)} duplicate row{?s} dropped."
      ),
      class = "enrichVolcano_input_error"
    )
  }
  invisible(enrich_df)
}

#' @keywords internal
#' @noRd
dedup_by_term <- function(enrich_df, cols) {
  if (!anyDuplicated(enrich_df[[cols$term]])) {
    return(enrich_df)
  }
  ord <- order(enrich_df[[cols$padj]], na.last = TRUE)
  enrich_df <- enrich_df[ord, , drop = FALSE]
  enrich_df[!duplicated(enrich_df[[cols$term]]), , drop = FALSE]
}

#' Validate the ring geometry knobs at the public boundary
#'
#' Catches mistuned `ring_radius` / `arc_height_range` before they silently
#' render an inverted or collapsed ring.
#' @keywords internal
#' @noRd
validate_ring_geometry <- function(ring_radius, volcano_radius, arc_height_range,
                                   label_headroom = 0.5) {
  if (!is.numeric(ring_radius) || length(ring_radius) != 1L ||
    !is.finite(ring_radius) || ring_radius <= 0) {
    ev_abort(
      "{.arg ring_radius} must be a single positive number; got {.val {ring_radius}}.",
      class = "enrichVolcano_param_error"
    )
  }
  if (!is.numeric(label_headroom) || length(label_headroom) != 1L ||
    !is.finite(label_headroom) || label_headroom < 0) {
    ev_abort(
      "{.arg label_headroom} must be a single non-negative number; got {.val {label_headroom}}.",
      class = "enrichVolcano_param_error"
    )
  }
  if (!is.numeric(arc_height_range) || length(arc_height_range) != 2L ||
    any(!is.finite(arc_height_range))) {
    ev_abort(
      "{.arg arc_height_range} must be a length-2 numeric {.code c(min, max)}.",
      class = "enrichVolcano_param_error"
    )
  }
  if (arc_height_range[1] < 0 || arc_height_range[2] < arc_height_range[1]) {
    ev_abort(
      c("{.arg arc_height_range} must be {.code c(min, max)} with `0 <= min <= max`.",
        "i" = "Got {.val {arc_height_range}}."
      ),
      class = "enrichVolcano_param_error"
    )
  }
  if (ring_radius < volcano_radius) {
    ev_warn(
      c(
        paste(
          "{.arg ring_radius} ({ring_radius}) is below {.arg volcano_radius}",
          "({volcano_radius}); the volcano will spill over the ring."
        ),
        "i" = "Raise {.arg ring_radius} or lower {.arg volcano_radius}."
      ),
      class = "enrichVolcano_param_warning"
    )
  }
  invisible(TRUE)
}

#' Validate the grid spacing knobs at the public boundary
#'
#' @keywords internal
#' @noRd
validate_grid_spacing <- function(panel_spacing, panel_margin, legend_width) {
  knobs <- list(
    panel_spacing = panel_spacing,
    panel_margin = panel_margin,
    legend_width = legend_width
  )
  for (nm in names(knobs)) {
    v <- knobs[[nm]]
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v < 0) {
      ev_abort(
        "{.arg {nm}} must be a single non-negative number (mm); got {.val {v}}.",
        class = "enrichVolcano_param_error"
      )
    }
  }
  invisible(TRUE)
}
