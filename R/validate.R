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
    cli::cli_warn(
      c(
        paste(
          "{.arg ring_radius} ({ring_radius}) is below {.arg volcano_radius}",
          "({volcano_radius}); the volcano will spill over the ring."
        ),
        "i" = "Raise {.arg ring_radius} or lower {.arg volcano_radius}."
      ),
      class = c("enrichVolcano_param_warning", "enrichVolcano_warning")
    )
  }
  invisible(TRUE)
}
