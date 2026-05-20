# Build subsampled tidy fixtures for cross-dataset integration tests.
# Skips datasets whose source CSV is absent.

library(readr)
library(dplyr)

src_paths <- list(
  yvo   = "/Users/dtl0018/Desktop/A_Proteomics_Analysis/A_YvO_2026/03_DEP/c_data/03_combined_results.csv",
  cvh   = c(
    "/Users/dtl0018/Desktop/A_Proteomics_Analysis/A_CvH_2026/03_DEP/c_data/03_combined_results_CR.csv",
    "/Users/dtl0018/Desktop/A_Proteomics_Analysis/A_CvH_2026/03_DEP/c_data/03_combined_results_CRvH.csv"
  ),
  hrvlr = "/Users/dtl0018/Desktop/A_Proteomics_Analysis/A_HRvLR_2026/03_DEP/c_data/03_combined_results.csv",
  mito  = "/Users/dtl0018/Desktop/A_Proteomics_Analysis/A_Mito_2026/03_DEP/c_data/03_combined_results.csv"
)

# Load enrichVolcano via devtools so ev_validate works
devtools::load_all(".")

dir.create("tests/testthat/fixtures", recursive = TRUE, showWarnings = FALSE)
set.seed(42)

for (ds in names(src_paths)) {
  paths <- src_paths[[ds]]
  existing <- paths[file.exists(paths)]
  if (length(existing) == 0) {
    message("Skipping ", ds, ": no source CSVs found.")
    next
  }
  wides <- lapply(existing, function(p) readr::read_csv(p, show_col_types = FALSE))
  wide_combined <- if (length(wides) > 1) {
    # Join by gene to merge wide-suffix tables that share gene IDs
    Reduce(function(a, b) {
      common_gene <- intersect(c("gene", "Gene", "name", "Protein"), intersect(colnames(a), colnames(b)))
      if (length(common_gene) > 0) {
        dplyr::full_join(a, b, by = common_gene[1])
      } else {
        dplyr::bind_rows(a, b)
      }
    }, wides)
  } else {
    wides[[1]]
  }
  tidy <- tryCatch(
    ev_validate(wide_combined),
    error = function(e) {
      message("Failed to validate ", ds, ": ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(tidy)) next
  # Subsample: keep top 200 proteins by |logFC| per contrast
  if ("contrast" %in% colnames(tidy)) {
    sub <- tidy |>
      dplyr::group_by(contrast) |>
      dplyr::arrange(dplyr::desc(abs(logFC))) |>
      dplyr::slice_head(n = 200) |>
      dplyr::ungroup()
  } else {
    sub <- tidy |>
      dplyr::arrange(dplyr::desc(abs(logFC))) |>
      dplyr::slice_head(n = 200)
  }
  out <- file.path("tests/testthat/fixtures", paste0(ds, "_tidy.rds"))
  saveRDS(sub, out)
  message("Wrote ", nrow(sub), " rows for ", ds, " to ", out)

  # Also bundle a smaller YvO subsample inside the package for vignette demos
  if (ds == "yvo") {
    bundle_dir <- "inst/extdata/examples"
    dir.create(bundle_dir, recursive = TRUE, showWarnings = FALSE)
    bundle <- if ("contrast" %in% colnames(tidy)) {
      tidy |>
        dplyr::group_by(contrast) |>
        dplyr::arrange(dplyr::desc(abs(logFC))) |>
        dplyr::slice_head(n = 100) |>
        dplyr::ungroup()
    } else {
      tidy |>
        dplyr::arrange(dplyr::desc(abs(logFC))) |>
        dplyr::slice_head(n = 100)
    }
    bundle_path <- file.path(bundle_dir, "yvo_tidy.rds")
    saveRDS(bundle, bundle_path)
    message("Bundled ", nrow(bundle), " rows for vignettes at ", bundle_path)
  }
}
