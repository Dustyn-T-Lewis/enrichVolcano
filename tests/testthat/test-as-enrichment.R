source(test_path("fixtures/make_toy.R"))

read_fixture <- function(name, ...) {
  utils::read.csv(test_path("fixtures", name), ...)
}

quietly <- function(expr) suppressMessages(expr)

test_that("an fgsea table becomes an enrichment with NES scores", {
  fg <- read_fixture("fgsea_small.csv")
  x <- quietly(as_enrichment(list(Aging = fg)))
  r <- x@results
  expect_identical(x@metadata$enrichment_test, "fgsea")
  expect_identical(x@metadata$score_type, "NES")
  expect_identical(unique(r$contrast), "Aging")
  expect_identical(r$term, fg$pathway)
  expect_equal(r$score, fg$NES)
  expect_equal(r$p, fg$pval)
  expect_equal(r$padj, fg$padj)
  expect_equal(r$size, fg$size)
  expect_identical(r$direction, ifelse(fg$NES > 0, "up", "down"))
  expect_true(all(is.na(r$database)))
  expect_identical(r$leading_edge[[1]], strsplit(fg$leadingEdge[1], ";")[[1]])
})

test_that("fgsea's list-column leadingEdge is read as is", {
  fg <- read_fixture("fgsea_small.csv")
  fg_list <- fg
  fg_list$leadingEdge <- strsplit(fg$leadingEdge, ";", fixed = TRUE)
  a <- quietly(as_enrichment(list(A = fg)))
  b <- quietly(as_enrichment(list(A = fg_list)))
  expect_identical(a@results$leading_edge, b@results$leading_edge)
})

test_that("the detected test is announced once", {
  fg <- read_fixture("fgsea_small.csv")
  expect_message(
    as_enrichment(list(A = fg)),
    class = "enrichVolcano_detected_test"
  )
})

test_that("a long table keeps its contrasts, databases and extra columns", {
  path <- system.file("extdata", "examples", "yvo_fgsea.csv.gz", package = "enrichVolcano")
  ex <- utils::read.csv(path)
  x <- quietly(as_enrichment(ex))
  r <- x@results
  expect_setequal(unique(r$contrast), c("Aging", "Interaction", "Training_Old", "Training_Young"))
  expect_setequal(unique(r$database), unique(ex$database))
  expect_true(all(c("dedup_status", "merged_into", "overlap_jaccard", "ES") %in% names(r)))
  expect_false(any(c("pathway", "NES", "pval", "leadingEdge") %in% names(r)))
})

test_that("rows fgsea could not score are dropped with a note", {
  fg <- read_fixture("fgsea_small.csv")
  fg$NES[2] <- NA
  expect_message(x <- as_enrichment(list(A = fg)), "1 row")
  expect_identical(nrow(x@results), nrow(fg) - 1L)
})

test_that("a database label stamps every row of an unlabelled input", {
  fg <- read_fixture("fgsea_small.csv")
  x <- quietly(as_enrichment(list(A = fg), database = "Reactome"))
  expect_true(all(x@results$database == "Reactome"))
})

test_that("a database label clashing with a database column is refused", {
  fg <- read_fixture("fgsea_small.csv")
  fg$database <- "Hallmark"
  expect_error(
    as_enrichment(list(A = fg), database = "Reactome"),
    class = "enrichVolcano_input_error"
  )
})

test_that("contrasts must be named", {
  fg <- read_fixture("fgsea_small.csv")
  expect_error(as_enrichment(fg), class = "enrichVolcano_input_error")
  expect_error(as_enrichment(list(fg, fg)), class = "enrichVolcano_input_error")
  expect_error(as_enrichment(list(A = fg, fg)), class = "enrichVolcano_input_error")
})

test_that("list names become contrasts and override a contrast column", {
  fg <- read_fixture("fgsea_small.csv")
  fg$contrast <- "ignored"
  x <- quietly(as_enrichment(list(Young = fg, Old = fg)))
  expect_identical(unique(x@results$contrast), c("Young", "Old"))
})

test_that("an unrecognised table asks for enrichment_test", {
  tbl <- data.frame(name = "x", value = 1)
  expect_error(as_enrichment(list(A = tbl)), "enrichment_test", class = "enrichVolcano_input_error")
})

test_that("an unknown enrichment_test is refused", {
  fg <- read_fixture("fgsea_small.csv")
  expect_error(
    as_enrichment(list(A = fg), enrichment_test = "gsva"),
    class = "enrichVolcano_param_error"
  )
})

test_that("inputs are not data frames or lists of them", {
  expect_error(as_enrichment(list(A = 1:3)), class = "enrichVolcano_input_error")
  expect_error(as_enrichment("fgsea"), class = "enrichVolcano_input_error")
})

test_that("a custom table maps its columns and derives direction", {
  tbl <- data.frame(
    Name = c("SET_A", "SET_B"),
    stat = c(1.7, -2.2),
    q = c(0.01, 0.03),
    n = c(20, 31),
    genes = c("G1/G2", "G3")
  )
  x <- as_enrichment(
    list(A = tbl),
    enrichment_test = "custom",
    term = "Name", score = "stat", padj = "q", size = "n", leading_edge = "genes"
  )
  r <- x@results
  expect_identical(x@metadata$score_type, "score")
  expect_identical(r$term, c("SET_A", "SET_B"))
  expect_identical(r$direction, c("up", "down"))
  expect_true(all(is.na(r$p)))
  expect_identical(r$leading_edge, list(c("G1", "G2"), "G3"))
})

test_that("a custom table can declare its score type", {
  tbl <- data.frame(term = "SET_A", score = 1.2, padj = 0.01, size = 10)
  x <- as_enrichment(list(A = tbl), enrichment_test = "custom", score_type = "NES")
  expect_identical(x@metadata$score_type, "NES")
})

test_that("a custom table missing a mapped column names it", {
  tbl <- data.frame(term = "SET_A", padj = 0.01, size = 10)
  expect_error(
    as_enrichment(list(A = tbl), enrichment_test = "custom"),
    "score",
    class = "enrichVolcano_column_error"
  )
})

test_that("list elements produced by different tests are refused", {
  fg <- read_fixture("fgsea_small.csv")
  fry <- read_fixture("limma_fry.csv", row.names = 1)
  expect_error(
    quietly(as_enrichment(list(A = fg, B = fry))),
    class = "enrichVolcano_input_error"
  )
})

test_that("fry output becomes signed -log10(FDR) scores", {
  fry <- read_fixture("limma_fry.csv", row.names = 1)
  x <- quietly(as_enrichment(list(Aging = fry)))
  r <- x@results
  expect_identical(x@metadata$enrichment_test, "fry")
  expect_identical(x@metadata$score_type, "signed -log10(FDR)")
  expect_identical(r$term, rownames(fry))
  expect_equal(r$score, ifelse(fry$Direction == "Up", 1, -1) * -log10(fry$FDR))
  expect_equal(r$p, fry$PValue)
  expect_equal(r$padj, fry$FDR)
  expect_equal(r$size, fry$NGenes)
  expect_identical(r$direction, tolower(fry$Direction))
  expect_identical(r$leading_edge, rep(list(character(0)), nrow(fry)))
})

test_that("mroast and camera-with-correlation are recognised by their columns", {
  mroast <- read_fixture("limma_mroast.csv", row.names = 1)
  camera <- read_fixture("limma_camera_cor.csv", row.names = 1)
  expect_identical(quietly(as_enrichment(list(A = mroast)))@metadata$enrichment_test, "mroast")
  expect_identical(quietly(as_enrichment(list(A = camera)))@metadata$enrichment_test, "camera")
})

test_that("default camera and cameraPR output must be named", {
  camera <- read_fixture("limma_camera.csv", row.names = 1)
  camera_pr <- read_fixture("limma_camerapr.csv", row.names = 1)
  expect_error(as_enrichment(list(A = camera)), "cameraPR", class = "enrichVolcano_input_error")
  expect_error(as_enrichment(list(A = camera_pr)), "cameraPR", class = "enrichVolcano_input_error")
  x <- as_enrichment(list(A = camera_pr), enrichment_test = "cameraPR")
  expect_identical(x@metadata$enrichment_test, "cameraPR")
  expect_identical(x@metadata$score_type, "signed -log10(FDR)")
})

test_that("an FDR of zero gives a finite score", {
  fry <- read_fixture("limma_fry.csv", row.names = 1)
  fry$FDR[1] <- 0
  x <- quietly(as_enrichment(list(A = fry)))
  expect_true(all(is.finite(x@results$score)))
})

test_that("a stated enrichment_test must match the columns it is given", {
  fry <- read_fixture("limma_fry.csv", row.names = 1)
  expect_error(
    as_enrichment(list(A = fry), enrichment_test = "fgsea"),
    class = "enrichVolcano_column_error"
  )
})

test_that("score_type cannot relabel a known test's score", {
  fg <- read_fixture("fgsea_small.csv")
  expect_error(
    quietly(as_enrichment(list(A = fg), score_type = "t")),
    class = "enrichVolcano_param_error"
  )
})
