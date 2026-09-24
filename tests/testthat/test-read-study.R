toy_dir <- function() test_path("fixtures", "toy_study")

# A temporary copy of the toy study with one sheet edited.
edited_study <- function(sheet = NULL, edit = identity, drop = character(0), env = parent.frame()) {
  dir <- withr::local_tempdir(.local_envir = env)
  file.copy(list.files(toy_dir(), full.names = TRUE), dir)
  unlink(file.path(dir, paste0(drop, ".csv")))
  if (!is.null(sheet)) {
    path <- file.path(dir, paste0(sheet, ".csv"))
    utils::write.csv(edit(utils::read.csv(path, check.names = FALSE)), path, row.names = FALSE)
  }
  dir
}

quietly <- function(expr) suppressMessages(expr)

test_that("a CSV study reads into validated parts", {
  s <- quietly(read_study(toy_dir()))
  expect_s3_class(s, "enrichVolcano_study")
  expect_s3_class(s$da, "enrichVolcano_da")
  expect_true(is.matrix(s$matrix) && is.numeric(s$matrix))
  expect_identical(dim(s$matrix), c(20L, 8L))
  expect_identical(colnames(s$matrix), s$samples$sample)
  expect_identical(dim(s$weights), dim(s$matrix))
  expect_identical(s$contrasts$name, "Training")
})

test_that("the Excel workbook reads the same as the CSV folder", {
  skip_if_not_installed("readxl")
  a <- quietly(read_study(toy_dir()))
  b <- quietly(read_study(test_path("fixtures", "toy_study.xlsx")))
  expect_equal(b$matrix, a$matrix)
  expect_equal(b$da, a$da)
  expect_equal(as.data.frame(b$samples), a$samples)
})

test_that("gzipped CSVs are read", {
  dir <- edited_study()
  for (f in list.files(dir, full.names = TRUE)) {
    writeLines(readLines(f), gz <- gzfile(paste0(f, ".gz"), "w"))
    close(gz)
    unlink(f)
  }
  expect_identical(dim(quietly(read_study(dir))$matrix), c(20L, 8L))
})

test_that("da_results alone is a study for fgsea and figures", {
  s <- quietly(read_study(edited_study(drop = c("matrix", "samples", "contrasts", "weights"))))
  expect_null(s$matrix)
  expect_identical(nrow(s$da), 20L)
})

test_that("the protein overlap is reported", {
  expect_message(read_study(toy_dir()), "20 of 20", class = "enrichVolcano_study_overlap")
})

test_that("a partial design is refused", {
  expect_error(
    read_study(edited_study(drop = "contrasts")),
    "contrasts",
    class = "enrichVolcano_input_error"
  )
  expect_error(
    read_study(edited_study(drop = c("matrix", "samples", "contrasts"))),
    "weights",
    class = "enrichVolcano_input_error"
  )
})

test_that("matrix columns and samples must name the same samples", {
  dir <- edited_study("samples", function(d) {
    d$sample[1] <- "S99"
    d
  })
  expect_error(read_study(dir), "S99", class = "enrichVolcano_input_error")
})

test_that("group names must be usable in contrast arithmetic", {
  dir <- edited_study("samples", function(d) {
    d$group[d$group == "Post"] <- "Post-train"
    d
  })
  expect_error(quietly(read_study(dir)), "Post-train", class = "enrichVolcano_input_error")
})

test_that("contrasts may only use existing groups", {
  dir <- edited_study("contrasts", function(d) {
    d$expression <- "Post - Baseline"
    d
  })
  expect_error(quietly(read_study(dir)), "Baseline", class = "enrichVolcano_input_error")
})

test_that("DA contrasts must match the contrasts sheet", {
  dir <- edited_study("da_results", function(d) {
    d$contrast <- "Train"
    d
  })
  expect_error(quietly(read_study(dir)), "Train", class = "enrichVolcano_input_error")
})

test_that("weights must match the matrix", {
  dir <- edited_study("weights", function(d) d[-1, ])
  expect_error(quietly(read_study(dir)), "weights", class = "enrichVolcano_input_error")
})

test_that("a matrix with duplicate or non-numeric rows is refused", {
  dup <- edited_study("matrix", function(d) rbind(d, d[1, ]))
  expect_error(quietly(read_study(dup)), "duplicate", class = "enrichVolcano_data_error")
  text <- edited_study("matrix", function(d) {
    d$S01 <- "high"
    d
  })
  expect_error(quietly(read_study(text)), "numeric", class = "enrichVolcano_data_error")
})

test_that("required columns and files are named when missing", {
  expect_error(
    quietly(read_study(edited_study("samples", function(d) d[c("sample", "subject")]))),
    "group",
    class = "enrichVolcano_column_error"
  )
  expect_error(
    quietly(read_study(edited_study("contrasts", function(d) d["name"]))),
    "expression",
    class = "enrichVolcano_column_error"
  )
  expect_error(read_study(edited_study(drop = "da_results")), "da_results", class = "enrichVolcano_input_error")
  expect_error(read_study(file.path(tempdir(), "no_such_study")), class = "enrichVolcano_input_error")
})

test_that("example_study() lists the shipped studies", {
  idx <- example_study()
  expect_setequal(idx$name, c("bfr_limpa", "mouse_pas", "yvo", "bfr_proteoda", "cvh", "hrvlr", "mito"))
  expect_true(all(c("species", "description") %in% names(idx)))
})

test_that("a shipped study loads with its species", {
  skip_if_not_installed("org.Mm.eg.db")
  s <- quietly(example_study("mouse_pas"))
  expect_identical(dim(s$matrix), c(1842L, 20L))
  expect_identical(dim(s$weights), dim(s$matrix))
  expect_identical(s$info$species, "Mus musculus")
  expect_true("Sdha" %in% s$da$gene)
})

test_that("a DA-only shipped study has no matrix", {
  skip_if_not_installed("org.Hs.eg.db")
  s <- quietly(example_study("cvh"))
  expect_null(s$matrix)
  expect_identical(length(unique(s$da$contrast)), 7L)
})

test_that("the YvO study ships its limpa matrix, weights and design", {
  skip_if_not_installed("org.Hs.eg.db")
  s <- quietly(example_study("yvo"))
  expect_identical(dim(s$matrix), c(2106L, 62L))
  expect_identical(dim(s$weights), dim(s$matrix))
  expect_identical(length(unique(s$samples$subject)), 32L)
  expect_setequal(s$contrasts$name, c("Training_Young", "Training_Old", "Aging", "Interaction"))
})

test_that("an unknown example study is refused", {
  expect_error(example_study("nope"), "nope", class = "enrichVolcano_param_error")
})

test_that("numeric sample IDs select matrix columns by name, not position", {
  base <- quietly(read_study(toy_dir()))
  dir <- edited_study("samples", function(s) transform(s, sample = rev(seq_along(sample))), drop = "weights")
  m <- utils::read.csv(file.path(dir, "matrix.csv"), check.names = FALSE)
  names(m)[-1] <- rev(seq_len(ncol(m) - 1))
  utils::write.csv(m, file.path(dir, "matrix.csv"), row.names = FALSE)
  s <- quietly(read_study(dir))
  expect_identical(colnames(s$matrix), as.character(rev(seq_len(ncol(m) - 1))))
  expect_equal(unname(s$matrix), unname(base$matrix))
})

test_that("a sample listed twice is refused", {
  dir <- edited_study("samples", function(s) rbind(s, s[1, ]), drop = "weights")
  expect_error(quietly(read_study(dir)), "twice", class = "enrichVolcano_input_error")
})
