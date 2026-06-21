test_that("build_registry.R reproduces the shipped db_registry.csv", {
  script <- testthat::test_path("..", "..", "data-raw", "build_registry.R")
  skip_if(!file.exists(script),
          "data-raw/build_registry.R not available (installed package)")
  shipped <- system.file("extdata", "db_registry.csv",
                         package = "enrichVolcano")
  script <- normalizePath(script)
  tmp <- withr::local_tempdir()
  withr::with_dir(tmp, sys.source(script, envir = new.env()))
  built <- file.path(tmp, "inst", "extdata", "db_registry.csv")
  expect_true(file.exists(built))
  expect_identical(readLines(built), readLines(shipped))
})

test_that("list_databases returns tibble with required schema", {
  db <- list_databases()
  expect_s3_class(db, "tbl_df")
  expect_true(all(c("name", "type", "species", "license") %in% colnames(db)))
  expect_true(nrow(db) >= 17)
  expect_true(all(nchar(db$license) > 0))
})

test_that("database_info prints details for a registered DB", {
  expect_message(database_info("hallmark"), regexp = "MSigDB Hallmark")
  expect_error(database_info("nope"), class = "ev_unknown_database")
})

test_that("bundled GMT files are valid when present", {
  for (db in c("mitocarta3", "corum")) {
    for (sp in c("human", "mouse", "rat")) {
      fn <- paste0(db, "_", sp, ".gmt")
      path <- system.file("extdata", "gmt", fn, package = "enrichVolcano")
      if (!nzchar(path)) {
        skip(paste("GMT not bundled yet:", fn))
      }
      paths <- fgsea::gmtPathways(path)
      expect_true(length(paths) > 0)
    }
  }
})

test_that("hpa_subcellular and uniprot_subcellular are registered", {
  reg <- list_databases()
  for (db in c("hpa_subcellular", "uniprot_subcellular")) {
    expect_true(db %in% reg$name)
    expect_identical(reg$type[reg$name == db], "organelle")
    expect_identical(reg$source[reg$name == db], "bundled_gmt")
  }
  expect_identical(reg$species[reg$name == "hpa_subcellular"], "human")
  expect_identical(reg$species[reg$name == "uniprot_subcellular"],
                   "human;mouse;rat")
})

test_that("HPA subcellular resolves for human only", {
  f <- system.file("extdata", "gmt", "hpa_subcellular_human.gmt",
                   package = "enrichVolcano")
  skip_if(!nzchar(f), "hpa_subcellular GMT not bundled")
  res <- ev_resolve_databases("hpa_subcellular", species = "human")
  expect_gt(length(res[["hpa_subcellular"]]), 30)
  expect_true(all(grepl("^HPA_", names(res[["hpa_subcellular"]]))))
  # human-only: mouse warns and resolves to nothing
  expect_warning(
    empty <- ev_resolve_databases("hpa_subcellular", species = "mouse"),
    class = "ev_species_unsupported"
  )
  expect_length(empty, 0)
})

test_that("UniProt subcellular resolves for human, mouse, and rat", {
  for (sp in c("human", "mouse", "rat")) {
    f <- system.file("extdata", "gmt",
                     paste0("uniprot_subcellular_", sp, ".gmt"),
                     package = "enrichVolcano")
    skip_if(!nzchar(f), paste("uniprot_subcellular GMT not bundled:", sp))
    res <- ev_resolve_databases("uniprot_subcellular", species = sp)
    p <- res[["uniprot_subcellular"]]
    expect_gt(length(p), 50)
    expect_true(all(grepl("^UNIPROT_", names(p))))
  }
})

test_that("ev_parse_uniprot_cc reduces CC text to top-level compartments", {
  script <- testthat::test_path("..", "..", "data-raw", "build_gmts.R")
  skip_if(!file.exists(script),
          "data-raw/build_gmts.R not available (installed package)")
  script <- normalizePath(script)
  env <- new.env()
  tmp <- withr::local_tempdir()
  # Sourced in an empty tempdir so the script's build_*() calls find no raw
  # inputs and no-op; we only want the parser definition.
  suppressMessages(
    withr::with_dir(tmp, sys.source(script, envir = env))
  )
  f <- env$ev_parse_uniprot_cc
  expect_identical(
    f(paste0("SUBCELLULAR LOCATION: Membrane {ECO:0000305}; Multi-pass ",
             "membrane protein {ECO:0000255}. Cytoplasm {ECO:0000305}. ",
             "Mitochondrion {ECO:0000269|PubMed:1}.")),
    c("Membrane", "Cytoplasm", "Mitochondrion")
  )
  # comma sub-locations collapse to the leading compartment; Note= dropped
  expect_identical(
    f("SUBCELLULAR LOCATION: Cytoplasm, cytoskeleton, spindle. Note=Free."),
    "Cytoplasm"
  )
  # isoform markers stripped, both locations kept
  expect_identical(
    f("SUBCELLULAR LOCATION: [Isoform 1]: Nucleus. [Isoform 2]: Cytoplasm."),
    c("Nucleus", "Cytoplasm")
  )
})

test_that("go_slim is registered with the canonical OBO loader", {
  reg <- list_databases()
  expect_true("go_slim" %in% reg$name)
  expect_identical(reg$source[reg$name == "go_slim"], "go_slim_obo")
})

test_that("go_slim OBO file is shipped and parses cleanly", {
  obo <- system.file("extdata", "goslim_generic.obo",
                     package = "enrichVolcano")
  skip_if(!nzchar(obo), "goslim_generic.obo not installed")
  bp <- enrichVolcano:::ev_parse_goslim_obo(obo)
  # The published GO generic slim BP subset has dozens of terms; we don't
  # pin a brittle count, only that the parser finds a reasonable number.
  expect_gt(nrow(bp), 30L)
  expect_true(all(grepl("^GO:\\d{7}$", bp$go_id)))
  expect_true(all(nzchar(bp$term)))
})

test_that("load_go_slim builds size-filtered GOSLIM_ sets per species", {
  for (sp in c("human", "mouse", "rat")) {
    orgdb <- switch(sp,
      human = "org.Hs.eg.db",
      mouse = "org.Mm.eg.db",
      rat   = "org.Rn.eg.db")
    skip_if_not_installed("GO.db")
    skip_if_not_installed("AnnotationDbi")
    skip_if_not_installed(orgdb)
    sets <- load_go_slim(species = sp)
    expect_gt(length(sets), 20L)
    expect_true(all(lengths(sets) >= 10L & lengths(sets) <= 500L))
    expect_true(all(grepl("^GOSLIM_", names(sets))))
  }
})
