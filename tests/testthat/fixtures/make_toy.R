# Toy-frame factory used by every Phase 3 test. Produces a 20-gene volcano
# + 5-pathway enrichment table using the package's conventional column
# names. Calling code controls the seed for reproducibility.

make_toy_volc <- function(seed = 1L, n = 20L) {
  set.seed(seed)
  half <- n %/% 2
  data.frame(
    gene = paste0("G", seq_len(n)),
    contrast = "toy",
    logFC = c(rnorm(half, 2, 0.4), rnorm(n - half, -2, 0.4)),
    P.Value = c(runif(half, 0, 0.01), runif(n - half, 0, 0.01)),
    padj = c(runif(half, 0, 0.04), runif(n - half, 0, 0.04)),
    pi_eq2 = c(runif(half, 0, 0.5), runif(n - half, 0, 0.5)),
    stringsAsFactors = FALSE
  )
}

make_toy_enrich <- function(seed = 1L) {
  set.seed(seed)
  data.frame(
    pathway = paste0("HALLMARK_TOY_", LETTERS[1:5]),
    NES = c(2.4, 1.8, -1.5, -2.1, 0.4),
    padj = c(0.001, 0.02, 0.01, 0.005, 0.5),
    size = c(40L, 30L, 25L, 35L, 18L),
    leading_edge = c(
      "G1;G2;G3;G4",
      "G5;G6",
      "G11;G12",
      "G13;G14;G15;G16",
      "G7;G17"
    ),
    stringsAsFactors = FALSE
  )
}

# Same enrichment table with fgsea-native list-column for leadingEdge.
make_toy_enrich_listcol <- function(seed = 1L) {
  d <- make_toy_enrich(seed = seed)
  d$leadingEdge <- strsplit(d$leading_edge, ";", fixed = TRUE)
  d$leading_edge <- NULL
  d
}

# clusterProfiler-style: core_enrichment as "/"-separated string.
make_toy_enrich_cp <- function(seed = 1L) {
  d <- make_toy_enrich(seed = seed)
  d$core_enrichment <- gsub(";", "/", d$leading_edge, fixed = TRUE)
  d$leading_edge <- NULL
  d
}

# enrichR-style: Genes column with ";" separator.
make_toy_enrich_enrichr <- function(seed = 1L) {
  d <- make_toy_enrich(seed = seed)
  d$Genes <- d$leading_edge
  d$leading_edge <- NULL
  d
}

# A valid `results` table for the enrichment class: two contrasts, one
# unlabelled-database row, leading edges pointing at make_toy_volc() genes.
make_toy_results <- function() {
  data.frame(
    contrast = rep(c("A", "B"), each = 3),
    database = c("Hallmark", "Hallmark", NA, "Hallmark", "Hallmark", NA),
    term = rep(c("HALLMARK_TOY_A", "HALLMARK_TOY_B", "SET_C"), 2),
    score = c(2.4, -1.5, 0.8, 1.9, -2.2, -0.3),
    p = c(0.001, 0.01, NA, 0.002, 0.004, 0.6),
    padj = c(0.004, 0.02, 0.4, 0.006, 0.008, 0.7),
    size = c(40L, 25L, 18L, 40L, 25L, 18L),
    direction = c("up", "down", "up", "up", "down", "down"),
    leading_edge = I(list(
      c("G1", "G2"), c("G11", "G12"), character(0),
      c("G1", "G3"), c("G13"), character(0)
    )),
    stringsAsFactors = FALSE
  )
}

make_toy_metadata <- function(enrichment_test = "fgsea", score_type = "NES") {
  list(enrichment_test = enrichment_test, score_type = score_type, dedup = NULL)
}
