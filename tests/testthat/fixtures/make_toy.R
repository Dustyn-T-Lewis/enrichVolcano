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

# The toy enrichment table as an enrichment object, one copy per contrast.
make_toy_ring_enrichment <- function(contrasts = "toy") {
  e <- make_toy_enrich()
  suppressMessages(as_enrichment(
    stats::setNames(rep(list(e), length(contrasts)), contrasts),
    enrichment_test = "custom", term = "pathway", score = "NES",
    leading_edge = "leading_edge", score_type = "NES"
  ))
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

# Gene sets with hand-computable overlaps (A vs B: J = 8/12, overlap = 8/10;
# A vs C: J = 4/10, overlap = 4/4; D is disjoint from everything).
make_toy_gene_sets <- function() {
  list(
    SET_A = paste0("G", 1:10),
    SET_B = paste0("G", c(1:8, 11:12)),
    SET_C = paste0("G", 1:4),
    SET_D = paste0("G", 20:29),
    SET_E = paste0("G", 1:10),
    SET_F = paste0("G", 1:10)
  )
}

# One contrast "K" in padj order A < B < C < D within Hallmark; E is a copy of
# A in another database; F is a copy of A that is not significant.
make_toy_dedup_enrichment <- function() {
  results <- data.frame(
    contrast = "K",
    database = c("Hallmark", "Hallmark", "Hallmark", "Hallmark", "Reactome", "Hallmark"),
    term = paste0("SET_", LETTERS[1:6]),
    score = c(2.5, 2.1, 1.9, -1.8, 2.4, 0.4),
    p = NA_real_,
    padj = c(0.001, 0.002, 0.003, 0.004, 0.0005, 0.5),
    size = 10,
    direction = c("up", "up", "up", "down", "up", "up"),
    stringsAsFactors = FALSE
  )
  results$leading_edge <- rep(list(character(0)), nrow(results))
  enrichment(results = results, metadata = make_toy_metadata())
}

# Two contrasts, 12 terms, significance at 0.05 designed by hand:
# Both = T01 T03 T04 T10 T12; A only = T02 T09; B only = T05 T11; NS = T06-T08.
# Significant quadrants: TR 4 (T01 T05 T09 T12), TL 1 (T04), BL 2 (T03 T10),
# BR 2 (T02 T11); 6 of 9 concordant. T12 is below the default label size.
make_toy_scatter_enrichment <- function() {
  term <- sprintf("T%02d", 1:12)
  score_a <- c(2.0, 1.5, -1.8, -2.2, 1.2, -1.0, 0.5, -0.4, 1.9, -1.6, 0.8, 2.4)
  score_b <- c(1.8, -1.2, -2.0, 1.4, 1.1, -0.9, 0.3, -0.6, 2.1, -1.9, -0.7, 2.2)
  padj_a <- c(.001, .01, .002, .02, .2, .3, .6, .7, .004, .03, .5, .0005)
  padj_b <- c(.002, .4, .001, .03, .01, .5, .8, .6, .3, .001, .04, .001)
  one <- function(contrast, score, padj) {
    data.frame(
      contrast = contrast, term = term, score = score, padj = padj,
      database = rep(c("Hallmark", "GO Slim"), c(8, 4)),
      size = c(rep(20, 11), 10), stringsAsFactors = FALSE
    )
  }
  tbl <- rbind(one("A", score_a, padj_a), one("B", score_b, padj_b))
  suppressMessages(as_enrichment(tbl, enrichment_test = "custom", score_type = "NES"))
}
