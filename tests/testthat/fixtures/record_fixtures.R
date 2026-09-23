# Regenerates the recorded fgsea and limma fixtures. Not run by the tests;
# rerun by hand if the upstream output format changes.

library(fgsea)
library(limma)

fixture <- function(name) here::here("tests", "testthat", "fixtures", name)

data(examplePathways, package = "fgsea")
data(exampleRanks, package = "fgsea")
set.seed(1)
fg <- fgsea(examplePathways[1:30], exampleRanks, minSize = 15, maxSize = 500)
fg <- fg[order(fg$padj)][1:8]
fg$leadingEdge <- vapply(fg$leadingEdge, paste, character(1), collapse = ";")
write.csv(fg, fixture("fgsea_small.csv"), row.names = FALSE)

set.seed(1)
y <- matrix(rnorm(100 * 6), 100, 6, dimnames = list(paste0("G", 1:100), NULL))
y[1:10, 4:6] <- y[1:10, 4:6] + 2
y[11:20, 4:6] <- y[11:20, 4:6] - 1
design <- cbind(Intercept = 1, Grp = c(0, 0, 0, 1, 1, 1))
idx <- list(SET_UP = 1:10, SET_DOWN = 11:20, SET_NULL = 41:60)
t_stat <- eBayes(lmFit(y, design))$t[, "Grp"]

write.csv(fry(y, idx, design), fixture("limma_fry.csv"))
write.csv(mroast(y, idx, design, nrot = 999), fixture("limma_mroast.csv"))
write.csv(camera(y, idx, design), fixture("limma_camera.csv"))
write.csv(camera(y, idx, design, inter.gene.cor = NA), fixture("limma_camera_cor.csv"))
write.csv(cameraPR(t_stat, idx), fixture("limma_camerapr.csv"))

t2g <- utils::stack(examplePathways[1:30])[, c("ind", "values")]
gsea <- clusterProfiler::GSEA(
  sort(exampleRanks, decreasing = TRUE),
  TERM2GENE = t2g, pvalueCutoff = 1, minGSSize = 15, maxGSSize = 500,
  seed = TRUE, verbose = FALSE, eps = 0
)
write.csv(head(gsea@result, 6), fixture("gsea_result_small.csv"), row.names = FALSE)
