# Regenerates tests/testthat/fixtures/toy_study/ and toy_study.xlsx: 20 human
# proteins (10 OXPHOS, 10 glycolysis), 4 subjects measured Pre and Post, OXPHOS
# raised after training. da_results come from the limma fit the package refits.

library(limma)

fixture <- function(...) here::here("tests", "testthat", "fixtures", ...)

proteins <- c(
  "P31040", "P21912", "O00483", "P13073", "P25705", "O75390", "P50213", "P40926", "P99999", "P31930",
  "P52789", "P14618", "P00338", "P06733", "P04406", "P08237", "P04075", "P60174", "P00558", "P18669"
)
samples <- data.frame(
  sample = sprintf("S%02d", 1:8),
  group = rep(c("Pre", "Post"), each = 4),
  subject = rep(sprintf("P%d", 1:4), 2)
)

set.seed(11)
subject_effect <- rnorm(4, sd = 0.5)
matrix <- matrix(rnorm(20 * 8, mean = 20), 20, 8, dimnames = list(proteins, samples$sample))
matrix <- sweep(matrix, 2, rep(subject_effect, 2), "+")
matrix[1:10, samples$group == "Post"] <- matrix[1:10, samples$group == "Post"] + 1
weights <- matrix(round(runif(20 * 8, 0.5, 1.5), 3), 20, 8, dimnames = dimnames(matrix))

design <- model.matrix(~ 0 + group, samples)
colnames(design) <- sub("^group", "", colnames(design))
contrasts <- data.frame(name = "Training", expression = "Post - Pre")
cm <- makeContrasts(contrasts = contrasts$expression, levels = design)
colnames(cm) <- contrasts$name
corr <- duplicateCorrelation(matrix, design, block = samples$subject)$consensus.correlation
fit <- eBayes(contrasts.fit(lmFit(matrix, design, block = samples$subject, correlation = corr), cm))
da <- topTable(fit, coef = "Training", number = Inf, sort.by = "none")
da_results <- data.frame(
  protein = rownames(da), contrast = "Training", logFC = da$logFC, t = da$t,
  p = da$P.Value, padj = da$adj.P.Val, abundance = da$AveExpr
)

sheets <- list(
  matrix = data.frame(protein = proteins, round(matrix, 4), check.names = FALSE),
  samples = samples,
  contrasts = contrasts,
  da_results = da_results,
  weights = data.frame(protein = proteins, weights, check.names = FALSE)
)
for (nm in names(sheets)) utils::write.csv(sheets[[nm]], fixture("toy_study", paste0(nm, ".csv")), row.names = FALSE)
openxlsx::write.xlsx(sheets, fixture("toy_study.xlsx"))
