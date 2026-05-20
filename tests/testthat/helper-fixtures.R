make_tidy_minimal <- function() {
  set.seed(42)
  tibble::tibble(
    gene = paste0("GENE", 1:50),
    contrast = rep(c("A_vs_B", "C_vs_D"), each = 25),
    logFC = rnorm(50),
    P.Value = runif(50),
    adj.P.Val = runif(50)
  )
}

make_mini_pathways <- function() {
  list(
    HALLMARK_APOPTOSIS = paste0("G", 1:20),
    HALLMARK_GLYCOLYSIS = paste0("G", 15:35),
    HALLMARK_HYPOXIA = paste0("G", 40:60)
  )
}

make_ranked_input <- function() {
  set.seed(123)
  tibble::tibble(
    gene = paste0("G", 1:100),
    contrast = "ctr",
    logFC = c(rnorm(20, mean = 3), rnorm(60), rnorm(20, mean = -3)),
    P.Value = c(runif(20, max = 0.01), runif(60), runif(20, max = 0.01)),
    adj.P.Val = runif(100),
    pi_eq2 = runif(100, max = 0.5)
  )
}
