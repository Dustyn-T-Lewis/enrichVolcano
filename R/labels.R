# Shared rule-driven label selection for the standalone volcano and the
# composite. Default mode "none" labels nothing.

#' Label text for a row: symbol, else accession, else gene.
#' @keywords internal
#' @noRd
ev_label_text <- function(df) {
  n <- nrow(df)
  symbol <- if ("symbol" %in% names(df)) df$symbol else rep(NA_character_, n)
  uniprot <- if ("uniprot" %in% names(df)) df$uniprot else rep(NA_character_, n)
  gene <- if ("gene" %in% names(df)) df$gene else rep(NA_character_, n)
  out <- ifelse(!is.na(symbol) & nzchar(symbol), symbol,
    ifelse(!is.na(uniprot) & nzchar(uniprot), uniprot, gene)
  )
  as.character(out)
}

#' Select which rows to label on a volcano
#'
#' @param df Single-contrast tibble with `logFC` and the significance column.
#' @param mode One of `"none"`, `"top_per_direction"`, `"top_total"`,
#'   `"all_significant"`, `"explicit"`.
#' @param n Cap for the top-* modes (per direction for `top_per_direction`).
#' @param rank_by `"significance"` (smaller `p_col` first) or `"logfc"`
#'   (larger `|logFC|` first).
#' @param genes For `"explicit"`: matched against symbol, accession, or gene.
#' @param p_col Significance column (e.g. `"pi_eq2"`, `"P.Value"`,
#'   `"adj.P.Val"`).
#' @param p_threshold,logfc_threshold Cutoffs defining the candidate pool.
#' @return Sub-tibble of rows to label, with a `label_text` column.
#' @keywords internal
#' @noRd
ev_select_labels <- function(df, mode, n, rank_by, genes,
                             p_col, p_threshold, logfc_threshold) {
  df$label_text <- ev_label_text(df)
  if (mode == "explicit") {
    in_genes <- function(col) {
      if (col %in% names(df)) df[[col]] %in% genes else rep(FALSE, nrow(df))
    }
    keep <- df$label_text %in% genes | in_genes("symbol") |
      in_genes("uniprot") | in_genes("gene")
    return(df[keep, , drop = FALSE])
  }
  if (mode == "none") {
    return(df[0, , drop = FALSE])
  }

  sig <- df[[p_col]] < p_threshold & abs(df$logFC) >= logfc_threshold
  pool <- df[sig & !is.na(sig), , drop = FALSE]
  if (nrow(pool) == 0) {
    return(pool)
  }

  ord <- if (rank_by == "logfc") {
    order(-abs(pool$logFC))
  } else {
    order(pool[[p_col]])
  }
  pool <- pool[ord, , drop = FALSE]

  switch(mode,
    all_significant = pool,
    top_total = utils::head(pool, n),
    top_per_direction = {
      # logFC == 0 belongs to neither direction and is intentionally excluded.
      up <- pool[pool$logFC > 0, , drop = FALSE]
      dn <- pool[pool$logFC < 0, , drop = FALSE]
      rbind(utils::head(up, n), utils::head(dn, n))
    },
    ev_abort("Unknown label mode {.val {mode}}.", class = "ev_bad_label_mode")
  )
}

# Pathway-name cleaning for ring display: strips DB prefixes, expands
# acronyms, shortens long phrases, and wraps for the arc labels.

#' Clean a pathway name for ring display
#'
#' Strips the database prefix, expands common acronyms, title-cases, wraps,
#' and applies MitoCarta-hierarchy shortening for `MITOCARTA_` pathways.
#'
#' @param name Character vector of raw pathway names.
#' @return Character vector of cleaned, wrapped labels.
#' @export
#' @examples
#' ev_clean_label(c("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "REACTOME_TCA_CYCLE"))
ev_clean_label <- function(name) {
  if (length(name) > 1) {
    return(vapply(name, ev_clean_label, character(1), USE.NAMES = FALSE))
  }
  if (is.na(name) || !nzchar(name)) {
    return(name)
  }

  if (startsWith(name, "MITOCARTA_")) {
    return(ev_clean_label_mitocarta(name))
  }

  out <- name |>
    sub("^HALLMARK_", "", x = _) |>
    sub("^REACTOME_", "", x = _) |>
    sub("^GOSLIM_", "", x = _) |>
    sub("^KEGG_(MEDICUS|LEGACY)_", "", x = _) |>
    sub("^KEGG_", "", x = _) |>
    sub("^GO(BP|CC|MF)_", "", x = _) |>
    sub("^WP_", "", x = _) |>
    sub("^CORUM_", "", x = _) |>
    gsub("_", " ", x = _) |>
    tolower() |>
    stringr::str_to_title()

  out <- ev_expand_acronyms(out)
  out <- ev_shorten_phrases(out)
  out <- ev_collapse_repeats(out)
  out <- stringr::str_wrap(out, width = 15)
  ev_post_wrap_overrides(trimws(out))
}

#' Collapse adjacent duplicate words left by overlapping MSigDB hierarchies
#' ("OXPHOS OXPHOS Subunits" -> "OXPHOS Subunits").
#' @keywords internal
#' @noRd
ev_collapse_repeats <- function(x) {
  repeat {
    y <- gsub("\\b([A-Za-z][A-Za-z0-9.]*)\\b\\s+\\1\\b", "\\1", x,
      perl = TRUE, ignore.case = TRUE
    )
    if (identical(y, x)) break
    x <- y
  }
  trimws(gsub("\\s+", " ", x))
}

#' Manual line-break overrides applied after wrapping, for a handful of names
#' that otherwise wrap awkwardly on the ring.
#' @keywords internal
#' @noRd
ev_post_wrap_overrides <- function(x) {
  x <- sub("Protein\nLocalization To\nPlasma Membrane",
    "Protein Localiz.\nto Plasma\nMem.", x,
    fixed = TRUE
  )
  x <- sub("(?s).*Maintenance.*Cell.*Polarity.*", "Maintenance\nof Polarity",
    x,
    perl = TRUE
  )
  x <- sub("^Heme Metabolism$", "Heme\nMetabolism", x)
  x <- sub("^tRNA Metabolism$", "tRNA\nMetabolism", x)
  x <- sub("^Mitotic Spindle$", "Mitotic\nSpindle", x)
  x <- sub("^MYC Targets V1$", "MYC Targets\nV1", x)
  x <- sub("^UV Response Dn$", "UV Response\nDn", x)
  x
}

#' @keywords internal
#' @noRd
ev_expand_acronyms <- function(x) {
  reps <- c(
    "\\bTca\\b" = "TCA", "\\bMapk(\\d?)\\b" = "MAPK\\1",
    "\\bPtk(\\d?)\\b" = "PTK\\1", "\\bRhoa\\b" = "RhoA",
    "\\bRhob\\b" = "RhoB", "\\bRhoc\\b" = "RhoC",
    "\\bGtpase\\b" = "GTPase", "\\bGtp\\b" = "GTP", "\\bGdp\\b" = "GDP",
    "\\bMrna\\b" = "mRNA", "\\bRna\\b" = "RNA", "\\bDna\\b" = "DNA",
    "\\bAtp\\b" = "ATP", "\\bNadh\\b" = "NADH", "\\bFadh\\b" = "FADH",
    "\\bRos\\b" = "ROS", "\\bIfn\\b" = "IFN", "\\bIl-?(\\d+)\\b" = "IL\\1",
    "\\bPi3k\\b" = "PI3K", "\\bAkt\\b" = "AKT", "\\bMtor\\b" = "mTOR",
    "\\bMtorc1\\b" = "MTORC1", "\\bE2f\\b" = "E2F", "\\bMyc\\b" = "MYC",
    "\\bKras\\b" = "KRAS", "\\bTnfa\\b" = "TNFA", "\\bNfkb\\b" = "NF-kB",
    "\\bG2m\\b" = "G2M", "\\bUv\\b" = "UV", "\\bWnt\\b" = "WNT",
    "\\bStat(\\d)\\b" = "STAT\\1", "\\bJak\\b" = "JAK", "\\bTgf\\b" = "TGF",
    "\\bPirnas?\\b" = "piRNAs", "\\bDgc\\b" = "DGC", "\\bMpc\\b" = "MPC",
    "\\bP53\\b" = "p53", "\\bSlc25a\\b" = "SLC25A",
    "Trna " = "tRNA "
  )
  for (pat in names(reps)) x <- gsub(pat, reps[[pat]], x, perl = TRUE)
  x
}

#' @keywords internal
#' @noRd
ev_shorten_phrases <- function(x) {
  reps <- c(
    # whole-name overrides for a few terms that stay unwieldy after shortening
    "Immunoregulatory Interactions Between A Lymphoid And A Non Lymphoid Cell" =
      "Lymphoid Cell Interactions",
    "Transcriptional And Post Translational Regulation Of Mitf M Expression And Activity" =
      "MITF-M Regulation",
    "Arrhythmogenic Right Ventricular Cardiomyopathy Arvc" = "ARVC",
    "Cargo Recognition For Clathrin Mediated Endocytosis" = "Clathrin Endocytosis",
    "Assembly Of Collagen Fibrils And Other Multimeric Structures" =
      "Collagen Fibril Assembly",
    "Collagen Chain Trimerization" = "Collagen Trimerization",
    "Processing Of Capped Intron Containing Pre mRNA" = "Pre-mRNA Processing",
    "Respiratory Chain Complex I (Holoenzyme), Mitochondrial" = "Respiratory Complex I",
    "Respiratory Chain Complex I, Mitochondrial" = "Respiratory Complex I",
    "Mitochondrial Ribosome, Large Subunit" = "Mitoribosome (Large)",
    "Mitochondrial Ribosome, Small Subunit" = "Mitoribosome (Small)",
    "Oxidative Phosphorylation" = "OXPHOS",
    "Unfolded Protein Response" = "UPR", "Fatty Acid" = "FA",
    "Amino Acid" = "AA", "Generation Of" = "Gen. of", " And " = " & ",
    "Mitochondrion" = "Mito.", "Mitochondrial" = "Mito.",
    "Ubiquinone" = "UQ", "Organization" = "Org.",
    "Cytoskeleton" = "Cytoskel.", "Microtubule" = "MT",
    "Respiratory" = "Resp.", "Electron Transport" = "ETC",
    "Synthesis Coupled" = "Synth.-Coupled",
    "Ubiquitin Dependent" = "Ub-Dep.",
    "Proteasome Mediated" = "Proteasome-Med.", "Proteasomal" = "Proteas.",
    "Phosphorylation" = "Phosph.",
    "Modification" = "Mod.", "Intracellular" = "Intracell.",
    "Regulation Of" = "Reg.", "Signaling Pathway" = "Signaling",
    "Biosynthetic Process" = "Biosynthesis",
    "Catabolic Process" = "Catabolism", "Metabolic Process" = "Metabolism",
    "Based Process" = "Process", "Response To" = "Resp. to",
    "Establishment Or Maintenance Of" = "Maintenance of",
    "Extracellular Matrix" = "ECM",
    "Epithelial Mesenchymal Transition" = "EMT",
    "Precursor Metabolites & Energy" = "Precursor Metabolites",
    "Citric Acid Cycle TCA Cycle" = "TCA Cycle",
    "Aerobic Respiration & Resp. ETC" = "Aerobic Resp. + ETC"
  )
  for (pat in names(reps)) x <- gsub(pat, reps[[pat]], x, fixed = TRUE)
  # variable-suffix MSigDB names that title-case to a long phrase
  x <- sub(
    "External Encapsulating Structure Or.*",
    "Extracellular Matrix Org.", x
  )
  x <- sub(
    "Enzyme Linked Receptor Protein Signaling.*",
    "Receptor Protein Signaling", x
  )
  x
}

#' @keywords internal
#' @noRd
ev_clean_label_mitocarta <- function(name) {
  n <- sub("^MITOCARTA_", "", name)
  parts <- strsplit(n, "__|>", perl = TRUE)[[1]]
  if (length(parts) == 1) parts <- strsplit(n, "_")[[1]]
  leaf <- utils::tail(parts, 1)
  leaf <- gsub("_", " ", leaf)
  leaf <- tools::toTitleCase(tolower(leaf))
  reps <- c(
    "\\bOxphos\\b" = "OXPHOS", "\\bTca\\b" = "TCA", "\\bAa\\b" = "AA",
    "\\bFa\\b" = "FA", "\\bRos\\b" = "ROS", "\\bImm\\b" = "IMM",
    "\\bOmm\\b" = "OMM", "\\bIms\\b" = "IMS",
    "\\bCi Subunits\\b" = "Complex I", "\\bCii Subunits\\b" = "Complex II",
    "\\bCiii Subunits\\b" = "Complex III", "\\bCiv Subunits\\b" = "Complex IV",
    "\\bCv Subunits\\b" = "Complex V", "\\bMitochondrial\\b" = "Mito.",
    "\\bMitochondrion\\b" = "Mito.", "\\bAmino Acid\\b" = "AA",
    "\\bFatty Acid\\b" = "FA"
  )
  for (pat in names(reps)) leaf <- gsub(pat, reps[[pat]], leaf, perl = TRUE)
  leaf <- trimws(gsub("\\s+", " ", leaf))
  stringr::str_wrap(leaf, width = 15)
}
