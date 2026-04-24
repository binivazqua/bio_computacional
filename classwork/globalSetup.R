# Shared setup for sequence analysis reports.

if (requireNamespace("knitr", quietly = TRUE)) {
  knitr::opts_chunk$set(
    echo = TRUE,
    warning = FALSE,
    message = FALSE
  )
}

load_if_available <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      sprintf(
        "Package '%s' is required but not installed. Install it before rendering.",
        pkg
      ),
      call. = FALSE
    )
  }

  suppressPackageStartupMessages(
    library(pkg, character.only = TRUE)
  )
}

get_script_path <- function() {
  frame_files <- vapply(sys.frames(), function(x) {
    if (!is.null(x$ofile)) x$ofile else NA_character_
  }, character(1))

  script_path <- tail(na.omit(frame_files), 1)

  if (length(script_path) == 0) {
    stop("Could not determine the path to globalSetup.R.", call. = FALSE)
  }

  normalizePath(script_path, winslash = "/", mustWork = TRUE)
}

setup_script <- get_script_path()
classwork_dir <- dirname(setup_script)
project_root <- dirname(classwork_dir)
assets_dir <- file.path(classwork_dir, "assets")
fasta_assets_dir <- file.path(assets_dir, "secuencias_fasta")
genbank_assets_dir <- file.path(assets_dir, "secuencias_gb")
wuhan_reference_file <- file.path(assets_dir, "wuhan.fasta")

if (!dir.exists(assets_dir)) {
  stop("The assets directory was not found at 'classwork/assets'.", call. = FALSE)
}

if (requireNamespace("knitr", quietly = TRUE)) {
  knitr::opts_knit$set(root.dir = project_root)
}

required_packages <- c(
  "seqinr",
  "ggplot2",
  "dplyr",
  "tidyr",
  "ape",
  "Biostrings",
  "DECIPHER",
  "msa"
)

invisible(lapply(required_packages, load_if_available))

optional_packages <- c(
  "DBI",
  "RSQLite",
  "adegenet",
  "ggtree",
  "viridis",
  "ggmsa"
)

optional_package_status <- vapply(
  optional_packages,
  requireNamespace,
  logical(1),
  quietly = TRUE
)

list_project_fasta_files <- function(include_combined = FALSE) {
  if (!dir.exists(fasta_assets_dir)) {
    stop(
      "The FASTA directory was not found at 'classwork/assets/secuencias_fasta'.",
      call. = FALSE
    )
  }

  fasta_files <- list.files(
    path = fasta_assets_dir,
    pattern = "\\.fasta$",
    full.names = TRUE
  )

  if (!include_combined) {
    fasta_files <- fasta_files[basename(fasta_files) != "all_arn_sequences.fasta"]
  }

  fasta_files
}

project_fasta_files <- list_project_fasta_files()
combined_fasta_file <- file.path(fasta_assets_dir, "all_arn_sequences.fasta")

if (length(project_fasta_files) == 0) {
  stop(
    "No individual FASTA files were found in 'classwork/assets/secuencias_fasta'.",
    call. = FALSE
  )
}

if (!file.exists(wuhan_reference_file)) {
  warning(
    "Reference file 'classwork/assets/wuhan.fasta' was not found.",
    call. = FALSE
  )
}

read_project_fastas <- function(include_combined = FALSE) {
  fasta_files <- list_project_fasta_files(include_combined = include_combined)
  dna_set <- Biostrings::readDNAStringSet(fasta_files)

  if (length(dna_set) == length(fasta_files)) {
    names(dna_set) <- tools::file_path_sans_ext(basename(fasta_files))
  }

  dna_set
}
