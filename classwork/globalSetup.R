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
assets_dir <- file.path(classwork_dir, "assets")
sars_assets_dir <- file.path(assets_dir, "secuencias_sars")

if (!dir.exists(assets_dir)) {
  stop("The assets directory was not found at 'classwork/assets'.", call. = FALSE)
}

if (requireNamespace("knitr", quietly = TRUE)) {
  knitr::opts_knit$set(root.dir = assets_dir)
}

required_packages <- c(
  "seqinr",
  "ggplot2",
  "dplyr",
  "tidyr",
  "ape"
)

invisible(lapply(required_packages, load_if_available))

optional_packages <- c(
  "DECIPHER",
  "Biostrings",
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

required_files <- c("classwork/assets/secuencias_fasta/mexico.fasta", "classwork/assets/wuhan.fasta", "classwork/assets/secuencias_fasta/francia.fasta")
has_sequence_file <- function(filename) {
  file.exists(file.path(assets_dir, filename)) ||
    file.exists(file.path(sars_assets_dir, filename))
}

missing_files <- required_files[!vapply(required_files, has_sequence_file, logical(1))]

if (length(missing_files) > 0) {
  stop(
    sprintf(
      "Missing required FASTA files in classwork/assets: %s",
      paste(missing_files, collapse = ", ")
    ),
    call. = FALSE
  )
}
