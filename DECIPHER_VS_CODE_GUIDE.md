# Complete Guide to DECIPHER in VS Code

## A Comprehensive Resource for Computational Biology Work

---

## Table of Contents

1. [Introduction](#introduction)
2. [Installation & Setup](#installation--setup)
3. [Core DECIPHER Functions](#core-decipher-functions)
4. [Working with DNA Sequences](#working-with-dna-sequences)
5. [Sequence Alignment](#sequence-alignment)
6. [Phylogenetic Analysis](#phylogenetic-analysis)
7. [VS Code Integration & Workflow](#vs-code-integration--workflow)
8. [Best Practices](#best-practices)
9. [Real Project Examples](#real-project-examples)
10. [Troubleshooting](#troubleshooting)
11. [Advanced Techniques](#advanced-techniques)
12. [Performance Optimization](#performance-optimization)

---

## Introduction

**DECIPHER** (Detect and Correct Errors in RNA and Protein) is a Bioconductor package for analyzing and manipulating DNA, RNA, and protein sequences. It's particularly powerful for:

- **Sequence alignment** (including gap opening/closing costs optimization)
- **Taxonomy classification** (via IDTAXA)
- **Polymorphism detection**
- **Community diversity analysis**
- **Sequence comparison and visualization**

### Why DECIPHER in your workflow?

Your project uses DECIPHER alongside:

- **seqinr** - Reading/writing FASTA files
- **ape** - Phylogenetic analysis
- **ggtree** - Advanced phylogenetic tree visualization
- **ggmsa** - Multiple sequence alignment visualization

---

## Installation & Setup

### 1. Installing DECIPHER

DECIPHER is available through **Bioconductor**, not CRAN:

```r
# Method 1: Using BiocManager (recommended)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DECIPHER")

# Method 2: Direct from Bioconductor
source("https://bioconductor.org/biocLite.R")
biocLite("DECIPHER")
```

### 2. Verifying Installation

```r
# Check if DECIPHER is installed
library(DECIPHER)

# Check version
packageVersion("DECIPHER")

# View citation
citation("DECIPHER")
```

### 3. VS Code Configuration

Create or update your `globalSetup.R` to handle DECIPHER gracefully:

```r
# Load DECIPHER with error handling
load_decipher <- function() {
  if (!requireNamespace("DECIPHER", quietly = TRUE)) {
    message("DECIPHER not found. Installing from Bioconductor...")
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("DECIPHER")
  }
  library(DECIPHER)
}

# Call this in your setup script
load_decipher()
```

---

## Core DECIPHER Functions

### Function Reference Table

| Function                    | Purpose                          | Input               | Output                 |
| --------------------------- | -------------------------------- | ------------------- | ---------------------- |
| `AlignSeqs()`               | Multiple sequence alignment      | DNAStringSet        | Aligned DNAStringSet   |
| `OrientNucleotides()`       | Orient sequences to reference    | DNAStringSet        | Oriented sequences     |
| `IdTaxa()`                  | Taxonomic classification         | DNAStringSet        | Classification results |
| `TreeLine()`                | Phylogenetic tree from alignment | DNAStringSet        | Dendrogram             |
| `Levenshtein()`             | Edit distance between sequences  | 2 character vectors | Numeric distance       |
| `AmplificationEfficiency()` | PCR efficiency prediction        | DNAString sequence  | Efficiency data        |
| `DesignProbes()`            | Probe design for sequences       | DNAStringSet        | Probe sequences        |
| `CountPatterns()`           | Count motif patterns             | DNAStringSet        | Pattern counts         |

---

## Working with DNA Sequences

### 1. Reading Sequences with DECIPHER

```r
library(DECIPHER)
library(Biostrings)

# Method 1: Using readDNAStringSet (DECIPHER-compatible)
sequences <- readDNAStringSet("path/to/sequences.fasta")

# Method 2: Using seqinr (then convert)
seq_list <- read.fasta("sequences.fasta", as.string = TRUE, forceDNAtolower = FALSE)
sequences <- DNAStringSet(seq_list)

# Method 3: From your project structure
mexico_seq <- readDNAStringSet("classwork/assets/secuencias_fasta/mexico.fasta")
wuhan_seq <- readDNAStringSet("classwork/assets/wuhan.fasta")
francia_seq <- readDNAStringSet("classwork/assets/secuencias_fasta/francia.fasta")
```

### 2. Basic Sequence Inspection

```r
# Get sequence information
length(mexico_seq)              # Number of sequences
width(mexico_seq)              # Sequence lengths
names(mexico_seq)              # Sequence names

# Calculate GC content
gc_content <- function(seq) {
  g_count <- letterFrequency(seq, "G")
  c_count <- letterFrequency(seq, "C")
  (g_count + c_count) / width(seq)
}

gc_mexico <- gc_content(mexico_seq[[1]])
gc_wuhan <- gc_content(wuhan_seq[[1]])
gc_francia <- gc_content(francia_seq[[1]])

# View sequence composition
letterFrequency(mexico_seq[[1]], letters = "ACGT")
```

### 3. Sequence Manipulation

```r
# Reverse complement
rev_comp_seq <- reverseComplement(mexico_seq)

# Extract subsequence
subseq <- subseq(mexico_seq[[1]], start = 100, end = 200)

# Replace ambiguous nucleotides
clean_seq <- chartr("RYSWKMBDHV", "AG", mexico_seq)

# Get unique sequences
unique_seqs <- unique(mexico_seq)
```

---

## Sequence Alignment

### 1. Multiple Sequence Alignment (MSA)

```r
library(DECIPHER)

# Combine sequences into one object
all_sequences <- c(mexico_seq, wuhan_seq, francia_seq)

# Perform alignment
aligned_sequences <- AlignSeqs(all_sequences,
                              verbose = TRUE,
                              iterations = 2)

# View alignment
print(aligned_sequences)
```

### 2. Alignment Parameters

```r
# Fine-tuning alignment (advanced)
aligned <- AlignSeqs(
  all_sequences,
  gapOpening = -25,          # Cost of opening a gap (more negative = costlier)
  gapExtension = -5,         # Cost of extending a gap
  useQuality = FALSE,        # Don't use quality scores
  verbose = TRUE,
  iterations = 3,            # More iterations = better alignment but slower
  anchor = NA,               # Don't use anchor sequences
  processors = NULL          # Use all available cores
)

# Calculate pairwise similarity
similarity <- AlignSeqs(
  all_sequences,
  returnAlignmentOnly = FALSE  # Returns alignment quality metrics too
)
```

### 3. Visualizing Alignments

```r
library(ggmsa)

# Convert to format compatible with ggmsa
fasta_file <- tempfile(fileext = ".fasta")
writeXStringSet(aligned_sequences, filepath = fasta_file)
aligned_df <- read.fasta(fasta_file, as.string = TRUE)

# Create visualization
ggmsa(fasta_file, start = 1, end = 500, color = "Clustal") +
  theme_minimal() +
  labs(title = "Multiple Sequence Alignment: Mexico, Wuhan, Francia")
```

---

## Phylogenetic Analysis

### 1. Build Trees from Alignment

```r
library(ape)

# Using DECIPHER alignment
aligned <- AlignSeqs(all_sequences)

# Convert to distance matrix
dist_matrix <- as.matrix(dist.hamming(aligned))

# Create phylogenetic tree (UPGMA method)
tree <- upgma(dist_matrix)

# Or use neighbor-joining
nj_tree <- nj(dist_matrix)

# Plot tree
plot(tree)
title("Phylogenetic Tree: SARS-CoV-2 Variants")

# Add bootstrap support
plot(nj_tree)
axisPhylo()
```

### 2. Advanced Tree Visualization with ggtree

```r
library(ggtree)

# Convert tree to compatible format
ggtree_obj <- ggtree(tree)

# Create visualization
ggtree_obj +
  geom_tiplab() +
  theme_tree2() +
  labs(title = "SARS-CoV-2 Phylogenetic Relationships")

# With node labels
ggtree(tree) +
  geom_tiplab() +
  geom_nodelab(aes(label = node)) +
  theme_tree()
```

---

## VS Code Integration & Workflow

### 1. Setting Up Your VS Code Environment

**Recommended Extensions:**

- `R` (Positron or R Client)
- `Git Lens` - Track code history
- `Markdown All in One` - For this guide and documentation
- `vscode-pdf` - View PDF outputs

**VS Code Settings (`settings.json`):**

```json
{
  "[r]": {
    "editor.defaultFormatter": "REditorSupport.r",
    "editor.formatOnSave": true,
    "editor.wordWrap": "on"
  },
  "r.sessionWatcher": true,
  "r.plot.useHttpgd": true,
  "r.rterm.mac": "/usr/local/bin/R"
}
```

### 2. Project Structure Best Practice

```
computational_biology/
├── classwork/
│   ├── globalSetup.R              # Shared configuration
│   ├── basicsDNAan.Rmd            # Analysis notebooks
│   ├── graficos.R                 # DECIPHER graphics work
│   ├── assets/
│   │   ├── wuhan.fasta
│   │   ├── secuencias_fasta/      # FASTA sequences
│   │   └── secuencias_gen/        # GenBank files
│   └── outputs/                   # Results directory
├── assignments/
│   ├── analisis1_sp.R             # Sequence analysis
│   └── tarea2.Rmd                 # Assignments
└── DECIPHER_VS_CODE_GUIDE.md      # This file!
```

### 3. Executing Code in VS Code

**Option A: Interactive R Terminal**

```
Ctrl+Shift+Enter  # Send line to R terminal
Cmd+Shift+Enter   # Send selection
```

**Option B: R Markdown in VS Code**

```
Ctrl+Alt+R        # Run current chunk
Ctrl+Shift+R      # Run all chunks
```

**Option C: Source entire script**

```r
source("classwork/graficos.R")
```

---

## Best Practices

### 1. Memory Management

```r
# For large FASTA files, read in chunks
process_large_fasta <- function(file_path, chunk_size = 100) {
  all_seqs <- readDNAStringSet(file_path)

  # Process in chunks to avoid memory overload
  for (i in seq(1, length(all_seqs), by = chunk_size)) {
    chunk <- all_seqs[i:min(i + chunk_size - 1, length(all_seqs))]
    # Process chunk here
    aligned <- AlignSeqs(chunk)
  }
}

# Clear unused variables
rm(old_alignment)
gc()  # Force garbage collection
```

### 2. Error Handling

```r
# Wrap DECIPHER operations with error handling
safe_align <- function(sequences) {
  tryCatch({
    AlignSeqs(sequences, verbose = TRUE)
  }, error = function(e) {
    message("Alignment failed: ", conditionMessage(e))
    NULL
  }, warning = function(w) {
    message("Warning during alignment: ", conditionMessage(w))
  })
}

# Check sequence validity before alignment
validate_sequences <- function(seqs) {
  # Check for empty sequences
  if (any(width(seqs) == 0)) {
    stop("Found empty sequences")
  }

  # Check for non-DNA characters
  valid_chars <- grepl("^[ACGTN-]*$", as.character(seqs))
  if (!all(valid_chars)) {
    warning("Found non-standard nucleotides")
  }

  return(TRUE)
}
```

### 3. Reproducibility

```r
# Always set seed for reproducible results
set.seed(42)

# Document versions
sessionInfo()

# Save your work with timestamps
save_analysis <- function(object, name) {
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  filename <- paste0("outputs/", name, "_", timestamp, ".RData")
  save(object, file = filename)
  message("Saved to: ", filename)
}
```

### 4. Code Organization

```r
# Organize functions at the top
source("classwork/globalSetup.R")

# Then analysis
mexico_seq <- readDNAStringSet("classwork/assets/secuencias_fasta/mexico.fasta")
wuhan_seq <- readDNAStringSet("classwork/assets/wuhan.fasta")

# Then visualization
aligned <- AlignSeqs(c(mexico_seq, wuhan_seq))
plot_alignment(aligned)
```

---

## Real Project Examples

### Example 1: Complete Variant Analysis (Your Project)

```r
library(DECIPHER)
library(Biostrings)
library(ape)
library(seqinr)
library(ggplot2)
library(dplyr)

# Setup
setwd("/Users/binivazquez/UniWorkspace/computational_biology")
source("classwork/globalSetup.R")

# 1. Load sequences
mexico_seq <- readDNAStringSet("classwork/assets/secuencias_fasta/mexico.fasta")
wuhan_seq <- readDNAStringSet("classwork/assets/wuhan.fasta")
francia_seq <- readDNAStringSet("classwork/assets/secuencias_fasta/francia.fasta")

# 2. Create sequence set
all_seqs <- c(mexico_seq, wuhan_seq, francia_seq)
names(all_seqs) <- c("Mexico", "Wuhan", "Francia")

# 3. Align sequences
aligned <- AlignSeqs(all_seqs, verbose = TRUE)

# 4. Analyze alignment
mutations <- data.frame()
for (i in 1:(length(aligned)-1)) {
  for (j in (i+1):length(aligned)) {
    seq1 <- aligned[[i]]
    seq2 <- aligned[[j]]
    diffs <- seq1 != seq2

    mutations <- rbind(mutations, data.frame(
      Pair = paste(names(aligned)[i], "vs", names(aligned)[j]),
      Mutations = sum(diffs),
      MutationRate = sum(diffs) / width(seq1)
    ))
  }
}

print(mutations)

# 5. Build phylogenetic tree
dist_mat <- as.matrix(dist.hamming(aligned))
tree <- nj(as.dist(dist_mat))

# 6. Visualize results
par(mfrow = c(1, 2))

# Tree
plot(tree)
title("Phylogenetic Tree")

# Mutation heatmap
heatmap(dist_mat, main = "Sequence Distances")

par(mfrow = c(1, 1))
```

### Example 2: Automated Variant Classification

```r
classify_variants <- function(sequences_file) {
  library(DECIPHER)

  # Load reference genomes
  refs <- readDNAStringSet("classwork/assets/reference_variants.fasta")

  # Load query sequences
  queries <- readDNAStringSet(sequences_file)

  # Align each query to references
  results <- data.frame()

  for (i in seq_len(length(queries))) {
    query <- queries[i]

    # Compare to each reference
    for (j in seq_len(length(refs))) {
      ref <- refs[j]
      aligned <- AlignSeqs(c(query, ref))

      # Calculate similarity
      similarity <- sum(aligned[[1]] == aligned[[2]]) / width(aligned[[1]])

      results <- rbind(results, data.frame(
        Query = names(queries)[i],
        Reference = names(refs)[j],
        Similarity = similarity
      ))
    }
  }

  # Get best match for each query
  best_matches <- results %>%
    group_by(Query) %>%
    slice_max(Similarity) %>%
    ungroup()

  return(best_matches)
}

# Usage
variant_classification <- classify_variants("new_sequences.fasta")
write.csv(variant_classification, "outputs/variant_classification.csv")
```

### Example 3: Mutation Hotspot Identification

```r
find_mutation_hotspots <- function(aligned_sequences, window_size = 100) {
  library(DECIPHER)
  library(ggplot2)

  # Calculate position-wise variation
  seq_array <- as.character(aligned_sequences)
  n_seqs <- length(aligned_sequences)
  seq_length <- width(aligned_sequences[1])

  # Sliding window analysis
  mutation_scores <- numeric(seq_length - window_size)

  for (i in 1:(seq_length - window_size)) {
    window <- substr(seq_array, i, i + window_size - 1)

    # Count polymorphic sites in window
    polymorphic <- 0
    for (pos in 1:window_size) {
      bases <- substr(window, pos, pos)
      unique_bases <- n_distinct(unique(bases))
      if (unique_bases > 1) {
        polymorphic <- polymorphic + 1
      }
    }

    mutation_scores[i] <- polymorphic / window_size
  }

  # Visualize
  df <- data.frame(
    Position = 1:length(mutation_scores),
    MutationDensity = mutation_scores
  )

  ggplot(df, aes(x = Position, y = MutationDensity)) +
    geom_line(color = "steelblue") +
    geom_smooth(method = "loess", color = "red", alpha = 0.2) +
    theme_minimal() +
    labs(
      title = "Mutation Hotspots in Alignment",
      x = "Genome Position",
      y = "Polymorphic Sites Density"
    )
}

# Usage
aligned <- AlignSeqs(all_seqs)
hotspot_plot <- find_mutation_hotspots(aligned, window_size = 200)
print(hotspot_plot)
```

---

## Troubleshooting

### Common Issues & Solutions

#### Issue 1: "could not find function 'AlignSeqs'"

```r
# Solution: DECIPHER not loaded
library(DECIPHER)

# Or check if installed
if (!requireNamespace("DECIPHER", quietly = TRUE)) {
  BiocManager::install("DECIPHER")
}
```

#### Issue 2: "Biostrings object is invalid"

```r
# Solution: Input format mismatch
# Make sure sequences are DNAStringSet, not character vectors

# Wrong:
bad_seq <- "ACGTACGT"

# Right:
good_seq <- DNAString("ACGTACGT")
good_set <- DNAStringSet(c("ACGTACGT", "TGCATGCA"))

# Or read from file:
seq_set <- readDNAStringSet("sequences.fasta")
```

#### Issue 3: Alignment too slow

```r
# Solution: Use multiprocessing and optimize parameters
aligned <- AlignSeqs(
  large_sequence_set,
  processors = 4,        # Use 4 cores
  iterations = 1,        # Reduce iterations
  verbose = TRUE         # Monitor progress
)

# Or align in batches
batch_align <- function(sequences, batch_size = 50) {
  n_batches <- ceiling(length(sequences) / batch_size)
  all_aligned <- list()

  for (i in 1:n_batches) {
    start_idx <- (i-1) * batch_size + 1
    end_idx <- min(i * batch_size, length(sequences))

    batch <- sequences[start_idx:end_idx]
    all_aligned[[i]] <- AlignSeqs(batch, verbose = FALSE)

    message(paste("Completed batch", i, "of", n_batches))
  }

  # Combine results
  do.call(c, all_aligned)
}
```

#### Issue 4: Memory overflow on large files

```r
# Solution: Stream processing
process_fasta_streaming <- function(fasta_file, chunk_size = 100) {
  # Read in chunks
  stream <- open(fasta_file)

  repeat {
    chunk <- readLines(stream, n = chunk_size * 2)  # FASTA lines
    if (length(chunk) == 0) break

    # Convert to sequences and process
    temp_file <- tempfile(fileext = ".fasta")
    writeLines(chunk, temp_file)

    sequences <- readDNAStringSet(temp_file)
    aligned <- AlignSeqs(sequences)

    # Save results
    saveRDS(aligned, paste0("chunk_", format(Sys.time(), "%s"), ".rds"))

    unlink(temp_file)
  }

  close(stream)
}
```

#### Issue 5: "No valid letters" warning

```r
# Solution: Clean sequences of ambiguous nucleotides
clean_sequences <- function(sequences) {
  # Replace ambiguous codes with N
  chartr("RYSWKMBDHV", "N", sequences)
}

# Or be more specific
strict_clean <- function(sequences) {
  # Keep only ACGT
  cleaned <- DNAStringSet(
    chartr("NRYSWKMBDHV", "A", as.character(sequences))
  )
  return(cleaned)
}

# Check before processing
invalid_chars <- function(sequences) {
  seq_str <- as.character(sequences)
  invalid <- grepl("[^ACGTN-]", seq_str)
  return(seq_str[invalid])
}
```

---

## Advanced Techniques

### 1. Custom Scoring Matrices

```r
# Create custom alignment scoring
custom_align <- function(sequences) {
  # Define scoring matrix
  score_matrix <- matrix(c(
    5, -4, -4, -4,   # A matches
    -4, 5, -4, -4,   # C matches
    -4, -4, 5, -4,   # G matches
    -4, -4, -4, 5    # T matches
  ), nrow = 4)

  rownames(score_matrix) <- c("A", "C", "G", "T")
  colnames(score_matrix) <- c("A", "C", "G", "T")

  # DECIPHER alignment with custom parameters
  AlignSeqs(
    sequences,
    gapOpening = -15,
    gapExtension = -5
  )
}
```

### 2. Taxonomic Classification

```r
# IDTAXA function for rapid taxonomic classification
classify_sequences <- function(sequences, training_set) {
  library(DECIPHER)

  # Load or create training set
  # training_set <- learnTaxa(reference_sequences)

  # Classify unknown sequences
  ids <- IdTaxa(
    sequences,
    training_set,
    strand = "top",
    threshold = 50
  )

  # Extract results
  for (i in seq_along(ids)) {
    cat(paste(names(sequences)[i], "\n"))
    print(ids[[i]])
  }

  return(ids)
}
```

### 3. Motif Analysis

```r
# Find and analyze motifs
analyze_motifs <- function(sequences, motif_pattern) {
  library(DECIPHER)

  # Count occurrences
  matches <- vmatchPattern(motif_pattern, sequences)

  # Extract statistics
  motif_stats <- data.frame(
    Sequence = names(sequences),
    MotifCount = sapply(matches, length),
    SequenceLength = width(sequences)
  )

  motif_stats$MotifDensity <-
    motif_stats$MotifCount / motif_stats$SequenceLength

  return(motif_stats)
}

# Usage
orf1_motif <- analyze_motifs(sequences, "ATG")
```

### 4. Consensus Sequence Generation

```r
# Generate consensus from aligned sequences
get_consensus <- function(aligned_sequences) {
  library(DECIPHER)

  # Convert alignment to matrix
  aligned_matrix <- as.character(aligned_sequences)

  consensus <- character(width(aligned_sequences[1]))

  for (i in 1:width(aligned_sequences[1])) {
    bases_at_pos <- substr(aligned_matrix, i, i)
    # Get most common base
    base_table <- table(bases_at_pos)
    consensus[i] <- names(base_table)[which.max(base_table)]
  }

  DNAString(paste(consensus, collapse = ""))
}
```

---

## Performance Optimization

### 1. Parallel Processing

```r
# Enable multicore processing
library(parallel)

# Auto-detect cores
n_cores <- detectCores() - 1  # Leave 1 core free

# Use in DECIPHER
aligned <- AlignSeqs(
  large_sequence_set,
  processors = n_cores,
  verbose = TRUE
)

# Manual parallel processing
run_analysis_parallel <- function(sequences, n_cores = 4) {
  library(parallel)

  # Split sequences
  seq_chunks <- split(sequences,
                      rep(1:n_cores,
                          length.out = length(sequences)))

  # Process in parallel
  results <- mclapply(seq_chunks, AlignSeqs, mc.cores = n_cores)

  # Combine results
  combined <- do.call(c, results)
  return(combined)
}
```

### 2. Caching Results

```r
# Cache alignment results to avoid recalculation
cache_alignment <- function(sequences, cache_dir = "cache/") {
  library(digest)

  # Create unique hash for this sequence set
  seq_hash <- digest(as.character(sequences), algo = "md5")
  cache_file <- file.path(cache_dir, paste0(seq_hash, ".rds"))

  if (file.exists(cache_file)) {
    message("Loading cached alignment...")
    return(readRDS(cache_file))
  }

  message("Computing alignment...")
  aligned <- AlignSeqs(sequences, verbose = TRUE)

  dir.create(cache_dir, showWarnings = FALSE)
  saveRDS(aligned, cache_file)

  return(aligned)
}
```

### 3. Progress Monitoring

```r
# Track long operations
monitor_alignment <- function(sequences) {
  library(progress)

  # Create progress bar
  pb <- progress_bar$new(
    format = "Aligning [:bar] :percent eta: :eta",
    total = length(sequences),
    clear = FALSE,
    width = 60
  )

  pb$tick(0)

  aligned <- AlignSeqs(
    sequences,
    verbose = TRUE,
    processors = 1
  )

  pb$tick(length(sequences))

  return(aligned)
}
```

### 4. Memory Profiling

```r
# Profile memory usage
profile_operation <- function(sequences) {
  library(profvis)

  profvis({
    aligned <- AlignSeqs(sequences)
    tree <- nj(as.dist(dist.hamming(aligned)))
  })
}

# Simpler memory check
check_memory <- function(sequences) {
  cat("Object sizes:\n")
  print(object.size(sequences), units = "MB")

  cat("\nAvailable RAM:\n")
  gc()  # Compact memory
}
```

---

## Quick Reference

### One-Liners for Common Tasks

```r
# Load and align
aligned <- AlignSeqs(readDNAStringSet("seqs.fasta"))

# Calculate sequence distance
dist <- dist.hamming(aligned)

# Build tree
tree <- nj(as.dist(dist))

# Plot tree
plot(tree); axisPhylo()

# Get GC content
gc <- rowSums(alphabetFrequency(sequences)[,1:2]) / width(sequences) * 100

# Find reverse complement
revcomp <- reverseComplement(sequences)

# Count patterns
counts <- vcountPattern("ATG", sequences)

# Extract alignment statistics
stats <- as.data.frame(letterFrequency(sequences, "ACGT", as.prob=TRUE))
```

---

## Resources & References

### Official Documentation

- **DECIPHER Manual**: https://www.bioconductor.org/packages/release/bioc/html/DECIPHER.html
- **Biostrings Vignette**: https://bioconductor.org/packages/release/bioc/vignettes/Biostrings/inst/doc/BiostringsQuickOverview.pdf
- **ape Package**: https://cran.r-project.org/web/packages/ape/

### Helpful Papers

- Wright et al. (2016) - DECIPHER: A Bioconductor package for sequence alignment, taxonomy classification, and phylogenetic reconstruction
- Schliep et al. (2017) - phangorn: phylogenetic analysis in R

### Community Resources

- Bioconductor Support: https://support.bioconductor.org/
- Stack Overflow: Tag `[r]` and `[bioconductor]`
- GitHub: https://github.com/grunwald/DECIPHER

---

## Tips for Effective VS Code + DECIPHER Workflow

### 1. Use R Markdown for Documentation

Create `.Rmd` files to combine code with markdown notes:

````r
---
title: "SARS-CoV-2 Variant Analysis"
output: html_document
---

## Introduction
Brief explanation of analysis...

## Methods
### Data Loading
```{r load-data}
sequences <- readDNAStringSet("...")
````

### Analysis

```{r analysis}
aligned <- AlignSeqs(sequences)
```

````

### 2. Create Reusable Function Libraries

```r
# functions/sequence_utils.R
align_and_analyze <- function(fasta_file) {
  seqs <- readDNAStringSet(fasta_file)
  aligned <- AlignSeqs(seqs)
  dist_mat <- as.matrix(dist.hamming(aligned))
  tree <- nj(as.dist(dist_mat))

  list(
    alignment = aligned,
    distance = dist_mat,
    tree = tree
  )
}

# Then source in your analysis:
source("functions/sequence_utils.R")
````

### 3. Version Control Your Analysis

```bash
git add *.R *.Rmd
git commit -m "Add DECIPHER-based variant analysis"
git push
```

### 4. Save Intermediate Results

```r
# Save alignments to avoid recalculation
saveRDS(aligned, "cache/alignment_2024.rds")

# Load later
aligned <- readRDS("cache/alignment_2024.rds")
```

---

## Final Checklist for DECIPHER Projects

- [ ] DECIPHER installed via BiocManager
- [ ] `globalSetup.R` configured with error handling
- [ ] Project structured with `classwork/`, `assignments/`, `outputs/` folders
- [ ] README documenting data sources
- [ ] .gitignore configured for large FASTA files
- [ ] Scripts use absolute/relative paths correctly
- [ ] Comments explaining key biological interpretations
- [ ] Results saved with timestamps
- [ ] Figures exported at publication quality (300 dpi)
- [ ] Code reproducible with `sessionInfo()` documented

---

**Last Updated**: April 17, 2026
**For your team's computational biology work**

Happy sequencing! 🧬
