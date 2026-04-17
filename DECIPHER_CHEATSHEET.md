# DECIPHER in VS Code - Quick Reference Cheat Sheet

## Installation (One-Time Setup)

```r
# Install from Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DECIPHER")

library(DECIPHER)  # Load the library
```

## Reading Sequences

```r
library(Biostrings)

# From FASTA file
sequences <- readDNAStringSet("path/to/file.fasta")

# Multiple files
seq1 <- readDNAStringSet("file1.fasta")
seq2 <- readDNAStringSet("file2.fasta")
all_seqs <- c(seq1, seq2)

# From GenBank
sequences <- readBStringSet("file.gb")  # For GenBank files
```

## Alignment

```r
# Basic alignment
aligned <- AlignSeqs(sequences, verbose = TRUE)

# Optimized alignment (longer but better)
aligned <- AlignSeqs(
  sequences,
  iterations = 2,
  gapOpening = -25,
  gapExtension = -5
)

# View alignment
print(aligned)
writeXStringSet(aligned, "output.fasta")
```

## Sequence Analysis

```r
# Get properties
width(sequences)                    # Lengths
names(sequences)                    # Names
length(sequences)                   # Number of sequences

# Base composition
letterFrequency(sequences, "ACGT")

# GC content
gc <- letterFrequency(sequences, "GC")
gc_percent <- rowSums(gc) / width(sequences) * 100

# Reverse complement
revcomp <- reverseComplement(sequences)

# Edit distance (Hamming distance)
dist <- dist.hamming(aligned)
```

## Tree Building

```r
library(ape)

# From aligned sequences
dist_mat <- as.matrix(dist.hamming(aligned))

# Build tree
nj_tree <- nj(as.dist(dist_mat))
upgma_tree <- upgma(dist_mat)

# Plot tree
plot(nj_tree)
title("Phylogenetic Tree")
axisPhylo()  # Add scale bar
```

## Visualization

```r
library(ggplot2)
library(ggtree)
library(ggmsa)

# Tree with ggtree
library(ggtree)
ggtree(tree) +
  geom_tiplab() +
  theme_tree2()

# Alignment visualization (write to FASTA first)
writeXStringSet(aligned, "aligned.fasta")
ggmsa("aligned.fasta", start = 1, end = 500)
```

## Pattern Searching

```r
# Find pattern
matches <- vmatchPattern("ATG", sequences)

# Count pattern
counts <- vcountPattern("ATG", sequences)

# Enumerate all patterns
freq <- oligonucleotideFrequency(sequences, width = 3)
```

## Common Workflows

### Workflow 1: Basic Variant Comparison

```r
library(DECIPHER)
library(Biostrings)
library(ape)

# Load
seqs <- readDNAStringSet("variants.fasta")

# Align
aligned <- AlignSeqs(seqs, verbose = TRUE)

# Distance & Tree
dist <- dist.hamming(aligned)
tree <- nj(as.dist(dist))

# Visualize
plot(tree)
```

### Workflow 2: Large Batch Analysis

```r
# Process multiple files in loop
fasta_files <- list.files("data/", pattern = "\\.fasta$", full.names = TRUE)

for (file in fasta_files) {
  seqs <- readDNAStringSet(file)
  aligned <- AlignSeqs(seqs)

  # Save results
  base_name <- sub(".fasta", "", basename(file))
  saveRDS(aligned, paste0("results/", base_name, "_aligned.rds"))
}
```

### Workflow 3: Taxonomic Classification

```r
library(DECIPHER)

# Classify sequences (requires trained classifier)
# ids <- IdTaxa(sequences, classifier, threshold = 50)
# For your project, focus on sequence homology instead
```

## Performance Tips

```r
# Use multiple cores
aligned <- AlignSeqs(sequences, processors = 4)

# For huge files, align in batches
big_seqs <- readDNAStringSet("huge_file.fasta")
batch_size <- 50
for (i in seq(1, length(big_seqs), by = batch_size)) {
  end <- min(i + batch_size - 1, length(big_seqs))
  batch <- big_seqs[i:end]
  aligned_batch <- AlignSeqs(batch)
  # Process batch
}

# Cache results
saveRDS(aligned, "cache/aligned.rds")
loaded <- readRDS("cache/aligned.rds")
```

## Troubleshooting

```r
# Check installation
library(DECIPHER)

# Check if sequences are valid
class(sequences)  # Should be "DNAStringSet"
length(sequences)
width(sequences)

# Clean ambiguous bases
cleaned <- chartr("NRYSWKMBDHV", "N", sequences)

# Memory check
gc()
object.size(sequences)

# Check for errors
tryCatch({
  aligned <- AlignSeqs(sequences)
}, error = function(e) {
  print(paste("Error:", e$message))
})
```

## VS Code Keyboard Shortcuts

| Action           | Shortcut (Mac)      |
| ---------------- | ------------------- |
| Run line in R    | Cmd + Enter         |
| Run chunk        | Cmd + Shift + Enter |
| Run all chunks   | Cmd + Shift + R     |
| Send to terminal | Ctrl + Alt + Enter  |
| Open terminal    | Ctrl + `            |

## Required Packages Summary

```r
# Required
library(DECIPHER)       # Sequence alignment & analysis
library(Biostrings)     # Sequence objects
library(ape)            # Phylogenetic trees

# Recommended for visualization
library(ggtree)         # Tree visualization
library(ggmsa)          # Alignment visualization
library(ggplot2)        # Plotting
library(seqinr)         # FASTA I/O (alternative)
```

## Your Project Structure

```
/Users/binivazquez/UniWorkspace/computational_biology/
├── classwork/
│   ├── globalSetup.R              ← Shared setup
│   ├── graficos.R                 ← Graphics & DECIPHER work
│   ├── basicsDNAan.Rmd            ← Analysis notebook
│   └── assets/
│       └── secuencias_fasta/      ← FASTA files
├── assignments/
│   └── analisis1_sp.R             ← Analysis script
├── DECIPHER_VS_CODE_GUIDE.md      ← Full guide
└── DECIPHER_CHEATSHEET.md         ← This file
```

## One-Liners

```r
# Load, align, distance, tree - all at once
plot(nj(as.dist(dist.hamming(AlignSeqs(readDNAStringSet("seqs.fasta"))))))

# Check what's in a sequence file
readDNAStringSet("file.fasta")

# Quick GC content
colMeans(letterFrequency(readDNAStringSet("seqs.fasta"), "GC")) / width(readDNAStringSet("seqs.fasta"))

# Find all "ATG" starts
vmatchPattern("ATG", readDNAStringSet("seqs.fasta"))
```

## When to Use DECIPHER vs Alternatives

| Task          | Tool    | DECIPHER         | Biostrings | seqinr      |
| ------------- | ------- | ---------------- | ---------- | ----------- |
| Read FASTA    | Yes     | readDNAStringSet | Yes        | read.fasta  |
| Alignment     | Yes     | **Best**         | No         | No          |
| Distance      | Yes     | dist.hamming     | No         | dist.seqinr |
| Trees         | No      | No               | Use `ape`  | Use `ape`   |
| Visualization | Partial | With ggtree      | With ggmsa | With seqinr |

---

**Keep this next to you while coding!** Print and post on monitor if desired.

For detailed explanations, see `DECIPHER_VS_CODE_GUIDE.md`
