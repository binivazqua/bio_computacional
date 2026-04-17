# DECIPHER FAQ & Common Solutions

## Installation & Setup

### Q: DECIPHER installation fails with "BiocManager not found"

**A:** You need to install BiocManager first:

```r
install.packages("BiocManager")
BiocManager::install("DECIPHER")
```

### Q: "Package 'DECIPHER' was built under R version X.X.X"

**A:** This is just a warning. It works fine. If you get actual errors:

```r
update.packages()  # Update R packages to match your R version
```

### Q: How do I check which version of DECIPHER I have?

**A:**

```r
library(DECIPHER)
packageVersion("DECIPHER")
citation("DECIPHER")  # Also see how to cite in publications
```

---

## Reading & Loading Sequences

### Q: Should I use `readDNAStringSet()` or `read.fasta()`?

**A:**

- **`readDNAStringSet()`** (from Biostrings) - For DECIPHER work ✅
- **`read.fasta()`** (from seqinr) - For basic I/O, then convert to DNAStringSet

Use this to convert seqinr output:

```r
# seqinr format
seq_list <- read.fasta("file.fasta", as.string = TRUE)

# Convert to DECIPHER format
seq_set <- DNAStringSet(unlist(seq_list))
```

### Q: My FASTA file has multiple sequences per entry - is that a problem?

**A:** No, DECIPHER handles it:

```r
sequences <- readDNAStringSet("file.fasta")
length(sequences)  # Shows number of sequences

# If you want to ensure they're separate
sequences <- readDNAStringSet("file.fasta", use.names = TRUE)
```

### Q: How do I load GenBank (.gb or .gbk) files?

**A:** DECIPHER supports GenBank through Biostrings:

```r
library(Biostrings)
sequences <- readBStringSet("file.gb")  # Generic biological sequences

# Or if you just want the sequences as DNA:
library(ape)
seq <- read.GenBank(accession_numbers)
```

### Q: Can I read from a URL directly?

**A:** Yes!

```r
url <- "https://example.com/sequences.fasta"
sequences <- readDNAStringSet(url)
```

---

## Alignment Issues

### Q: Alignment is running very slowly. Can I speed it up?

**A:** Yes! Several options:

```r
# Option 1: Use parallel processing
aligned <- AlignSeqs(sequences, processors = 4)  # Use 4 CPU cores

# Option 2: Reduce iterations (faster but less accurate)
aligned <- AlignSeqs(sequences, iterations = 1)  # Default is 2

# Option 3: Align in smaller batches
# (See Performance section in main guide)

# Option 4: Check available cores first
detectCores()  # See how many cores your Mac has
# Use one less than total to keep system responsive
```

### Q: "Error: could not find function 'AlignSeqs'"

**A:** DECIPHER is not loaded:

```r
library(DECIPHER)

# Or check if it's installed
library(DECIPHER, verbose = TRUE)  # Shows loading details
```

### Q: Alignment output looks wrong - sequences have gaps everywhere

**A:** Possible causes:

```r
# 1. Check if sequences are already aligned
width(sequences)  # All should be same width if pre-aligned

# 2. Check sequence quality
# Make sure there are no invalid characters
invalid_chars <- sequences[!grepl("^[ACGTN-]*$", as.character(sequences))]

# 3. Re-run without iterations (might help):
aligned <- AlignSeqs(sequences, iterations = 1)

# 4. Try different gap penalties
aligned <- AlignSeqs(sequences, gapOpening = -50, gapExtension = -10)
```

### Q: Can I use a reference sequence for alignment?

**A:** DECIPHER doesn't support this directly, but you can:

```r
# Include reference as first sequence
reference <- readDNAStringSet("reference.fasta")
queries <- readDNAStringSet("queries.fasta")
all_seqs <- c(reference, queries)

aligned <- AlignSeqs(all_seqs)

# Reference will be first in alignment
```

### Q: How do I view my alignment after AlignSeqs()?

**A:**

```r
# Method 1: Print (shows in console)
print(aligned)

# Method 2: Save and view as FASTA
writeXStringSet(aligned, "alignment.fasta")

# Method 3: View in ggmsa (graphical)
writeXStringSet(aligned, "temp.fasta")
library(ggmsa)
ggmsa("temp.fasta", start = 1, end = 100)

# Method 4: Extract as character
aligned_char <- as.character(aligned)
```

---

## Phylogenetics & Distance

### Q: How do I calculate distance between aligned sequences?

**A:** Use `dist.hamming()` from ape:

```r
library(ape)

aligned <- AlignSeqs(sequences)
dist_matrix <- dist.hamming(aligned)  # Hamming distance
dist_matrix <- dist.dna(aligned, model = "JC69")  # Evolutionary model

# Convert to numeric matrix for visualization
dist_numeric <- as.matrix(dist_matrix)
```

### Q: What's the difference between Hamming and other distances?

**A:**

- **Hamming** - Simple: count differences
- **Jukes-Cantor (JC69)** - Accounts for multiple hits at same site
- **Kimura 2-parameter (K80)** - Distinguishes transitions/transversions
- **Tamura-Nei (TN93)** - Most complex, usually most accurate

```r
library(ape)
# Compare methods
hamming <- dist.hamming(aligned)
jc69 <- dist.dna(aligned, model = "JC69")
tn93 <- dist.dna(aligned, model = "TN93")
```

### Q: How do I build a phylogenetic tree from my alignment?

**A:**

```r
library(ape)

aligned <- AlignSeqs(sequences)
dist_mat <- dist.hamming(aligned)

# Method 1: Neighbor-Joining (most common)
tree <- nj(as.dist(dist_mat))

# Method 2: UPGMA (assumes clock)
tree <- upgma(dist_mat)

# Plot it
plot(tree)
axisPhylo()  # Add scale bar

# Get statistics
tree$tip.label   # Sequence names
tree$edge        # Tree structure
```

### Q: How do I add bootstrap support to my tree?

**A:** Use ape's `boot.phylo()`:

```r
library(ape)

# Define function for bootstrap
boot_func <- function(x) {
  d <- dist.hamming(x)
  nj(as.dist(d))
}

# Run bootstrap (100 replicates)
boot_trees <- boot.phylo(tree, aligned, boot_func, rooted = TRUE)

# Add to plot
plot(tree)
nodelabels(boot_trees, cex = 0.7)  # Shows bootstrap values at nodes
```

---

## Visualization Issues

### Q: How do I make better-looking trees with ggtree?

**A:**

```r
library(ggtree)

# Basic tree
ggtree(tree) +
  geom_tiplab() +
  theme_tree2()

# Enhanced tree with colors and shapes
ggtree(tree, layout = "rectangular") +
  geom_tiplab(aes(color = ifelse(grepl("mexico", label), "red", "blue"))) +
  geom_nodepoint(aes(color = bootstrap > 95), size = 3) +
  scale_color_manual(values = c("TRUE" = "green", "FALSE" = "gray")) +
  theme_tree2() +
  theme(legend.position = "right")
```

### Q: How do I visualize my multiple sequence alignment?

**A:**

```r
library(ggmsa)

# Save alignment to temp file
writeXStringSet(aligned, "temp_alignment.fasta")

# Visualize with ggmsa
p <- ggmsa("temp_alignment.fasta",
           start = 1,
           end = 200,
           char_width = 0.5,
           color = "Clustal")

p + theme_minimal()

# Clean up
unlink("temp_alignment.fasta")
```

### Q: How do I make a publication-quality figure?

**A:**

```r
# Save at high resolution
png("figures/tree.png", width = 3000, height = 2000, res = 300)

# Or PDF (better for publications)
pdf("figures/tree.pdf", width = 10, height = 10)

# Create plot
ggtree(tree) +
  geom_tiplab(size = 4, fontface = "italic") +
  geom_nodelab(aes(label = bootstrap), hjust = -0.3) +
  theme_tree2() +
  xlim(0, 0.1)

dev.off()  # Close and save

cat("Saved to: figures/tree.pdf\n")
```

---

## Performance & Large Files

### Q: My FASTA file is HUGE (>1GB). How do I work with it?

**A:** Process in chunks:

```r
process_large_fasta <- function(file_path) {
  # Read in batches
  con <- file(file_path, open = "r")

  batch_size <- 100  # Number of sequences per batch
  batch_num <- 1

  while (TRUE) {
    # Read batch
    seqs <- readDNAStringSet(con, nrec = batch_size)

    if (length(seqs) == 0) break  # End of file

    # Process batch
    aligned <- AlignSeqs(seqs, verbose = FALSE)

    # Save results
    saveRDS(aligned, paste0("aligned_batch_", batch_num, ".rds"))

    batch_num <- batch_num + 1
    message(paste("Processed batch", batch_num))
  }

  close(con)
}
```

### Q: How do I check memory usage?

**A:**

```r
# Check available memory
memory.size()  # Windows
gc()           # Garbage collection summary

# Monitor object size
object.size(sequences)  # In bytes
print(object.size(sequences), units = "MB")

# Typical sizes
# 1 million bp sequence = ~1 MB
# 1000 sequences × 30 kb each = ~30 MB
```

### Q: Can I use my GPU to speed up alignment?

**A:** DECIPHER doesn't support GPU acceleration, but you can:

- Use multiple CPU cores (as shown earlier)
- Use other tools like MUSCLE or MAFFT from command line
- Then load result back into R

---

## Common Errors

### Q: "Error: no valid letters in x"

**A:** You have non-DNA characters in your sequence. Fix:

```r
# View problematic character
invalid <- sequences[!grepl("^[ACGTN-]*$", as.character(sequences))]

# Clean the sequences
cleaned <- chartr("NRYSWKMBDHV", "N", sequences)  # Ambiguous → N
sequences <- cleaned  # Replace original

# Or be strict - only keep ACGT
strict_clean <- chartr("nryswkmbdhv", "A", sequences)
```

### Q: "Error in AlignSeqs: subscript out of bounds"

**A:** Likely empty sequences or mismatched data. Check:

```r
# Verify sequences
length(sequences)   # Should be > 1
width(sequences)    # Should all be > 0
any(width(sequences) == 0)  # Remove if TRUE

# Try with debug info
result <- tryCatch(
  AlignSeqs(sequences, verbose = TRUE),
  error = function(e) {
    print(e)
    NULL
  }
)
```

### Q: "Error: cannot coerce ... to DNAStringSet"

**A:** Wrong input format. Check:

```r
class(sequences)  # Should be "DNAStringSet"

# Fix by converting
if (class(sequences) == "list") {
  sequences <- DNAStringSet(sequences)
}

if (class(sequences) == "character") {
  sequences <- DNAStringSet(sequences)
}
```

### Q: R freezes during alignment

**A:** Usually means insufficient RAM. Solutions:

```r
# Check available memory
memory.limit()

# Restart R and clear memory
rm(list = ls())
gc()

# Try again with processors = 1
aligned <- AlignSeqs(sequences, processors = 1)

# Or reduce data size
sequences <- sequences[1:50]  # Test with subset first
```

---

## File Organization

### Q: Where should I put my FASTA files?

**A:** Your structure:

```
computational_biology/
├── classwork/
│   └── assets/
│       ├── wuhan.fasta
│       └── secuencias_fasta/
│           ├── mexico.fasta
│           ├── francia.fasta
│           └── [other sequences]
└── outputs/               # Save results here
    └── alignments/
```

### Q: How do I organize multiple projects?

**A:**

```
computational_biology/
├── projects/
│   ├── variant_analysis/
│   │   ├── data/
│   │   ├── scripts/
│   │   └── outputs/
│   ├── phylogenetics/
│   └── [other projects]
└── shared_resources/
    ├── globalSetup.R
    └── functions/
```

### Q: How do I ensure my analysis is reproducible?

**A:**

```r
# 1. Always save sessionInfo()
writeLines(
  capture.output(sessionInfo()),
  "outputs/session_info.txt"
)

# 2. Document your parameters
params <- list(
  alignment_method = "AlignSeqs",
  gap_opening = -25,
  gap_extension = -5,
  iterations = 2,
  date_run = Sys.Date()
)

# 3. Save everything
save(aligned, tree, params, file = "outputs/analysis.RData")

# 4. Later, reproduce with
load("outputs/analysis.RData")
```

---

## For Your Specific Project

### Q: I'm analyzing SARS-CoV-2 variants. Any special tips?

**A:** Yes! SARS-CoV-2 specific:

```r
# 1. The genome is ~29.9 kb
genome_size <- 29903

# 2. Known variant reference genomes
variants <- list(
  wuhan = "NC_045512.2",      # Original Wuhan
  alpha = "MW989814.1",        # B.1.1.7
  delta = "AY.100",            # B.1.617.2
  omicron = "OM429650.1"       # BA.1
)

# 3. Common features to analyze
orf1 <- c(1, 21563)            # ORF1a/b
spike <- c(21563, 25384)       # S protein (common mutations)
nucleocapsid <- c(28274, 29533) # N protein

# 4. Use these to extract regions
spike_seqs <- subseq(sequences, spike[1], spike[2])

# 5. Check common mutations
variants_markers <- list(
  alpha = "69-70del",         # S gene dropout
  delta = "L452R_P681R",
  omicron = "S373P_S375F"     # Part of BA.2 signature
)
```

### Q: What alignment parameters work best for variants?

**A:** For SARS-CoV-2 with small variations:

```r
# Conservative alignment (best for variants)
aligned <- AlignSeqs(
  sequences,
  gapOpening = -50,      # High penalty for gaps
  gapExtension = -10,    # Extension also costly
  iterations = 2,
  verbose = TRUE
)

# Compare variants only
dist_matrix <- dist.hamming(aligned)  # Hamming is fine for variants

# Check mutations per position
mutations_per_site <- colSums(aligned != aligned[[1]])
hist(mutations_per_site, main = "Mutations per position")
```

---

## Team Collaboration Tips

### Q: How do I share my analysis with the team?

**A:**

```bash
# Version control
git add classwork/
git commit -m "Add DECIPHER variant analysis"
git push

# Create readable output
# Save both .RData (for R) and CSV (for everyone)
write.csv(mutation_results, "outputs/mutations.csv", row.names = FALSE)
write.csv(tree_data, "outputs/tree_data.csv", row.names = FALSE)

# Document with Markdown
# Create analysis report as .Rmd or .md
```

### Q: How do I write analysis notes others can understand?

**A:** Use R Markdown and comments:

```r
# ===== ANALYSIS SECTION =====
# Purpose: Compare Mexico and Wuhan sequences
# Author: [Your name]
# Date: 2024-04-17

# 1. Load sequences
mexico <- readDNAStringSet("mexico.fasta")

# 2. Describe what you're doing
# We expect high similarity between variants within 2021
# This section identifies mutations

# 3. Show results clearly
print(paste("Found", nrow(mutations), "mutations"))
```

---

## Additional Resources

### Online Help

```r
# Within R
?AlignSeqs       # Function documentation
help(AlignSeqs)  # More detailed help

# Online
# Bioconductor: https://www.bioconductor.org/packages/DECIPHER
# Support: https://support.bioconductor.org/
# DECIPHER manual:
```

### Key Functions Quick List

```r
# Sequence I/O
readDNAStringSet()      # Load FASTA
writeXStringSet()       # Save FASTA

# Analysis
AlignSeqs()             # Multiple alignment
dist.hamming()          # Pairwise distances
nj()                    # Neighbor-joining tree

# Conversion
as.character()          # To text
as.matrix()             # To matrix format
DNAStringSet()          # Create DNAStringSet
```

---

**Still have questions?**

1. Check `DECIPHER_VS_CODE_GUIDE.md` for detailed examples
2. Check `DECIPHER_CHEATSHEET.md` for quick reference
3. Run `?AlignSeqs` in R for function help
4. Post to Bioconductor support: https://support.bioconductor.org/

Good luck with your analysis! 🧬
