# VS Code Setup Guide for DECIPHER Work

## Prerequisites

- VS Code installed
- R installed (`/usr/local/bin/R` on your Mac)
- Command line tools/XCode

## Step 1: Install VS Code Extensions

### Required Extensions

1. **R (Positron)** or **R Client**
   - Search: "R" in Extensions
   - Install the official R extension

2. **Git Lens** (optional but recommended)
   - For tracking your code changes

3. **Markdown All in One**
   - For better markdown editing

### Installation Steps in VS Code:

```
Cmd + Shift + X  →  Search "R"  →  Install
```

## Step 2: Configure VS Code Settings

Open `Code > Preferences > Settings` or press `Cmd + ,`

Add these settings (search for "R" in settings):

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

## Step 3: Verify R Installation in Terminal

```bash
which R              # Should show: /usr/local/bin/R
R --version         # Shows R version
Rscript --version   # Shows Rscript version
```

## Step 4: Initial R Setup

When you first use R in VS Code:

1. Open Command Palette: `Cmd + Shift + P`
2. Type: "R: Create R Terminal"
3. This starts an interactive R session

## Step 5: Install Required Bioconductor Packages

In the R terminal:

```r
# Install BiocManager first
install.packages("BiocManager")

# Install DECIPHER from Bioconductor
BiocManager::install("DECIPHER")

# Install other bioinformatics packages
BiocManager::install(c("Biostrings", "ape", "ggtree", "ggmsa"))

# Install CRAN packages
install.packages(c("seqinr", "ggplot2", "dplyr", "tidyr"))

# Verify installations
library(DECIPHER)
library(Biostrings)
packageVersion("DECIPHER")
```

## Step 6: Project Structure Setup

In your workspace folder:

```bash
mkdir -p classwork/assets/secuencias_fasta classwork/assets/secuencias_gen
mkdir -p assignments
mkdir -p outputs
mkdir -p cache
mkdir -p functions
touch README.md
```

## Step 7: Create .gitignore

Create `.gitignore` in your project root:

```
# Large data files
*.fasta
*.fa
*.gb
*.gbk

# R output
*.png
*.pdf
*.html
.Rhistory
.RData

# Cache
cache/
outputs/temp/

# IDE
.Rproj.user/
.vscode/

# OS
.DS_Store
Thumbs.db
```

## Step 8: Create globalSetup.R

This should already exist, but ensure it has:

```r
# globalSetup.R - Shared configuration for all analyses

# Load packages safely
load_if_available <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' not installed. Install with: BiocManager::install('%s')", pkg, pkg))
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# Load required packages
load_if_available("seqinr")
load_if_available("ape")
load_if_available("ggplot2")
load_if_available("dplyr")
load_if_available("Biostrings")

# Load optional packages
tryCatch(load_if_available("DECIPHER"), error = function(e) warning(e))
tryCatch(load_if_available("ggtree"), error = function(e) warning(e))
tryCatch(load_if_available("ggmsa"), error = function(e) warning(e))

# Set working directory
if (!exists("assets_dir")) {
  assets_dir <- file.path(dirname(getwd()), "assets")
}

# Print session info
cat("\n=== Computational Biology Setup ===\n")
cat("Working directory:", getwd(), "\n")
cat("R version:", R.version$version.string, "\n")
```

## Step 9: Create Your First Analysis File

Create `classwork/test_decipher.R`:

```r
# Test DECIPHER installation and basic workflow
setwd("/Users/binivazquez/UniWorkspace/computational_biology")
source("classwork/globalSetup.R")

cat("\n=== Testing DECIPHER ===\n")

# Load test sequences
mexico <- readDNAStringSet("classwork/assets/secuencias_fasta/mexico.fasta")
wuhan <- readDNAStringSet("classwork/assets/wuhan.fasta")

cat("Sequences loaded:\n")
cat("Mexico:", length(mexico), "sequences\n")
cat("Wuhan:", length(wuhan), "sequences\n")

# Combine and align
all_seqs <- c(mexico, wuhan)
names(all_seqs) <- c("Mexico", "Wuhan")

cat("\nAligning sequences...\n")
aligned <- AlignSeqs(all_seqs, verbose = TRUE)

cat("\nAlignment complete!\n")
cat("Alignment width:", width(aligned[1]), "bp\n")

# Save result
saveRDS(aligned, "outputs/test_alignment.rds")
cat("Alignment saved to: outputs/test_alignment.rds\n")
```

## Step 10: Run Your First Script

In VS Code:

1. Open `test_decipher.R`
2. Press `Cmd + Shift + Enter` to run the entire script
3. Or use `Cmd + Enter` to run line by line

You should see output in the R terminal showing:

- Sequences loaded
- Alignment progress
- Results saved

## VS Code Keyboard Shortcuts for R

| Action                    | Shortcut                                 |
| ------------------------- | ---------------------------------------- |
| Create R Terminal         | Cmd + Shift + P → "R: Create R Terminal" |
| Run current line          | Cmd + Enter                              |
| Run selection             | Cmd + Enter (with text selected)         |
| Run entire script         | Cmd + Shift + Enter                      |
| Run from current to end   | Cmd + Shift + Alt + E                    |
| Run from start to current | Cmd + Shift + Alt + B                    |
| Insert code section       | Cmd + Shift + I                          |

## Useful VS Code Commands for R

Open Command Palette (`Cmd + Shift + P`) and type:

- `R: Create R Terminal` - Start R session
- `R: Create R Script` - New R file
- `R: Run Command` - Run custom R command
- `R: Show Help` - Access R documentation
- `R: Open Html Preview` - View HTML output

## Common Issues & Fixes

### Issue: "R command not found"

```bash
# Check R installation path
which R

# If not found, install R from: https://cran.r-project.org/bin/macosx/
# After installation, restart VS Code
```

### Issue: Package not found when running script

```r
# In R terminal:
.libPaths()  # Check library paths
library()    # List all installed packages
```

### Issue: Permission denied on Mac

```bash
# Fix R permissions
sudo chown -R $(whoami) /usr/local/bin/R
```

### Issue: DECIPHER alignment is very slow

```r
# Use parallel processing (adjust to your system)
aligned <- AlignSeqs(sequences, processors = 2, iterations = 1)
```

## Tips for Efficient Workflow

### Tip 1: Use R Markdown for Analysis

Create `classwork/analysis.Rmd`:

````
---
title: "SARS-CoV-2 Analysis"
author: "Your Team"
output: html_document
---

```{r setup}
source("globalSetup.R")
```

## Load Data

```{r load}
sequences <- readDNAStringSet("path/to/file.fasta")
```

## Align

```{r align}
aligned <- AlignSeqs(sequences)
```
````

Then render with: `Cmd + Shift + K`

### Tip 2: Keep Helper Functions in Separate File

Create `functions/alignment_helpers.R`:

```r
# Quick alignment and tree building
quick_analysis <- function(fasta_file) {
  seqs <- readDNAStringSet(fasta_file)
  aligned <- AlignSeqs(seqs, verbose = TRUE)
  dist_mat <- as.matrix(dist.hamming(aligned))
  tree <- nj(as.dist(dist_mat))
  list(aligned = aligned, tree = tree, dist = dist_mat)
}
```

Then use in your analysis:

```r
source("functions/alignment_helpers.R")
results <- quick_analysis("seqs.fasta")
```

### Tip 3: Save Frequently

```r
# Save results immediately after computation
saveRDS(aligned, "outputs/aligned_sequences.rds")

# Quick save function
save_result <- function(object, name) {
  filename <- paste0("outputs/", name, "_", Sys.Date(), ".rds")
  saveRDS(object, filename)
  message("Saved to: ", filename)
}

# Usage
save_result(aligned, "mexico_wuhan_alignment")
```

### Tip 4: Use Version Control

```bash
cd /Users/binivazquez/UniWorkspace/computational_biology

git init
git add .
git commit -m "Initial commit: DECIPHER analysis setup"

# Regular commits
git add classwork/
git commit -m "Added DECIPHER alignment analysis"
git push origin main
```

## Monthly Maintenance

1. **Update packages** (once a month):

```r
update.packages(ask = FALSE)
BiocManager::install()  # Update Bioconductor packages
```

2. **Clean cache**:

```bash
rm -rf cache/*
```

3. **Archive old results**:

```bash
mkdir -p archive/$(date +%Y-%m)
mv outputs/* archive/$(date +%Y-%m)/
```

## Quick Reference for Your Setup

**Your R is at:** `/usr/local/bin/R`
**Your project:** `/Users/binivazquez/UniWorkspace/computational_biology`
**DECIPHER install command:** `BiocManager::install("DECIPHER")`

## Next Steps

1. ✅ Install DECIPHER and dependencies
2. ✅ Run `test_decipher.R` to verify setup
3. Read `DECIPHER_VS_CODE_GUIDE.md` for detailed examples
4. Check `DECIPHER_CHEATSHEET.md` while coding
5. Start your analysis!

---

**Everything is now configured. You're ready to work!** 🧬

Questions? Check:

- Terminal output for specific error messages
- Package documentation: `?AlignSeqs`
- Bioconductor support: https://support.bioconductor.org/
