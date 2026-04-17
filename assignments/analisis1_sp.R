source("classwork/globalSetup.R")

# Load FASTA files (globalSetup.R sets assets_dir)
mexico <- read.fasta(file.path(assets_dir, "mexico.fasta"))
wuhan <- read.fasta(file.path(assets_dir, "wuhan.fasta"))
francia <- read.fasta(file.path(assets_dir, "francia.fasta"))

print("Longitud de la secuencia del virus de Mexico: ")
length(mexico[[1]])

print("Longitud de la secuencia del virus de Wuhan: ")
length(wuhan[[1]])

print("Longitud de la secuencia del virus de Francia: ")
length(francia[[1]])

compara <- function(mex, wuh, fra) {
    par(mfrow = c(3, 1))
    barplot(table(mex), col = 1:4, name = "mex")
    barplot(table(wuh), col = 1:4, name = "wuh")
    barplot(table(fra), col = 1:4, name = "fra")
}
compara(mexico[[1]], wuhan[[1]], francia[[1]])
