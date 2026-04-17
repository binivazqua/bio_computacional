setwd("/Users/binivazquez/UniWorkspace/computational_biology")
source("classwork/globalSetup.R")

# Load FASTA files (globalSetup.R sets assets_dir)
mexico <- read.fasta(file.path(sars_assets_dir, "mexico.fasta"))
wuhan <- read.fasta(file.path(sars_assets_dir, "wuhanog.fasta"))
francia <- read.fasta(file.path(sars_assets_dir, "francia.fasta"))

print("Longitud de la secuencia del virus de Mexico: ")
length(mexico[[1]])

print("Longitud de la secuencia del virus de Wuhan: ")
length(wuhan[[1]])

print("Longitud de la secuencia del virus de Francia: ")
length(francia[[1]])

compara <- function(mex, wuh, fra) {
    par(mfrow = c(3, 1))
    barplot(table(mex), col = 1:4, main = "Mexico")
    barplot(table(wuh), col = 1:4, main = "Wuhan")
    barplot(table(fra), col = 1:4, main = "Francia")
}
compara(mexico[[1]], wuhan[[1]], francia[[1]])
