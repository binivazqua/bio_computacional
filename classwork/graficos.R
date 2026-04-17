library(seqinr)

setwd("classwork/assets/")

mexico <- read.fasta("classwork/assets/mexico.fasta")
wuhan <- read.fasta("classwork/assets/wuhan.fasta")
francia <- read.fasta("classwork/assets/francia.fasta")

library(adegenet)
library(ape)
library(ggtree) ### https://bioconductor.org/packages/release/bioc/html/ggtree.html
library(DECIPHER)
library(viridis)
library(ggmsa) ### http://www.bioconductor.org/packages/release/bioc/html/ggmsa.html
library(ggplot2)
library(dplyr)

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
