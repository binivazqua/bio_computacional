setwd("/Users/binivazquez/UniWorkspace/computational_biology")
source("classwork/globalSetup.R")


seqs_gb <- read.GenBank("/Users/binivazquez/UniWorkspace/computational_biology/classwork/assets/secuencias_gen/all_arn_sequences.gb")
# Esto lee SOLO el archivo.gb
# Si tiene múltiples registros dentro, los lee todos
# Si tiene una sola secuencia, lee esa

seqs_gb_dna <- DNAStringSet(seqs_gb$sequence)
# Extrae las secuencias y las pone en un DNAStringSet
# Si archivo.gb tiene 3 secuencias, aquí crea un conjunto de 3

aligned_gb <- AlignSeqs(seqs_gb_dna)
# Alinea lo que haya (1, 3, o N secuencias)
