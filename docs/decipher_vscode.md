Autoría: Biniza Vázquez y OpenAI Codex  
Abril de 2026

# DECIPHER en VS Code

Esta guía describe el flujo de trabajo recomendado para usar DECIPHER en VS Code dentro de este repositorio. Toma como documento rector `docs/guia_que_si_funciona.md` y corrige una fuente común de errores: leer archivos GenBank como si fueran FASTA.

## Alcance

Esta guía no reemplaza la configuración básica de VS Code. Su propósito es dejar claro qué archivo abrir, qué lector usar y cómo interpretar el error más frecuente observado en la terminal de R.

## Rutas válidas en este proyecto

```text
classwork/assets/secuencias_fasta/all_arn_sequences.fasta
classwork/assets/secuencias_gb/all_sequences.gen
```

## Regla de decisión antes de correr código

- Si el archivo empieza con `>`, se lee como FASTA con `readDNAStringSet()`.
- Si el archivo empieza con `LOCUS`, se trata como GenBank y se carga con `Seqs2DB(..., type = "GenBank")`.

## Flujo 1: trabajar con FASTA en VS Code

Úsalo para alineamiento directo y árboles sencillos.

```r
library(Biostrings)
library(DECIPHER)
library(ape)

fasta_file <- "classwork/assets/secuencias_fasta/all_arn_sequences.fasta"

seqs <- readDNAStringSet(fasta_file)
aligned <- AlignSeqs(seqs, verbose = TRUE)
dist_mat <- dist.hamming(aligned)
tree <- nj(as.dist(dist_mat))

plot(tree)
```

## Flujo 2: trabajar con GenBank en VS Code

Úsalo cuando conservas la estructura GenBank y quieres importar desde una base temporal con DECIPHER.

```r
library(DBI)
library(RSQLite)
library(DECIPHER)
library(ape)

gb_file <- "classwork/assets/secuencias_gb/all_sequences.gen"
dbConn <- dbConnect(RSQLite::SQLite(), ":memory:")

Seqs2DB(
  gb_file,
  type = "GenBank",
  dbConn = dbConn,
  identifier = "coronavirus"
)

dna <- SearchDB(dbConn)
d <- DistanceMatrix(dna, correction = "Jukes-Cantor")
tree <- TreeLine(myDistMatrix = d)

plot(as.phylo(tree))
dbDisconnect(dbConn)
```

## Error observado en la terminal de R

```text
Error in .Call2("fasta_index", filexp_list, nrec, skip, seek.first.rec, :
  reading FASTA file /Users/binivazquez/UniWorkspace/computational_biology/classwork/assets/secuencias_gb/all_sequences.gen: ">" expected at beginning of line 1
```

## Qué significa ese error

- VS Code y R estaban intentando usar un lector FASTA.
- El archivo `all_sequences.gen` no es FASTA.
- La primera línea del archivo es `LOCUS ...`, no `>...`.
- La corrección no es editar la terminal, sino usar el lector correcto.

## Comprobaciones rápidas en la terminal de VS Code

```r
file.exists("classwork/assets/secuencias_gb/all_sequences.gen")
file.exists("classwork/assets/secuencias_fasta/all_arn_sequences.fasta")
```

```bash
sed -n '1p' classwork/assets/secuencias_gb/all_sequences.gen
sed -n '1p' classwork/assets/secuencias_fasta/all_arn_sequences.fasta
```

Interpretación:

- Si ves `LOCUS`, sigue el flujo GenBank.
- Si ves `>`, sigue el flujo FASTA.

## Qué conviene usar en clase

- Para mostrar un alineamiento rápido: FASTA.
- Para respetar registros GenBank combinados: `Seqs2DB(..., type = "GenBank")`.
- Para depurar errores: inspeccionar la primera línea del archivo antes de cambiar funciones.

## Recomendación final

En este repositorio no conviene asumir que `.gen` significa automáticamente FASTA ni que `.gb` es el único contenedor válido de GenBank. Lo importante es verificar el contenido y luego elegir la función adecuada.
