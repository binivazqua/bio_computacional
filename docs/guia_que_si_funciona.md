# Guía que sí funciona para DECIPHER

Documento elaborado por Biniza Vázquez y OpenAI Codex.

Esta guía resume el flujo que sí funcionó en este proyecto para construir un árbol filogenético con DECIPHER a partir de archivos GenBank y, cuando conviene, con FASTA. La regla principal es sencilla: el contenido del archivo determina cómo debe leerse en R. La extensión por sí sola no basta.

## Objetivo

Trabajar con las secuencias de SARS-CoV-2 del repositorio sin mezclar dos flujos que DECIPHER trata de manera distinta:

- Flujo GenBank: para archivos como `classwork/assets/secuencias_gb/all_sequences.gen`
- Flujo FASTA: para archivos como `classwork/assets/secuencias_fasta/all_arn_sequences.fasta`

## Estructura del proyecto que se toma como válida

```text
classwork/assets/
├── secuencias_fasta/
│   └── all_arn_sequences.fasta
├── secuencias_gb/
│   ├── alphaus.gb
│   ├── betaghana.gb
│   ├── deltagermany.gb
│   ├── gammaitaly.gb
│   ├── wuhanog.gb
│   └── all_sequences.gen
└── wuhan.fasta
```

## Regla principal de lectura

- Si el archivo contiene encabezados FASTA que comienzan con `>`, úsalo con `readDNAStringSet()`.
- Si el archivo contiene registros GenBank con bloques `LOCUS`, `DEFINITION`, `ORIGIN` y cierre `//`, úsalo con `Seqs2DB(..., type = "GenBank")`.
- No intentes leer un archivo GenBank con una función pensada para FASTA.

## Flujo recomendado para GenBank

Este es el flujo más estable para el archivo combinado en `classwork/assets/secuencias_gb/all_sequences.gen`.

### 1. Unir los archivos GenBank

Desde la terminal:

```bash
cat classwork/assets/secuencias_gb/*.gb > classwork/assets/secuencias_gb/all_sequences.gen
```

### 2. Validar la estructura mínima

Cada registro GenBank debe terminar con `//`.

```bash
grep -c '^//$' classwork/assets/secuencias_gb/all_sequences.gen
```

El resultado debe coincidir con el número de archivos `.gb` incluidos.

### 3. Cargar y procesar en R

```r
library(DBI)
library(RSQLite)
library(DECIPHER)
library(ape)

gb_file <- "classwork/assets/secuencias_gb/all_sequences.gen"

dbConn <- dbConnect(RSQLite::SQLite(), ":memory:")

Seqs2DB(
  "classwork/assets/secuencias_gb/all_sequences.gen",
  type = "GenBank",
  dbConn,
  "coronavirus"
)

dna <- SearchDB(dbConn)
length(dna)

d <- DistanceMatrix(dna, correction = "Jukes-Cantor")
tree <- TreeLine(myDistMatrix = d)

plot(as.phylo(tree))

dbDisconnect(dbConn)
```

### 4. Verificaciones útiles en R

```r
file.exists(gb_file)
length(dna)
class(tree)
```

## Flujo alternativo para FASTA

Si tu archivo es FASTA, entonces el flujo correcto es otro.

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

## Error real que sí apareció y cómo interpretarlo

En la terminal de R apareció este error:

```text
Error in .Call2("fasta_index", filexp_list, nrec, skip, seek.first.rec, :
  reading FASTA file /Users/binivazquez/UniWorkspace/computational_biology/classwork/assets/secuencias_gb/all_sequences.gen: ">" expected at beginning of line 1
```

La interpretación correcta es la siguiente:

- `readDNAStringSet()` o una ruta de lectura FASTA intentó abrir `all_sequences.gen`.
- El archivo no era FASTA.
- La primera línea empezaba con `LOCUS` y no con `>`.
- DECIPHER falló porque estaba usando el lector equivocado.

## Qué no hacer

No mezcles estos dos casos:

```r
readDNAStringSet("classwork/assets/secuencias_gb/all_sequences.gen")
```

Eso falla si el archivo contiene GenBank.

Tampoco uses esto sobre un FASTA:

```r
Seqs2DB("classwork/assets/secuencias_fasta/all_arn_sequences.fasta",
        type = "GenBank",
        dbConn = dbConn,
        identifier = "coronavirus")
```

Eso falla porque el archivo no contiene registros GenBank.

## Diagnóstico rápido

### Caso 1. El archivo es GenBank

Primera línea esperada:

```text
LOCUS ...
```

Lector correcto:

```r
Seqs2DB(file, type = "GenBank", dbConn = dbConn, identifier = "coronavirus")
```

### Caso 2. El archivo es FASTA

Primera línea esperada:

```text
>identificador
```

Lector correcto:

```r
readDNAStringSet(file)
```

## Recomendación práctica para este repositorio

- Usa `classwork/assets/secuencias_gb/all_sequences.gen` cuando quieras conservar el flujo GenBank y cargar secuencias con `Seqs2DB`.
- Usa `classwork/assets/secuencias_fasta/all_arn_sequences.fasta` cuando quieras trabajar directamente con `readDNAStringSet()` y `AlignSeqs()`.
- Si un archivo produce errores de lectura, revisa la primera línea antes de cambiar el código.

## Resumen operativo

1. Decide si tu archivo es GenBank o FASTA.
2. Usa el lector correcto para ese formato.
3. No tomes la extensión `.gen` como garantía de contenido GenBank.
4. En este proyecto, el archivo combinado de GenBank se procesa con `Seqs2DB(..., type = "GenBank")`.
5. En este proyecto, el archivo combinado FASTA se procesa con `readDNAStringSet()`.
