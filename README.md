# Análisis de secuencias con DECIPHER en VS Code

Documento de referencia del repositorio para trabajar con secuencias de SARS-CoV-2 en R y DECIPHER.

## Autoría

Contenido elaborado por Biniza Vázquez y OpenAI Codex.

## Documento principal

La referencia operativa de este proyecto es [docs/guia_que_si_funciona.md](docs/guia_que_si_funciona.md). Si algún otro archivo de documentación contradice esa guía, debe prevalecer la guía principal.

## Documentación disponible

- [docs/guia_que_si_funciona.md](docs/guia_que_si_funciona.md): flujo validado del proyecto
- [docs/decipher_vscode.md](docs/decipher_vscode.md): uso de DECIPHER dentro de VS Code
- [docs/setup_vscode.md](docs/setup_vscode.md): configuración mínima de VS Code para R

## Estructura relevante del repositorio

```text
classwork/assets/
├── secuencias_fasta/
│   └── all_arn_sequences.fasta
├── secuencias_gb/
│   └── all_sequences.gen
└── wuhan.fasta
```

## Regla principal sobre formatos

- `readDNAStringSet()` es para archivos FASTA.
- `Seqs2DB(..., type = "GenBank")` es para archivos GenBank.
- La extensión del archivo no reemplaza la validación del contenido.

## Flujos soportados

### Flujo FASTA

```r
library(Biostrings)
library(DECIPHER)
library(ape)

seqs <- readDNAStringSet("classwork/assets/secuencias_fasta/all_arn_sequences.fasta")
aligned <- AlignSeqs(seqs, verbose = TRUE)
dist_mat <- dist.hamming(aligned)
tree <- nj(as.dist(dist_mat))
plot(tree)
```

### Flujo GenBank

```r
library(DBI)
library(RSQLite)
library(DECIPHER)
library(ape)

dbConn <- dbConnect(RSQLite::SQLite(), ":memory:")

Seqs2DB(
  "classwork/assets/secuencias_gb/all_sequences.gen",
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

## Error que motivó esta documentación

Si aparece este mensaje:

```text
reading FASTA file .../all_sequences.gen: ">" expected at beginning of line 1
```

la causa es que se intentó leer un archivo GenBank como si fuera FASTA.

## Recomendación práctica

Antes de ejecutar un análisis, inspecciona la primera línea del archivo:

- `LOCUS` indica flujo GenBank
- `>` indica flujo FASTA

## Entrega Evidencia 1

En `deliverables/` se preparó una versión compatible con los archivos solicitados para la entrega:

- `FJ26_BT1013_201_BinizaVazquez_Evidencia1.Rmd`
- `FJ26_BT1013_201_BinizaVazquez_Evidencia1.html`
- `FJ26_BT1013_201_BinizaVazquez_Ev1-Inputs.zip`

La copia entregable del `.Rmd` está pensada para abrirse en RStudio o un entorno equivalente sin depender de VS Code, siempre que el archivo `.Rmd` y la carpeta de inputs descomprimida se mantengan en el mismo directorio.

Repositorio de respaldo:

- <https://github.com/binivazqua/bio_computacional>
