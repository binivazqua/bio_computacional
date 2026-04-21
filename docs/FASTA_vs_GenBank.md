# FASTA vs GenBank: Comparativa y Decisión

## ¿Se puede usar FASTA?

Con la nueva doc de DECIPHER, sí, más directo y rápido.

## Comparativa

| Aspecto                | FASTA                                      | GenBank (.gen)                            |
| ---------------------- | ------------------------------------------ | ----------------------------------------- |
| **Lectura**            | `readDNAStringSet()`                       | `Seqs2DB()` requiere base de datos SQLite |
| **Velocidad**          | Más rápida (directo en memoria)            | Más lenta (requiere BD)                   |
| **Metadata**           | Solo descripción simple                    | LOCUS, DEFINITION, ORIGIN, AUTHORS, etc.  |
| **Alineamiento**       | `AlignSeqs()` directo                      | `SearchDB()` luego `AlignSeqs()`          |
| **Árbol filogenético** | `Treeline()` o `nj()`                      | `Treeline()` o `nj()`                     |
| **Complejidad**        | Bajo: 5 líneas de código                   | Medio: 10+ líneas                         |
| **Ideal para**         | Análisis rápidos, educación, visualización | Investigación seria, traceo de genes      |

## Recomendación por uso

### Opción 1: USA FASTA (recomendado para usarlo en clase)

- Más simple y directo
- Menos requisitos de dependencias
- Árboles igualmente buenos
- Mejor para visualización con headers limpios
- **Ideal para: Demos, trabajos en clase, presentaciones**

```r
library(Biostrings)
library(DECIPHER)
library(ape)

fasta_file <- "classwork/assets/secuencias_fasta/all_arn_sequences_CLEAN.fasta"

seqs <- readDNAStringSet(fasta_file)
aligned <- AlignSeqs(seqs, verbose = TRUE)
d <- DistanceMatrix(aligned, correction="K80", verbose=FALSE)
tree <- DECIPHER::Treeline(myXStringSet=aligned, method="NJ", cutoff=0.05, showPlot=TRUE)

plot(tree)
```

### Opción 2: USA GenBank (.gen)

- Preserva metadata completa
- Mejor para análisis downstream
- Rastreable a registros originales
- **Ideal para: Publicaciones, estudios serios, replicabilidad**

```r
library(DECIPHER)
library(ape)

gb_file <- "classwork/assets/secuencias_gb/all_sequences.gen"
dbConn <- dbConnect(RSQLite::SQLite(), ":memory:")

Seqs2DB(gb_file, type = "GenBank", dbConn = dbConn, identifier = "coronavirus")
dna <- SearchDB(dbConn)
d <- DistanceMatrix(dna, correction="K80", verbose=FALSE)
tree <- DECIPHER::Treeline(myXStringSet=dna, method="NJ", cutoff=0.05, showPlot=TRUE)

plot(tree)
dbDisconnect(dbConn)
```

## ¿Por qué usa GenBank en la documentación oficial?

**Razón técnica:** Los registros GenBank de NCBI contienen anotaciones que sirven para análisis posteriores (identificación de mutaciones, regiones codificantes, etc.).

**Razón pedagógica:** Enseña a los estudiantes a trabajar con formatos estándar de bioinformática.

## Decisión para tu proyecto

**Recomendación: USA FASTA + headers limpios**

Porque:

1. **Más simple** para estudiantes
2. **Árboles igual de buenos** (la secuencia es la misma)
3. **Headers legibles** con el script de limpieza
4. **Menos dependencias** (no necesitas RSQLite/DBI)
5. **Más rápido** de ejecutar

**GenBank si:**

- Necesitas anotar mutaciones específicas
- Quieres análisis downstream (CDS, ORF, etc.)
- Planejas publicar el análisis

## Conclusión

**No es obligación usar .gen.** Podemos usar FASTA sin problema. La diferencia es:

- **FASTA = rápido y simple**
- **GenBank = completo y trazable**
