Autor: Biniza Vázquez  
Realizado con apoyo de Claude Haiku 4.5  
Abril, 2026

---

# Guía Completa de DECIPHER en VS Code

Referencia técnica para el análisis de secuencias de ADN usando DECIPHER en entorno VS Code.

---

## Tabla de Contenidos

1. [Introducción](#introducción)
2. [Instalación y configuración](#instalación-y-configuración)
3. [Funciones principales](#funciones-principales)
4. [Trabajo con secuencias](#trabajo-con-secuencias)
5. [Alineamiento de secuencias](#alineamiento-de-secuencias)
6. [Análisis filogenético](#análisis-filogenético)
7. [Integración con VS Code](#integración-con-vs-code)
8. [Buenas prácticas](#buenas-prácticas)
9. [Ejemplos de análisis](#ejemplos-de-análisis)
10. [Solución de problemas](#solución-de-problemas)
11. [Técnicas avanzadas](#técnicas-avanzadas)
12. [Optimización de rendimiento](#optimización-de-rendimiento)

---

## Introducción

DECIPHER es un paquete de R para el análisis comparativo de secuencias de ADN que facilita:

- Alineamiento múltiple de secuencias usando algoritmos de programación dinámica
- Detección de mutaciones y variaciones entre muestras
- Cálculo de distancias evolutivas mediante múltiples modelos estadísticos
- Construcción de árboles filogenéticos para visualizar relaciones evolutivas
- Identificación de patrones y motivos en secuencias de ADN

En el contexto de tu proyecto, DECIPHER permite comparar variantes de SARS-CoV-2 y construir árboles filogenéticos que visualizan las relaciones evolutivas entre muestras.

DECIPHER se integra con otros paquetes de bioinformática:

- seqinr: lectura y escritura de archivos FASTA
- ape: métodos filogenéticos estadísticos
- ggtree: visualización especializada de árboles filogenéticos
- ggmsa: visualización de alineamientos múltiples

---

## Instalación y configuración

### Instalación de DECIPHER

DECIPHER se distribuye a través del repositorio Bioconductor, requiriendo BiocManager:

```r
# Instalación inicial de BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Instalación de DECIPHER desde Bioconductor
BiocManager::install("DECIPHER")
```

### Verificación de la instalación

```r
# Carga el paquete
library(DECIPHER)

# Consulta la versión instalada
packageVersion("DECIPHER")

# Información de citación
citation("DECIPHER")
```

### Configuración en VS Code

En el archivo `globalSetup.R`, implementa una función de carga robusta:

```r
# Función para carga segura de DECIPHER
load_decipher <- function() {
  if (!requireNamespace("DECIPHER", quietly = TRUE)) {
    message("DECIPHER no encontrado. Instalando desde Bioconductor...")
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("DECIPHER")
  }
  library(DECIPHER)
}

# Carga
load_decipher()
```

---

## Funciones principales

Las funciones más utilizadas en análisis de secuencias:

| Función | Propósito | Entrada | Salida |
|---------|-----------|---------|--------|
| `AlignSeqs()` | Alineamiento múltiple | DNAStringSet | Secuencias alineadas |
| `dist.hamming()` | Distancia de Hamming | Alineamiento | Matriz de distancias |
| `nj()` | Árbol neighbor-joining | Matriz distancias | Árbol filogenético |
| `readDNAStringSet()` | Lectura FASTA | Archivo .fasta | DNAStringSet |
| `writeXStringSet()` | Escritura de alineamiento | DNAStringSet | Archivo FASTA |
| `letterFrequency()` | Composición de bases | Secuencias | Tabla de frecuencias |
| `reverseComplement()` | Complemento inverso | DNAStringSet | Secuencia complementaria |
| `vmatchPattern()` | Búsqueda de motivos | Secuencias | Posiciones encontradas |

---

## Trabajo con secuencias

### Carga de secuencias

```r
library(DECIPHER)
library(Biostrings)

# Método recomendado usando readDNAStringSet
sequences <- readDNAStringSet("path/to/sequences.fasta")

# Conversión desde formato seqinr
seq_list <- read.fasta("sequences.fasta", as.string = TRUE)
sequences <- DNAStringSet(seq_list)

# Carga de archivos del proyecto
mexico_seq <- readDNAStringSet("classwork/assets/secuencias_fasta/mexico.fasta")
wuhan_seq <- readDNAStringSet("classwork/assets/wuhan.fasta")
francia_seq <- readDNAStringSet("classwork/assets/secuencias_fasta/francia.fasta")
```

### Inspección de secuencias

```r
# Información estructural
length(mexico_seq)              # Número de secuencias
width(mexico_seq)               # Longitud de cada secuencia
names(mexico_seq)               # Identificadores

# Cálculo de contenido GC
gc_content <- function(seq) {
  g_count <- letterFrequency(seq, "G")
  c_count <- letterFrequency(seq, "C")
  (g_count + c_count) / width(seq)
}

gc_mexico <- gc_content(mexico_seq[[1]])

# Composición de bases
letterFrequency(mexico_seq[[1]], letters = "ACGT")
```

### Manipulación de secuencias

```r
# Cálculo de complemento inverso
rev_comp_seq <- reverseComplement(mexico_seq)

# Extracción de región específica
spike_region <- subseq(mexico_seq[[1]], start = 21563, end = 25384)

# Limpieza de bases ambiguas
cleaned <- chartr("NRYSWKMBDHV", "N", mexico_seq)

# Identificación de secuencias únicas
unique_seqs <- unique(mexico_seq)
```

---

## Alineamiento de secuencias

El alineamiento múltiple organiza secuencias de forma que las posiciones homólogas se correspondan, permitiendo comparación directa. DECIPHER ejecuta esto mediante algoritmos de programación dinámica.

### Alineamiento básico

```r
library(DECIPHER)

# Combinación de secuencias para alineamiento
all_sequences <- c(mexico_seq, wuhan_seq, francia_seq)

# Alineamiento con reporte de progreso
aligned_sequences <- AlignSeqs(all_sequences,
                              verbose = TRUE,
                              iterations = 2)

# Visualización del resultado
print(aligned_sequences)
```

### Alineamiento personalizado

Para control más granular sobre los parámetros:

```r
aligned <- AlignSeqs(
  all_sequences,
  gapOpening = -25,      # Penalidad por apertura de gap
  gapExtension = -5,     # Penalidad por extensión de gap
  useQuality = FALSE,    # Sin ponderación por calidad
  verbose = TRUE,        # Reporte en consola
  iterations = 3,        # Iteraciones (mayor precisión)
  anchor = NA,           # Sin secuencias ancla
  processors = NULL      # Todos los procesadores
)
```

### Visualización del alineamiento

```r
library(ggmsa)

# Visualización en consola
print(aligned)

# Exportación a formato FASTA
writeXStringSet(aligned, "alignment.fasta")

# Visualización gráfica
writeXStringSet(aligned, "temp.fasta")
ggmsa("temp.fasta", start = 1, end = 500, color = "Clustal") +
  theme_minimal() +
  labs(title = "Alineamiento múltiple")
```

---

## Análisis filogenético

Los árboles filogenéticos representan gráficamente las relaciones evolutivas entre secuencias. La proximidad en el árbol refleja el grado de similitud entre muestras.

### Cálculo de distancias

```r
library(ape)

# Generación de matriz de distancias
aligned <- AlignSeqs(all_sequences)
dist_matrix <- as.matrix(dist.hamming(aligned))

# Métodos alternativos:
# dist.hamming(): distancia basada en diferencias directas
# dist.dna(): modelos evolutivos (JC69, TN93, etc.)
```

### Construcción del árbol

```r
# Método neighbor-joining (más común)
nj_tree <- nj(as.dist(dist_matrix))

# Método UPGMA (asume reloj molecular)
upgma_tree <- upgma(dist_matrix)

# Visualización básica
plot(nj_tree)
title("Árbol filogenético")
axisPhylo()
```

### Visualización avanzada con ggtree

```r
library(ggtree)

# Árbol rectangular con etiquetas
ggtree(nj_tree) +
  geom_tiplab() +
  theme_tree2() +
  labs(title = "Relaciones filogenéticas")

# Árbol con codificación de color
ggtree(nj_tree, layout = "rectangular") +
  geom_tiplab(aes(color = ifelse(grepl("mexico", label), "red", "blue"))) +
  geom_nodepoint(size = 3) +
  theme_tree2() +
  theme(legend.position = "right")
```

---

## Integración con VS Code

### Extensiones requeridas

Accede a la tienda de extensiones mediante Cmd + Shift + X:

- R (oficial): soporte de lenguaje y ejecución de código
- Git Lens: control de versiones integrado
- Markdown All in One: edición avanzada de documentación

### Configuración de VS Code

Accede a preferencias mediante Cmd + , y aplica:

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

### Ejecución de código

| Acción | Atajo (Mac) |
|--------|------------|
| Ejecutar línea actual | Cmd + Enter |
| Ejecutar bloque seleccionado | Cmd + Shift + Enter |
| Ejecutar script completo | Cmd + Shift + R |

---

## Buenas prácticas

### Gestión de memoria

Para archivos grandes, implementa procesamiento en lotes:

```r
process_large_fasta <- function(file_path, chunk_size = 100) {
  all_seqs <- readDNAStringSet(file_path)
  
  for (i in seq(1, length(all_seqs), by = chunk_size)) {
    chunk <- all_seqs[i:min(i + chunk_size - 1, length(all_seqs))]
    aligned <- AlignSeqs(chunk)
    # Procesamiento específico del lote
  }
}

# Liberación de memoria
rm(old_object)
gc()
```

### Manejo de errores

```r
# Ejecución robusta con manejo de excepciones
safe_align <- function(sequences) {
  tryCatch({
    AlignSeqs(sequences, verbose = TRUE)
  }, error = function(e) {
    message("Error durante alineamiento: ", conditionMessage(e))
    NULL
  })
}

# Validación previa de secuencias
validate_sequences <- function(seqs) {
  if (any(width(seqs) == 0)) {
    stop("Secuencias vacías detectadas")
  }
  if (!all(grepl("^[ACGTN-]*$", as.character(seqs)))) {
    warning("Bases ambiguas detectadas")
  }
  return(TRUE)
}
```

### Reproducibilidad

```r
# Fijación de semilla para resultados reproducibles
set.seed(42)

# Documentación de versiones
sessionInfo()

# Guardado con marca temporal
save_analysis <- function(object, name) {
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  filename <- paste0("outputs/", name, "_", timestamp, ".RData")
  save(object, file = filename)
  message("Guardado en: ", filename)
}
```

### Organización del código

```r
# Carga de configuración inicial
source("classwork/globalSetup.R")

# Lectura de datos
mexico_seq <- readDNAStringSet("classwork/assets/secuencias_fasta/mexico.fasta")
wuhan_seq <- readDNAStringSet("classwork/assets/wuhan.fasta")

# Análisis
aligned <- AlignSeqs(c(mexico_seq, wuhan_seq))

# Visualización
plot_results(aligned)
```

---

## Ejemplos de análisis

### Análisis integrado de variantes

```r
library(DECIPHER)
library(Biostrings)
library(ape)
library(ggplot2)

# Inicialización
setwd("/Users/binivazquez/UniWorkspace/computational_biology")
source("classwork/globalSetup.R")

# Carga de secuencias
mexico_seq <- readDNAStringSet("classwork/assets/secuencias_fasta/mexico.fasta")
wuhan_seq <- readDNAStringSet("classwork/assets/wuhan.fasta")
francia_seq <- readDNAStringSet("classwork/assets/secuencias_fasta/francia.fasta")

# Combinación
all_seqs <- c(mexico_seq, wuhan_seq, francia_seq)
names(all_seqs) <- c("Mexico", "Wuhan", "Francia")

# Alineamiento
aligned <- AlignSeqs(all_seqs, verbose = TRUE)

# Análisis de mutaciones
mutations <- data.frame()
for (i in 1:(length(aligned)-1)) {
  for (j in (i+1):length(aligned)) {
    seq1 <- aligned[[i]]
    seq2 <- aligned[[j]]
    diffs <- seq1 != seq2
    
    mutations <- rbind(mutations, data.frame(
      Pair = paste(names(aligned)[i], "vs", names(aligned)[j]),
      Mutations = sum(diffs),
      MutationRate = sum(diffs) / width(seq1)
    ))
  }
}

# Construcción del árbol
dist_mat <- as.matrix(dist.hamming(aligned))
tree <- nj(as.dist(dist_mat))

# Visualización
par(mfrow = c(1, 2))
plot(tree)
title("Árbol filogenético")
heatmap(dist_mat, main = "Matriz de distancias")
par(mfrow = c(1, 1))

# Almacenamiento de resultados
saveRDS(aligned, "outputs/alignment.rds")
write.csv(mutations, "outputs/mutations.csv")
```

### Procesamiento por lotes

```r
# Procesamiento automático de múltiples archivos
fasta_files <- list.files("data/", pattern = "\\.fasta$", full.names = TRUE)

for (file in fasta_files) {
  cat("Procesando:", file, "\n")
  
  seqs <- readDNAStringSet(file)
  aligned <- AlignSeqs(seqs)
  
  base_name <- sub(".fasta", "", basename(file))
  saveRDS(aligned, paste0("results/", base_name, "_aligned.rds"))
}

cat("Procesamiento completado\n")
```

---

## Solución de problemas

### AlignSeqs no disponible

```r
# Asegúrate de cargar DECIPHER
library(DECIPHER)

# Verifica instalación
if (!requireNamespace("DECIPHER", quietly = TRUE)) {
  BiocManager::install("DECIPHER")
}
```

### Objeto Biostrings inválido

```r
# Verifica el tipo de objeto
class(sequences)  # Debe mostrar "DNAStringSet"

# Incorrecto
bad <- "ACGTACGT"

# Correcto
good <- DNAStringSet(c("ACGTACGT", "TGCATGCA"))
```

### Alineamiento lento

```r
# Opción 1: Procesamiento paralelo
aligned <- AlignSeqs(sequences, processors = 4)

# Opción 2: Reducción de iteraciones
aligned <- AlignSeqs(sequences, iterations = 1)

# Opción 3: Procesamiento por lotes
batch_size <- 50
batches <- split(sequences, rep(1:ceiling(length(sequences)/batch_size), 
                               length.out = length(sequences)))
results <- lapply(batches, AlignSeqs)
```

### Caracteres inválidos en secuencia

```r
# Identificación de bases ambiguas
invalid <- sequences[!grepl("^[ACGTN-]*$", as.character(sequences))]

# Limpieza
cleaned <- chartr("NRYSWKMBDHV", "N", sequences)
```

### R deja de responder

```r
# Libera memoria
rm(list = ls())
gc()

# Reinicia con parámetros más conservadores
aligned <- AlignSeqs(sequences, processors = 1)

# Prueba con subset más pequeño
test_seqs <- sequences[1:10]
```

---

## Técnicas avanzadas

### Búsqueda de patrones

```r
# Localización de codones START
start_codons <- vmatchPattern("ATG", sequences)

# Conteo de ocurrencias
counts <- vcountPattern("ATG", sequences)

# Frecuencia de trinucleótidos
freq <- oligonucleotideFrequency(sequences, width = 3)
```

### Secuencia consenso

```r
# Generación de secuencia consenso
get_consensus <- function(aligned_sequences) {
  aligned_matrix <- as.character(aligned_sequences)
  consensus <- character(width(aligned_sequences[1]))
  
  for (i in 1:width(aligned_sequences[1])) {
    bases <- substr(aligned_matrix, i, i)
    base_table <- table(bases)
    consensus[i] <- names(base_table)[which.max(base_table)]
  }
  
  DNAString(paste(consensus, collapse = ""))
}
```

### Bootstrap y soporte nodal

```r
library(ape)

# Función de remuestreo
boot_func <- function(x) {
  d <- dist.hamming(x)
  nj(as.dist(d))
}

# Generación de 100 réplicas bootstrap
boot_trees <- boot.phylo(tree, aligned, boot_func, rooted = TRUE)

# Visualización con valores de soporte
plot(tree)
nodelabels(boot_trees, cex = 0.7)
```

---

## Optimización de rendimiento

### Procesamiento paralelo

```r
library(parallel)

# Detección de núcleos disponibles
n_cores <- detectCores() - 1

# Alineamiento paralelo
aligned <- AlignSeqs(sequences, 
                    processors = n_cores, 
                    verbose = TRUE)
```

### Almacenamiento en caché

```r
# Guardado del alineamiento
saveRDS(aligned, "cache/alignment.rds")

# Recuperación posterior
aligned <- readRDS("cache/alignment.rds")
```

### Monitoreo de memoria

```r
# Consulta de tamaño de objeto
object.size(sequences)
print(object.size(sequences), units = "MB")

# Compactación de memoria
gc()
```

---

## Referencia rápida

Comandos esenciales para análisis típicos:

```r
# Instalación y carga
if (!requireNamespace("BiocManager")) install.packages("BiocManager")
BiocManager::install("DECIPHER")
library(DECIPHER)

# Flujo de análisis
seqs <- readDNAStringSet("sequences.fasta")
aligned <- AlignSeqs(seqs, verbose = TRUE)
library(ape)
tree <- nj(as.dist(dist.hamming(aligned)))
plot(tree)
saveRDS(aligned, "resultado.rds")
```

---

## Recursos de referencia

**Documentación oficial:**
- [DECIPHER Manual](https://www.bioconductor.org/packages/DECIPHER)
- [Biostrings](https://bioconductor.org/packages/Biostrings)
- [ape Package](https://cran.r-project.org/web/packages/ape/)

**Soporte técnico:**
- [Bioconductor Support](https://support.bioconductor.org/)
- Stack Overflow con tags `[r]` y `[bioconductor]`

---

**Última actualización:** Abril 2026
