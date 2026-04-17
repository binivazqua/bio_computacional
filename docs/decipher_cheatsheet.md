# DECIPHER en VS Code - Guía Rápida para Estudiantes

**¡Dale! Aquí va todo lo que necesitas para no andar perdido.**

## Instalación (Solo una vez)

```r
# Instala desde Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DECIPHER")

library(DECIPHER)  # Carga la librería
```

---

## Cargar Secuencias

```r
library(Biostrings)

# De un archivo FASTA
sequences <- readDNAStringSet("path/to/file.fasta")

# Varios archivos a la vez
seq1 <- readDNAStringSet("file1.fasta")
seq2 <- readDNAStringSet("file2.fasta")
all_seqs <- c(seq1, seq2)

# De GenBank
sequences <- readBStringSet("file.gb")
```

---

## Alineamiento (Lo más importante)

```r
# Alineamiento rápido
aligned <- AlignSeqs(sequences, verbose = TRUE)

# Alineamiento bien hecho (tarda más pero vale la pena)
aligned <- AlignSeqs(
  sequences,
  iterations = 2,
  gapOpening = -25,
  gapExtension = -5
)

# Mira el resultado
print(aligned)
writeXStringSet(aligned, "output.fasta")
```

---

## Análisis de Secuencias

```r
# Propiedades básicas
width(sequences)                    # Longitudes
names(sequences)                    # Nombres
length(sequences)                   # Cuántas hay

# Composición de bases
letterFrequency(sequences, "ACGT")

# Contenido GC (importante en biología)
gc <- letterFrequency(sequences, "GC")
gc_percent <- rowSums(gc) / width(sequences) * 100

# Secuencia invertida complementaria
revcomp <- reverseComplement(sequences)

# Distancia entre secuencias
dist <- dist.hamming(aligned)
```

---

## Construir Árbol Filogenético

```r
library(ape)

# De las secuencias alineadas
dist_mat <- as.matrix(dist.hamming(aligned))

# Construye el árbol (elige uno)
nj_tree <- nj(as.dist(dist_mat))      # Neighbor-Joining (lo más común)
upgma_tree <- upgma(dist_mat)         # UPGMA (asume reloj molecular)

# Dibuja el árbol
plot(nj_tree)
title("Árbol Filogenético")
axisPhylo()  # Agrega barra de escala
```

---

## Visualización (Gráficos Bonitos)

```r
library(ggplot2)
library(ggtree)
library(ggmsa)

# Árbol chulo con ggtree
ggtree(tree) +
  geom_tiplab() +
  theme_tree2()

# Ver alineamiento (primero guarda como FASTA)
writeXStringSet(aligned, "aligned.fasta")
ggmsa("aligned.fasta", start = 1, end = 500)
```

---

## Buscar Patrones en Secuencias

```r
# Encuentra un patrón (ej: codon START "ATG")
matches <- vmatchPattern("ATG", sequences)

# Cuenta cuántas veces aparece
counts <- vcountPattern("ATG", sequences)

# Calcula frecuencia de patrones de 3 bases
freq <- oligonucleotideFrequency(sequences, width = 3)
```

---

## Flujos de Trabajo Típicos

### Flujo 1: Comparar Variantes (Lo básico)

```r
library(DECIPHER)
library(Biostrings)
library(ape)

# Carga todo
seqs <- readDNAStringSet("variants.fasta")

# Alinea
aligned <- AlignSeqs(seqs, verbose = TRUE)

# Calcula distancias
dist <- dist.hamming(aligned)

# Crea árbol
tree <- nj(as.dist(dist))

# Dibuja
plot(tree)
```

### Flujo 2: Procesar Muchos Archivos

```r
# Procesa varios archivos automáticamente
fasta_files <- list.files("data/", pattern = "\\.fasta$", full.names = TRUE)

for (file in fasta_files) {
  seqs <- readDNAStringSet(file)
  aligned <- AlignSeqs(seqs)
  
  # Guarda cada resultado
  base_name <- sub(".fasta", "", basename(file))
  saveRDS(aligned, paste0("results/", base_name, "_aligned.rds"))
}
```

### Flujo 3: Análisis de Variantes SARS-CoV-2

```r
# Para tu proyecto específicamente
aligned <- AlignSeqs(
  sequences,
  gapOpening = -50,      # Parámetros buenos para variantes
  gapExtension = -10,
  iterations = 2
)

# Mira mutaciones por posición
mutations_per_site <- colSums(aligned != aligned[[1]])
hist(mutations_per_site, main = "Mutaciones por posición")
```

---

## Tips para que Todo Vaya Más Rápido

```r
# Usa varios núcleos del computador
aligned <- AlignSeqs(sequences, processors = 4)

# Para archivos gigantes, alinea por partes
big_seqs <- readDNAStringSet("huge_file.fasta")
batch_size <- 50
for (i in seq(1, length(big_seqs), by = batch_size)) {
  end <- min(i + batch_size - 1, length(big_seqs))
  batch <- big_seqs[i:end]
  aligned_batch <- AlignSeqs(batch)
}

# Guarda resultados (así no repites el trabajo)
saveRDS(aligned, "cache/aligned.rds")
loaded <- readRDS("cache/aligned.rds")
```

---

## Si Algo Sale Mal

```r
# Verifica que esté instalado
library(DECIPHER)

# Revisa si las secuencias están bien formateadas
class(sequences)        # Debe ser "DNAStringSet"
length(sequences)       # Cuántas tienes
width(sequences)        # Longitud de cada una

# Limpia bases raras o inválidas
cleaned <- chartr("NRYSWKMBDHV", "N", sequences)

# Chequea memoria disponible
gc()
object.size(sequences)

# Atrapa errores sin que todo explote
tryCatch({
  aligned <- AlignSeqs(sequences)
}, error = function(e) {
  print(paste("Error:", e$message))
})
```

---

## Atajos de Teclado en VS Code (Mac)

| Qué hacer | Teclas |
|--------|---|
| Ejecutar una línea | Cmd + Enter |
| Ejecutar bloque de código | Cmd + Shift + Enter |
| Ejecutar todo el script | Cmd + Shift + R |
| Enviar a terminal | Ctrl + Alt + Enter |
| Abrir/cerrar terminal | Ctrl + ` |
| Buscar | Cmd + F |

---

## Librerías que Necesitas

```r
# OBLIGATORIO - Sin estas no funciona nada
library(DECIPHER)       # Alineamiento & análisis de secuencias
library(Biostrings)     # Objetos para manejar secuencias
library(ape)            # Árboles filogenéticos

# RECOMENDADO - Para que los gráficos se vean lindos
library(ggtree)         # Dibuja árboles bonitos
library(ggmsa)          # Visualiza alineamientos
library(ggplot2)        # Gráficos en general
library(seqinr)         # Lee FASTA (alternativa)
```

---

## Estructura del Proyecto (Cómo Está Organizado)

```
/Users/binivazquez/UniWorkspace/computational_biology/
├── classwork/
│   ├── globalSetup.R              ← Configuración compartida
│   ├── graficos.R                 ← Gráficos & DECIPHER
│   ├── basicsDNAan.Rmd            ← Análisis en RMarkdown
│   └── assets/
│       ├── wuhan.fasta            ← Secuencia de referencia
│       └── secuencias_fasta/      ← Tus secuencias
├── assignments/
│   └── analisis1_sp.R             ← Tu análisis
├── outputs/                       ← Donde guardar resultados
├── DECIPHER_VS_CODE_GUIDE.md      ← Guía completa (en inglés)
└── DECIPHER_CHEATSHEET_ES.md      ← Este archivo
```

---

## Trucos de Una Línea (Para Impresionar)

```r
# Todo en una sola línea: carga, alinea, calcula distancia, dibuja árbol
plot(nj(as.dist(dist.hamming(AlignSeqs(readDNAStringSet("seqs.fasta"))))))

# Mira qué hay en un archivo
readDNAStringSet("file.fasta")

# Calcula GC content rapidísimo
colMeans(letterFrequency(readDNAStringSet("seqs.fasta"), "GC")) / width(readDNAStringSet("seqs.fasta"))

# Encuentra todos los START codones
vmatchPattern("ATG", readDNAStringSet("seqs.fasta"))
```

---

## Cuándo Usar Qué Herramienta

| Necesitas | Usa | Alternativa |
|------|----------|-----------|
| Leer FASTA | `readDNAStringSet()` ✅ | `read.fasta()` |
| Alinear secuencias | **DECIPHER** ⭐ | MUSCLE, MAFFT |
| Calcular distancia | `dist.hamming()` | `dist.dna()` |
| Hacer árboles | `nj()`, `upgma()` (ape) | Figtree (programa externo) |
| Gráficos lindos | `ggtree()`, `ggmsa()` | Seqinr básico |

---

## Quick Reference (Lo Más Importante)

```r
# 1. Carga
seqs <- readDNAStringSet("file.fasta")

# 2. Alinea
aligned <- AlignSeqs(seqs, verbose = TRUE)

# 3. Distancia
dist <- dist.hamming(aligned)

# 4. Árbol
tree <- nj(as.dist(dist))

# 5. Dibuja
plot(tree)
```

---

## Pro Tips de Estudiante

💾 **Guarda frecuentemente**: Cada 5 minutos, hazlo
```r
save(aligned, tree, file = "my_analysis.RData")
```

🔄 **Reutiliza código**: Copia-pega lo que funcione

📝 **Comenta todo**: Explica qué hace cada línea

🐛 **Errores son normales**: Busca en `DECIPHER_FAQ_REFERENCE_ES.md`

🧬 **Los árboles filogenéticos son cool**: Disfruta lo que haces

---

## Cuando Estés Atascado

1. **Busca primero en:** 
   - `DECIPHER_FAQ_REFERENCE_ES.md` (preguntas frecuentes)
   - `DECIPHER_VS_CODE_GUIDE.md` (guía completa)

2. **Ejecuta esto:**
   ```r
   ?AlignSeqs  # Lee la ayuda de cualquier función
   ```

3. **Pregunta en:**
   - [Bioconductor Support](https://support.bioconductor.org/)
   - Stack Overflow con tags `[r]` `[bioconductor]`

---

**¡Imprime esto y pégalo en tu monitor!** 

Cada sección está diseñada para ser copiar-pega. Modifica las rutas de los archivos según TUS carpetas y listo.

**¡Dale, a codificar!** 🧬🚀
