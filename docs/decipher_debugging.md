# DECIPHER FAQ - Preguntas que Todos Hacen

**¿Problema? ¡Dale, está aquí la solución!**

---

## 🔧 Instalación y Setup

### P: "No puedo instalar DECIPHER, dice que falta BiocManager"

**R:** Primero instala BiocManager:
```r
install.packages("BiocManager")
BiocManager::install("DECIPHER")
```

### P: "¿Qué versión de DECIPHER tengo?"

**R:** 
```r
library(DECIPHER)
packageVersion("DECIPHER")
```

---

## 📂 Cargar Secuencias

### P: "¿Debo usar `readDNAStringSet()` o `read.fasta()`?"

**R:** 
- **`readDNAStringSet()`** → Para DECIPHER ✅
- **`read.fasta()`** → Solo si vienes de seqinr

Si usas `read.fasta()`, convierte así:
```r
seq_list <- read.fasta("file.fasta", as.string = TRUE)
seq_set <- DNAStringSet(unlist(seq_list))
```

### P: "Mi archivo FASTA tiene muchas secuencias. ¿Está bien?"

**R:** Claro, DECIPHER lo maneja sin problemas:
```r
sequences <- readDNAStringSet("file.fasta")
length(sequences)  # Te dice cuántas hay
```

### P: "¿Puedo leer archivos directamente de internet?"

**R:** Sí, funciona igual:
```r
url <- "https://example.com/sequences.fasta"
sequences <- readDNAStringSet(url)
```

---

## ⚡ El Alineamiento (Lo Importante)

### P: "El alineamiento tarda DEMASIADO. ¿Qué hago?"

**R:** Tienes varias opciones:

```r
# Opción 1: Usa múltiples núcleos (lo más rápido)
aligned <- AlignSeqs(sequences, processors = 4)

# Opción 2: Menos iteraciones (menos preciso pero más rápido)
aligned <- AlignSeqs(sequences, iterations = 1)

# Opción 3: Alinea por partes
# (Mira DECIPHER_VS_CODE_GUIDE.md)
```

### P: "Error: 'could not find function AlignSeqs'"

**R:** No cargaste DECIPHER:
```r
library(DECIPHER)
```

### P: "El alineamiento sale raro. Hay gaps (guiones) por todo lados"

**R:** Probables causas:

```r
# 1. Verifica que las secuencias sean válidas
width(sequences)  # Todas diferentes = problema

# 2. Limpia secuencias raras
cleaned <- chartr("NRYSWKMBDHV", "N", sequences)

# 3. Intenta con menos iteraciones
aligned <- AlignSeqs(sequences, iterations = 1)
```

### P: "¿Puedo usar una secuencia de referencia?"

**R:** DECIPHER no lo hace directamente, pero puedes:
```r
# Pon la referencia primero
reference <- readDNAStringSet("reference.fasta")
queries <- readDNAStringSet("queries.fasta")
all_seqs <- c(reference, queries)

aligned <- AlignSeqs(all_seqs)
# La referencia estará primera en el resultado
```

### P: "¿Cómo veo el alineamiento después de hacerlo?"

**R:** Tienes varias formas:

```r
# Forma 1: En la consola
print(aligned)

# Forma 2: Guardarlo y verlo en FASTA
writeXStringSet(aligned, "alignment.fasta")

# Forma 3: Gráfico bonito (primero guarda)
library(ggmsa)
ggmsa("alignment.fasta", start = 1, end = 100)
```

---

## 🌳 Árboles Filogenéticos

### P: "¿Cómo calculo distancia entre mis secuencias?"

**R:** Usa la función de distancia:
```r
library(ape)

aligned <- AlignSeqs(sequences)
dist_matrix <- dist.hamming(aligned)  # Lo más simple

# O si quieres ser más fancy:
dist_matrix <- dist.dna(aligned, model = "JC69")
```

### P: "¿Qué diferencia hay entre Hamming y JC69?"

**R:**
- **Hamming** - Solo cuenta diferencias (simple pero rápido)
- **Jukes-Cantor (JC69)** - Considera mutaciones múltiples (mejor)
- **Kimura** - Distingue transiciones de transversiones
- **Tamura-Nei** - La más precisa pero lenta

Para tu proyecto SARS-CoV-2, **Hamming está bien**.

### P: "¿Cómo hago un árbol de mis secuencias?"

**R:**
```r
library(ape)

aligned <- AlignSeqs(sequences)
dist_mat <- dist.hamming(aligned)

# Neighbor-Joining (lo más común)
tree <- nj(as.dist(dist_mat))

# O UPGMA (si asumes un reloj molecular)
tree <- upgma(dist_mat)

# Dibuja
plot(tree)
axisPhylo()  # Agrega escala
```

### P: "¿Cómo agrego valores de bootstrap al árbol?"

**R:**
```r
library(ape)

boot_func <- function(x) {
  d <- dist.hamming(x)
  nj(as.dist(d))
}

boot_trees <- boot.phylo(tree, aligned, boot_func, rooted = TRUE)

plot(tree)
nodelabels(boot_trees, cex = 0.7)  # Muestra valores en nodos
```

---

## 🎨 Visualización (Gráficos)

### P: "¿Cómo hago un árbol bonito con ggtree?"

**R:**
```r
library(ggtree)

ggtree(tree) +
  geom_tiplab() +                    # Nombres en puntitas
  theme_tree2() +                    # Tema limpio
  labs(title = "Mi árbol")           # Título
```

O más elaborado:
```r
ggtree(tree, layout = "rectangular") +
  geom_tiplab(aes(color = ifelse(grepl("mexico", label), "red", "blue"))) +
  theme_tree2() +
  theme(legend.position = "right")
```

### P: "¿Cómo visualizo mi alineamiento?"

**R:**
```r
library(ggmsa)

writeXStringSet(aligned, "temp_alignment.fasta")

ggmsa("temp_alignment.fasta", 
      start = 1, 
      end = 200,
      color = "Clustal")

unlink("temp_alignment.fasta")  # Limpia
```

### P: "¿Cómo guardo una figura para mi reporte?"

**R:**
```r
# Como PNG (resolución alta)
png("my_tree.png", width = 3000, height = 2000, res = 300)
plot(tree)
dev.off()

# O como PDF (mejor para publicar)
pdf("my_tree.pdf", width = 10, height = 10)
ggtree(tree) + geom_tiplab() + theme_tree2()
dev.off()
```

---

## 🚀 Performance (Cuando Todo es MUY Lento)

### P: "Tengo un archivo ENORME (>1GB). ¿Qué hago?"

**R:** Procesa por partes:
```r
big_seqs <- readDNAStringSet("huge_file.fasta")
batch_size <- 50

for (i in seq(1, length(big_seqs), by = batch_size)) {
  end <- min(i + batch_size - 1, length(big_seqs))
  batch <- big_seqs[i:end]
  aligned <- AlignSeqs(batch)
  
  # Guarda este batch
  saveRDS(aligned, paste0("batch_", i, ".rds"))
}
```

### P: "¿Cuánta memoria está usando?"

**R:**
```r
# Ve el tamaño
object.size(sequences)      # En bytes
print(object.size(sequences), units = "MB")

# Limpia memoria
gc()

# Regla: 1 millón bp ≈ 1 MB
```

---

## 🐛 Errores Comunes

### P: "Error: 'no valid letters in x'"

**R:** Hay caracteres raros en tus secuencias:
```r
# Mira cuál es el problema
invalid <- sequences[!grepl("^[ACGTN-]*$", as.character(sequences))]

# Limpia todo
cleaned <- chartr("NRYSWKMBDHV", "N", sequences)

# O sé más estricto (solo ACGT)
strict <- chartr("nryswkmbdhv", "A", sequences)
```

### P: "Error: 'subscript out of bounds' en AlignSeqs"

**R:** Probablemente secuencias vacías o raras:
```r
# Verifica
length(sequences)           # Al menos 2
width(sequences)            # Todas > 0
any(width(sequences) == 0)  # No debe haber vacías

# Intenta con debug
tryCatch({
  AlignSeqs(sequences, verbose = TRUE)
}, error = function(e) print(e))
```

### P: "Error: 'cannot coerce... to DNAStringSet'"

**R:** Formato incorrecto:
```r
class(sequences)  # Debe ser "DNAStringSet"

# Convierte si falta
if (class(sequences) == "list") {
  sequences <- DNAStringSet(sequences)
}

if (class(sequences) == "character") {
  sequences <- DNAStringSet(sequences)
}
```

### P: "R se congela durante el alineamiento"

**R:** Se quedó sin memoria:
```r
# Reinicia R limpio
rm(list = ls())
gc()

# Intenta con 1 núcleo (menos RAM)
aligned <- AlignSeqs(sequences, processors = 1)

# O reduce el tamaño
sequences <- sequences[1:50]  # Prueba primero
```

---

## 📊 Organización (Dónde Poner Las Cosas)

### P: "¿Dónde meto mis archivos FASTA?"

**R:** Así está organizado:
```
computational_biology/
├── classwork/
│   └── assets/
│       ├── wuhan.fasta          ← FASTA aquí
│       └── secuencias_fasta/    ← O aquí
└── outputs/                     ← Resultados aquí
    └── alignments/
```

### P: "Tengo múltiples proyectos. ¿Cómo los organizo?"

**R:**
```
computational_biology/
├── proyecto1/
│   ├── data/
│   ├── scripts/
│   └── outputs/
├── proyecto2/
│   └── ...
└── shared/
    ├── globalSetup.R
    └── functions/
```

### P: "¿Cómo aseguro que mi análisis sea reproducible?"

**R:**
```r
# Guarda la sesión
writeLines(
  capture.output(sessionInfo()),
  "outputs/session_info.txt"
)

# Documenta parámetros
params <- list(
  method = "AlignSeqs",
  gap_opening = -25,
  iterations = 2,
  date = Sys.Date()
)

# Guarda todo junto
save(aligned, tree, params, file = "outputs/analysis.RData")
```

---

## 🧬 Tu Proyecto Específico (SARS-CoV-2)

### P: "¿Hay parámetros especiales para variantes de COVID?"

**R:** Sí, usa estos:
```r
# Para variantes (secuencias muy similares)
aligned <- AlignSeqs(
  sequences,
  gapOpening = -50,      # Caro abrir gaps
  gapExtension = -10,    # Caro extenderlos
  iterations = 2
)

# Usa Hamming distance (es suficiente)
dist <- dist.hamming(aligned)
```

### P: "¿Qué regiones del genoma debería analizar?"

**R:** El genoma de SARS-CoV-2 tiene zonas clave:
```r
# Regiones importantes
orf1 <- c(1, 21563)              # ORF1a/b
spike <- c(21563, 25384)         # S protein (mutaciones comunes)
nucleocapsid <- c(28274, 29533)  # N protein

# Extrae y analiza
spike_seqs <- subseq(sequences, spike[1], spike[2])
aligned_spike <- AlignSeqs(spike_seqs)
```

### P: "¿Cómo identifico mutaciones específicas?"

**R:**
```r
# Compara contra la referencia (Wuhan)
ref <- aligned[[1]]
query <- aligned[[2]]

# Encuentra diferencias
diffs <- which(ref != query)

# Ve qué cambió
data.frame(
  Position = diffs,
  Reference = ref[diffs],
  Query = query[diffs]
)
```

---

## 👥 Trabajo en Equipo

### P: "¿Cómo comparto mi análisis con el equipo?"

**R:**
```bash
# Via Git
git add classwork/
git commit -m "Análisis de variantes con DECIPHER"
git push

# Guarda en formatos compartibles
write.csv(results, "outputs/results.csv", row.names = FALSE)
```

### P: "¿Cómo documento el análisis para que otros lo entiendan?"

**R:** Usa R Markdown y comenta bien:
```r
# ===== ANÁLISIS DE VARIANTES =====
# Propósito: Comparar secuencias de COVID
# Autor: Tu nombre
# Fecha: 2024-04-17

# 1. Cargamos datos
seqs <- readDNAStringSet("data.fasta")

# 2. Alineamos (esto toma tiempo)
aligned <- AlignSeqs(seqs, verbose = TRUE)

# 3. Calculamos distancia
dist <- dist.hamming(aligned)

# Resultado: Encontramos X mutaciones diferentes
```

---

## 🎓 Tips Finales

### "¿Por dónde empiezo si soy principiante?"

1. Lee `DECIPHER_CHEATSHEET_ES.md` (este archivo)
2. Copia los ejemplos y prueba
3. Cuando se rompa, busca aquí en el FAQ
4. Pregunta en clase (no hay vergüenza)

### "¿Cuál es el flujo típico de un análisis?"

1. **Carga** → `readDNAStringSet()`
2. **Limpia** → Verifica secuencias válidas
3. **Alinea** → `AlignSeqs()`
4. **Calcula** → `dist.hamming()`
5. **Árbol** → `nj()` o `upgma()`
6. **Visualiza** → `plot()` o `ggtree()`
7. **Guarda** → `saveRDS()` y/o `write.csv()`

### "¿Qué si todo se va a la mierda?"

1. Respira hondo
2. Lee el error completo
3. Busca aquí en el FAQ
4. Si no está, googlea: `[r] [bioconductor] "tu error"`
5. Pregunta en clase

---

## 📚 Recursos Rápidos

**Dentro de R:**
```r
?AlignSeqs           # Ayuda de cualquier función
help(AlignSeqs)
sessionInfo()        # Info de tu versión
```

**Online:**
- [Bioconductor Support](https://support.bioconductor.org/)
- [Stack Overflow con tags [r] [bioconductor]](https://stackoverflow.com/questions/tagged/bioconductor)
- [DECIPHER Manual](https://www.bioconductor.org/packages/DECIPHER)

---

## 🚨 Resumen: Lo Más Importante

| Necesitas | Código |
|-----------|--------|
| Cargar FASTA | `readDNAStringSet("file.fasta")` |
| Alinear | `AlignSeqs(seqs, verbose=TRUE)` |
| Distancia | `dist.hamming(aligned)` |
| Árbol | `nj(as.dist(dist))` |
| Dibujar | `plot(tree)` |
| Guardar | `saveRDS(object, "file.rds")` |

---

**¿Tu problema no está aquí?**

1. Chequea `DECIPHER_VS_CODE_GUIDE.md` (guía completa)
2. Googlea el error
3. Pregunta a tus compañeros
4. Pregunta al profe

**¡Mucha suerte con tu análisis!** 🧬🚀

*Última actualización: Abril 2026*
*Para: Estudiantes de Biología Computacional*
