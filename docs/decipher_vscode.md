Autoría: Biniza Vázquez y Claude Haiku 4.5
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

## ¿FASTA o GenBank? Cuál usar

| Aspecto         | FASTA                       | GenBank                          |
| --------------- | --------------------------- | -------------------------------- |
| **Complejidad** | Baja (5 líneas de código)   | Media (10+ líneas)               |
| **Velocidad**   | Rápida                      | Más lenta (requiere BD)          |
| **Metadata**    | Solo secuencia              | LOCUS, DEFINITION, AUTHORS, etc. |
| **Ideal para**  | Clase, demos, visualización | Investigación, publicaciones     |

**Recomendación:** Para trabajos en clase y visualización de árboles, **usa FASTA**. Es más simple, igual de efectivo, y genera árboles igual de buenos.

Ver documento completo: [FASTA_vs_GenBank.md](FASTA_vs_GenBank.md)

## Flujo 1: trabajar con FASTA en VS Code (RECOMENDADO)

Úsalo para alineamiento directo y árboles sencillos.
Obtenido de la documentación de DECIPHER actual + ayuda de Claude Code Haiku 4.5

```r
library(Biostrings)
library(DECIPHER)
library(ape)

# fasta_file <- "classwork/assets/secuencias_fasta/all_arn_sequences.fasta"

fasta_file <- "classwork/assets/secuencias_fasta/all_arn_sequences_CLEAN.fasta"

seqs <- readDNAStringSet(fasta_file)
aligned <- AlignSeqs(seqs, verbose = TRUE)
d <- DistanceMatrix(aligned, correction="K80", verbose=FALSE)  # ← Aquí está
tree <- nj(as.dist(d))



plot(tree)
```

## Hacer un Árbol estilo MSA sin que tarde años en correr...

```r
library(Biostrings)
library(DECIPHER)
library(ape)

# Configuración de archivos (recomendación en docs para guardar png)
fasta_file <- "classwork/assets/secuencias_fasta/all_arn_sequences.fasta"
output_dir <- "outputs"

# Crear directorio si no existe (failsafe, en cualquier contexto)
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}
# ÁRBOL FILOGENÉTICO (Neighbor-Joining + K80)

# Mejor para: Entender relaciones evolutivas entre variantes
# Interpreta: Tiempo evolutivo, ramas reflejan distancia genética

cat("\n========== OPCIÓN 1: ÁRBOL FILOGENÉTICO ==========\n")

# Preparar datos desde FASTA

seqs1 <- readDNAStringSet(fasta_file)

# Alinear primero

aligned1 <- AlignSeqs(seqs1, verbose = FALSE)

# Calcular distancias

d1 <- DistanceMatrix(aligned1, correction = "K80", verbose = FALSE)

# Construir árbol filogenético usando hclust para mejor control

hc1 <- hclust(as.dist(d1), method = "average")
tree_phylo <- as.phylo(hc1)

# Visualizar

png(file.path(output_dir, "arbol_filogenetico_asmsa.png"), width = 1400, height = 700)

# CALVE: "rightwards" lo hace estilo MSA
plot.phylo(tree_phylo,
    direction = "rightwards", cex = 0.7,
    main = "Árbol Filogenético (UPGMA + K80)\n28 variantes SARS-CoV-2"
)
dev.off() #sólo para MACOs

```

## Flujo 2: trabajar con GenBank en VS Code

Úsalo cuando conservas la estructura GenBank y quieres importar desde una base temporal con DECIPHER.

#### Ojo a la línea `tree <- DECIPHER::Treeline()`

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
d <- DistanceMatrix(dna, correction="K80", verbose=FALSE)
tree <- DECIPHER::Treeline(myXStringSet=dna, method="NJ", cutoff=0.05, showPlot=TRUE, verbose=FALSE)

plot(tree)
dbDisconnect(dbConn)
```

## Parámetro `correction="K80"` en `DistanceMatrix()`

La línea `d <- DistanceMatrix(dna, correction="K80", verbose=FALSE)` es fundamental. Calcula las **distancias evolutivas** entre secuencias antes de construir el árbol.

### ¿Qué es K80?

**K80** (Kimura 2-Parameter) es un **modelo de evolución molecular** que asume dos tasas diferentes:

- Transiciones (A-G, C-T): cambios dentro del mismo grupo químico
- Transversiones (A-C, A-T, G-C, G-T): cambios entre grupos químicos distintos

Las transversiones son menos frecuentes biológicamente, así que K80 las pondera más en el cálculo de distancia.

### Modelos alternativos

| Modelo                  | Comportamiento                               | Cuándo usarlo                               |
| ----------------------- | -------------------------------------------- | ------------------------------------------- |
| `"JC69"` (Jukes-Cantor) | Asume tasas iguales para todos los cambios   | Datos muy divergentes o cuando no hay sesgo |
| `"raw"`                 | Cuenta diferencias directas (sin corrección) | Solo para comparaciones rápidas             |
| `"NJ"`                  | Balanceado, ramas desiguales permitidas      | Por defecto, mejor para datos reales        |

### Impacto en el árbol

El modelo K80 produce distancias más **realistas**, lo que resulta en:

- **Árboles más precisos** en la topología
- **Ramas mejor calibradas** que reflejan tiempo evolutivo real
- **Clusters de variantes más confiables** (ej.: Ómicron BA.2 vs Ómicron Alemania)

Para SARS-CoV-2 en clase, **K80 es el estándar** porque:

- Respeta la biología molecular (transiciones vs transversiones)
- Funciona bien con secuencias no muy divergentes
- Es ampliamente aceptado en filogeografía viral

### Combinación recomendada

```r
# Paso 1: Calcular distancias con K80
d <- DistanceMatrix(dna, correction="K80", verbose=FALSE)

# Paso 2: Construir árbol con NJ y cutoff=0.05
tree <- DECIPHER::Treeline(
  myXStringSet=dna,
  method="NJ",
  cutoff=0.05,
  showPlot=TRUE,
  verbose=FALSE
)
```

Esta combinación (K80 + NJ + cutoff=0.05) es el **flujo estándar** para este proyecto.

## Parámetros de `TreeLine()` y su impacto visual

La línea `tree <- DECIPHER::Treeline()` controla la construcción y apariencia del árbol filogenético, si no usamos esta sintaxis, es muy probable que no corra adecuadamente.

### `method=""` - Algoritmo de construcción

Es el **parámetro clave** para la estructura del árbol:

| Método                    | Resultado Visual                                         | Cuándo usarlo                        |
| ------------------------- | -------------------------------------------------------- | ------------------------------------ |
| `"NJ"` (Neighbor-Joining) | Árbol **desbalanceado**, ramas con longitudes desiguales | Por defecto, mejor para datos reales |
| `"UPGMA"`                 | Árbol **muy simétrico**, ramas casi iguales              | Si asuces reloj molecular            |
| `"NJst"`                  | Árbol para datos **multisecuencia**                      | Casos especiales                     |

**Impacto visual:**

- `method="NJ"`: Las ramas pueden ser muy cortas o largas → árbol con estructura irregular
- `method="UPGMA"`: Las ramas son más simétricas → árbol visualmente ordenado

### `cutoff=0.05` - Umbral de agrupamiento (collapse)

Este parámetro **collapsa ramas cortas** e agrupa secuencias cercanas:

| Valor         | Efecto                   | Árbol resultante                                      |
| ------------- | ------------------------ | ----------------------------------------------------- |
| `cutoff=0.01` | Muy restrictivo          | **Más nodos separados**, ramificación densa           |
| `cutoff=0.05` | Intermedio (recomendado) | Agrupa variantes **casi idénticas**, estructura clara |
| `cutoff=0.1`  | Permisivo                | **Menos nodos**, más ramas colapsadas                 |

**Ejemplo visual:**

```
cutoff=0.01 (detallado)    cutoff=0.05 (balanceado)    cutoff=0.1 (simplificado)
    ├─ VarianteA             ├─ VarianteA┐               ├─ Grupo1┐
    ├─ VarianteB             ├─ VarianteB┤ (colapsadas)  ├─ Grupo2┤ (más colapsadas)
    ├─ VarianteC             ├─ VarianteC
    ├─ VarianteD             ├─ ...
```

### `showPlot=TRUE/FALSE` - Renderización inmediata

- `showPlot=TRUE`: **Dibuja el árbol inmediatamente** en la ventana gráfica de VS Code
- `showPlot=FALSE`: Calcula el árbol pero no lo dibuja (útil si planeas custom-plotearlo con `ggtree` después)

### `verbose=TRUE/FALSE` - Información de progreso

- `verbose=TRUE`: Imprime mensajes en la consola ("Aligning sequences...", "Computing distances...")
- `verbose=FALSE`: Ejecución silenciosa (mejor para scripts y documentación limpia)

## Ejemplos prácticos para la situación problema

### Opción 1: Árbol **detallado** (máxima ramificación)

Úsalo si necesitas ver todas las diferencias pequeñas entre variantes:

```r
tree <- DECIPHER::Treeline(
  myXStringSet=dna,
  method="NJ",
  cutoff=0.01,        # ← Bajo: máxima granularidad
  showPlot=TRUE,
  verbose=FALSE
)
```

### Opción 2: Árbol **balanceado** (recomendado para 28 variantes)

Úsalo como configuración estándar. Muestra grupos claros sin exceso de ruido:

```r
tree <- DECIPHER::Treeline(
  myXStringSet=dna,
  method="NJ",
  cutoff=0.05,        # ← Estándar: buen balance
  showPlot=TRUE,
  verbose=FALSE
)
```

### Opción 3: Árbol **simplificado** (mínima ramificación)

Úsalo si los datos son muy ruidosos o si necesitas destacar solo grupos principales:

```r
tree <- DECIPHER::Treeline(
  myXStringSet=dna,
  method="NJ",
  cutoff=0.1,         # ← Alto: menos ruido
  showPlot=TRUE,
  verbose=FALSE
)
```

### Opción 4: Calcular + custom-ploteo

Si quieres mejorar la visualización con `ggtree`:

```r
tree <- DECIPHER::Treeline(
  myXStringSet=dna,
  method="NJ",
  cutoff=0.05,
  showPlot=FALSE
)

library(ggtree)
ggtree(tree) +
  geom_tiplab(size=3) +
  geom_treescale() +
  theme_tree2()
```

## Recomendación para SARS-CoV-2

Con **28 variantes**, recomendamos:

- `method="NJ"`: estándar filogénico
- `cutoff=0.05`: agrupa variantes casi idénticas pero mantiene diferencias importantes visibles
- `showPlot=TRUE`: visualiza inmediatamente
- `verbose=FALSE`: ejecución limpia

Este balance permite ver clusters de variantes (p. ej., Ómicron BA.2 vs Ómicron Alemania) sin perder la estructura general del árbol.

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

---

## Fuentes y referencias oficiales

Esta documentación se basa en las fuentes oficiales de DECIPHER y paquetes relacionados:

### DECIPHER (Bioconductor)

- **Página oficial en Bioconductor:**  
  https://bioconductor.org/packages/release/bioc/html/DECIPHER.html

- **Manual de referencia (PDF):**  
  https://bioconductor.org/packages/release/bioc/manuals/DECIPHER/man/DECIPHER.pdf

- **Viñetas y tutoriales:**  
  https://bioconductor.org/packages/release/bioc/vignettes/DECIPHER/inst/doc/

- **GitHub (desarrollo):**  
  https://github.com/gaospecial/DECIPHER

### Paquetes relacionados

**Biostrings** (lectura de secuencias):

- https://bioconductor.org/packages/release/bioc/html/Biostrings.html
- https://bioconductor.org/packages/release/bioc/manuals/Biostrings/man/Biostrings.pdf

**ape** (árboles filogenéticos):

- https://CRAN.R-project.org/package=ape
- https://emmanuelparadis.github.io/ape/ape-introduction.pdf

**RSQLite** (bases de datos SQLite en R):

- https://CRAN.R-project.org/package=RSQLite
- https://rsqlite.r-dbi.org/

**ggtree** (visualización avanzada de árboles):

- https://bioconductor.org/packages/release/bioc/html/ggtree.html
- https://yulab-smu.top/treedata-book/

### Publicaciones académicas

- **DECIPHER: Accelerating automated sequence analysis for emerging human infectious diseases** (2013)  
  Wright, E. S. (2016). _Bioinformatics_, 32(12), 1877-1878.  
  https://doi.org/10.1093/bioinformatics/btw053

- **Phylogenetic Trees Made Easy: A How-To Manual**  
  Hall, B. G. (2013). _Sinauer Associates_.  
  https://www.sinauer.com/phylogenetic-trees-made-easy.html

### Recursos sobre métodos filogenéticos

- **Neighbor-Joining y algoritmos de construcción:**  
  Saitou, N., & Nei, M. (1987). The neighbor-joining method: A new method for reconstructing phylogenetic trees. _Molecular Biology and Evolution_, 4(4), 406-425.

- **Modelos de distancia (K80, JC69):**  
  Kimura, M. (1980). A simple method for estimating evolutionary rates of base substitutions through comparative studies of nucleotide sequences. _Journal of Molecular Evolution_, 16(2), 111-120.

### Herramientas complementarias

- **NCBI GenBank:**  
  https://www.ncbi.nlm.nih.gov/genbank/

- **FASTA format specification:**  
  https://en.wikipedia.org/wiki/FASTA_format

- **GenBank format specification:**  
  https://www.ncbi.nlm.nih.gov/genbank/

---

**Nota:** Esta documentación fue elaborada con apoyo de Claude Haiku 4.5 y está diseñada para complementar (no reemplazar) la documentación oficial de DECIPHER y sus paquetes dependientes.
