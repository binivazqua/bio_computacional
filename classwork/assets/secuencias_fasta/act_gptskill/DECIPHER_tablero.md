# DECIPHER: Tablero de Categorías Clave

**Referencia rápida para flujos de trabajo filogenéticos**  
_Actualizado: Abril 2026 | Recurso: Bioconductor DECIPHER + Biostrings + ape_
_Por: Biniza Vázquez A01737294 @Claude Haiku 4.5_

---

## 1. LECTURA DE DATOS

| Formato        | Función              | Entrada          | Salida                     | Complejidad |
| -------------- | -------------------- | ---------------- | -------------------------- | ----------- |
| **FASTA**      | `readDNAStringSet()` | `.fasta` / `.fa` | `DNAStringSet`             | ⭐ Baja     |
| **GenBank**    | `Seqs2DB()`          | `.gb` / `.gen`   | BD SQLite + `DNAStringSet` | ⭐⭐⭐ Alta |
| **Validación** | `file.exists()`      | Ruta de archivo  | Lógico (TRUE/FALSE)        | ⭐ Baja     |

**Decisión rápida:**

- ¿Archivo comienza con `>`? → **FASTA**
- ¿Archivo comienza con `LOCUS`? → **GenBank**

---

## 2. ALINEAMIENTO DE SECUENCIAS

| Función               | Parámetros Clave     | Salida                  | Velocidad | Cuándo usarlo                   |
| --------------------- | -------------------- | ----------------------- | --------- | ------------------------------- |
| `AlignSeqs()`         | `verbose=TRUE/FALSE` | `DNAStringSet` alineado | Rápida    | Siempre (estándar)              |
| `ConsensusSequence()` | `threshold=0.5`      | Secuencia consenso      | Rápida    | Después de alineamiento         |
| `RemoveGaps()`        | `removeGaps="all"`   | Sin gaps                | Rápida    | Si necesitas secuencias limpias |

**Flujo típico:**

```
Lectura → AlignSeqs() → [opcional] RemoveGaps() → Distancias
```

---

## 3. CÁLCULO DE DISTANCIAS

| Modelo                  | Ecuación       | Ventaja                 | Desventaja                         | Ideal para             |
| ----------------------- | -------------- | ----------------------- | ---------------------------------- | ---------------------- |
| **K80** (Kimura)        | 2 parámetros   | Biológicamente realista | Requiere divergencia baja          | ⭐⭐⭐ **RECOMENDADO** |
| **JC69** (Jukes-Cantor) | 1 parámetro    | Simple, uniforme        | Ignora transiciones/transversiones | Datos muy divergentes  |
| **raw**                 | Conteo directo | Inmediato               | Sin corrección evolutiva           | Comparaciones rápidas  |

**Función:** `DistanceMatrix(alignment, correction="K80")`

---

## 4. CONSTRUCCIÓN DEL ÁRBOL

| Método                    | Estructura     | Uso                   | Velocidad | Topología                  |
| ------------------------- | -------------- | --------------------- | --------- | -------------------------- |
| **NJ** (Neighbor-Joining) | Desbalanceado  | ⭐⭐⭐ Estándar       | Rápida    | Realista (ramas variables) |
| **UPGMA**                 | Muy simétrico  | Asume reloj molecular | Rápida    | Simétrico (artificial)     |
| **ML** (Max. Likelihood)  | Probabilístico | Más preciso (lento)   | Muy lenta | Altamente realista         |
| **hclust**                | Dendrograma    | Control total         | Rápida    | Flexible                   |

**Flujo recomendado:**

```r
d <- DistanceMatrix(alignment, correction="K80")
hc <- hclust(as.dist(d), method="average")
tree <- as.phylo(hc)
```

---

## 5. PARÁMETROS DE VISUALIZACIÓN

### `plot.phylo()` - Opciones principales

| Parámetro      | Valores                     | Efecto                | Ejemplo                  |
| -------------- | --------------------------- | --------------------- | ------------------------ |
| `direction`    | `"rightwards"`, `"upwards"` | Orientación del árbol | `direction="rightwards"` |
| `cex`          | Número 0.5–2                | Tamaño de etiquetas   | `cex=0.7`                |
| `main`         | String                      | Título del gráfico    | `main="Mi árbol"`        |
| `label.offset` | Número                      | Distancia label-rama  | `label.offset=0.01`      |
| `edge.width`   | Número                      | Grosor de ramas       | `edge.width=1.5`         |

### Treeline() en DECIPHER (alternativa)

| Parámetro  | Rango                       | Efecto                 | Recomendación       |
| ---------- | --------------------------- | ---------------------- | ------------------- |
| `method`   | `"NJ"`, `"UPGMA"`, `"NJst"` | Algoritmo              | `"NJ"`              |
| `cutoff`   | 0.01–0.1                    | Colapsar ramas cortas  | `0.05` (equilibrio) |
| `showPlot` | `TRUE`, `FALSE`             | Dibujar inmediatamente | `TRUE`              |
| `verbose`  | `TRUE`, `FALSE`             | Mensajes en consola    | `FALSE`             |

**Ejemplos visuales:**

- `cutoff=0.01` → Árbol detallado, muchos nodos
- `cutoff=0.05` → Árbol balanceado **(recomendado)**
- `cutoff=0.1` → Árbol simplificado, pocos nodos

---

## 6. FUNCIONES AUXILIARES

| Función          | Propósito                     | Entrada             | Salida          |
| ---------------- | ----------------------------- | ------------------- | --------------- |
| `SearchDB()`     | Recuperar secuencias de BD    | BD SQLite           | `DNAStringSet`  |
| `dbConnect()`    | Conectar BD SQLite            | Parámetros          | Conexión activa |
| `dbDisconnect()` | Cerrar BD                     | Conexión            | —               |
| `as.phylo()`     | Convertir dendrograma a árbol | `hclust`            | `phylo`         |
| `nj()` (ape)     | Neighbor-Joining directo      | Matriz de distancia | `phylo`         |

---

## 7. FLUJOS COMPLETOS

### Flujo A: FASTA (RÁPIDO y SIMPLE) ⭐ RECOMENDADO

```r
library(Biostrings)
library(DECIPHER)
library(ape)

# 1. Lectura
seqs <- readDNAStringSet("ruta/al/archivo.fasta")

# 2. Alineamiento
alignment <- AlignSeqs(seqs, verbose=FALSE)

# 3. Distancias
d <- DistanceMatrix(alignment, correction="K80", verbose=FALSE)

# 4. Árbol
hc <- hclust(as.dist(d), method="average")
tree <- as.phylo(hc)

# 5. Visualización
plot.phylo(tree, direction="rightwards", cex=0.7, main="Mi árbol")
```

**Tiempo:** ~10 segundos (28 variantes)

---

### Flujo B: GenBank (COMPLETO y TRAZABLE)

```r
library(DBI)
library(RSQLite)
library(DECIPHER)
library(ape)

# 1. Conexión BD
dbConn <- dbConnect(RSQLite::SQLite(), ":memory:")

# 2. Cargar GenBank
Seqs2DB("ruta/archivo.gen", type="GenBank", dbConn=dbConn, identifier="id")

# 3. Recuperar secuencias
dna <- SearchDB(dbConn)

# 4. Distancias
d <- DistanceMatrix(dna, correction="K80", verbose=FALSE)

# 5. Árbol
tree <- DECIPHER::Treeline(myXStringSet=dna, method="NJ", cutoff=0.05, showPlot=TRUE)

# 6. Limpiar
dbDisconnect(dbConn)
```

**Tiempo:** ~20 segundos (28 variantes)

---

### Flujo C: Visualización Avanzada (con `ggtree`)

```r
library(ggtree)

# (usar Flujo A o B hasta obtener 'tree')

ggtree(tree) +
  geom_tiplab(size=3) +
  geom_treescale() +
  theme_tree2() +
  labs(title="Árbol SARS-CoV-2")
```

---

## 8. MATRIZ DE DECISIÓN

```
¿Qué archivo tengo?
│
├─→ FASTA (.fasta, .fa)
│   ├─→ readDNAStringSet()
│   ├─→ AlignSeqs()
│   ├─→ DistanceMatrix(correction="K80")
│   └─→ hclust() + as.phylo() ✓ RÁPIDO
│
└─→ GenBank (.gb, .gen)
    ├─→ dbConnect(RSQLite::SQLite())
    ├─→ Seqs2DB(type="GenBank")
    ├─→ SearchDB()
    ├─→ DistanceMatrix(correction="K80")
    └─→ Treeline() ✓ METADATOS
```

---

## 9. ERRORES COMUNES Y SOLUCIONES

| Error                                 | Causa                              | Solución                                           |
| ------------------------------------- | ---------------------------------- | -------------------------------------------------- |
| `">" expected at beginning of line 1` | Usando FASTA en archivo GenBank    | Cambiar a `Seqs2DB(..., type="GenBank")`           |
| `object 'tree' not found`             | `tree` no se creó                  | Revisar que `Treeline()` o `as.phylo()` se ejecutó |
| Etiquetas cortadas en gráfico         | Falta `direction="rightwards"`     | Agregar parámetro a `plot.phylo()`                 |
| Árbol muy amontonado                  | `cex` muy grande o sin `direction` | Reducir `cex` o cambiar dirección                  |
| BD no se conecta                      | Puerto SQLite bloqueado            | Usar `:memory:` en lugar de archivo                |

---

## 10. CHECKLIST DE CALIDAD

- [ ] Archivo leído correctamente (`length(seqs)` > 0)
- [ ] Alineamiento sin errores (`class(alignment)` = "DNAStringSet")
- [ ] Distancias calculadas (`class(d)` = "dist")
- [ ] Árbol generado sin NA (`all(!is.na(tree$edge))`)
- [ ] Gráfico legible (`cex ≤ 0.7`, `direction="rightwards"`)
- [ ] Títulos descriptivos agregados
- [ ] Guardado en archivo PNG/PDF (opcional)

---

## 11. REFERENCIAS RÁPIDAS

**Librerías necesarias:**

```r
install.packages("ape")
BiocManager::install("Biostrings")
BiocManager::install("DECIPHER")
install.packages("RSQLite")
BiocManager::install("ggtree")
```

**Documentación oficial:**

- DECIPHER: https://bioconductor.org/packages/release/bioc/html/DECIPHER.html
- ape: https://cran.r-project.org/package=ape
- Biostrings: https://bioconductor.org/packages/release/bioc/html/Biostrings.html

**Archivos del proyecto:**

```
classwork/assets/
├── secuencias_fasta/
│   ├── all_arn_sequences_CLEAN.fasta
│   └── act_gptskill/sequences_gptskill.fasta
└── secuencias_gb/
    ├── all_sequences.gen
    └── *.gb (individuales)
```

---

**Última actualización:** 24 de abril de 2026  
**Responsable:** Biniza Vázquez  
**Basado en:** DECIPHER Bioconductor + documentación oficial
