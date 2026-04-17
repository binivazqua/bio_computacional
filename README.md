# 🧬 DECIPHER en Español - Para Estudiantes

**¡Bienvenido! Aquí está TODO lo que necesitas.**

---

## Documentos Disponibles

Referencia al trabajar en vs code:

| Documento                  | Para qué                       | Dificultad |
| -------------------------- | ------------------------------ | ---------- |
| **decipher_cheatsheet.md** | Lo más importante, copiar-pega | ⭐ Fácil   |
| **decipher_debugging.md**  | Problemas y cómo solucionarlos | ⭐ Fácil   |

---

## STARTER GUIDE

### Si NUNCA has usado DECIPHER:

1. **Lee esto primero:** `DECIPHER_CHEATSHEET_ES.md` (5 minutos)
2. **Luego copia código:** Prueba los ejemplos básicos
3. **Si algo falla:** Busca en `DECIPHER_FAQ_ES.md`
4. **¿Todavía no?** Pregunta a tus compañeros o al profe

### Si NECESITAS INSTALAR DECIPHER:

1. Lee el primer bloque de `DECIPHER_CHEATSHEET_ES.md`
2. Sigue paso a paso
3. Verifica que funcione

### Si TENGO UN ERROR:

1. Lee el mensaje de error completo
2. Busca tu error en `DECIPHER_FAQ_ES.md`
3. Copia la solución
4. Prueba de nuevo

---

## 💡 ¿Qué es DECIPHER? (Resumen Rápido)

DECIPHER es una herramienta de R para:

✅ **Alinear secuencias de ADN** - Compara múltiples secuencias  
✅ **Calcular distancias** - Qué tan diferentes son  
✅ **Hacer árboles filogenéticos** - Ver cómo están relacionadas  
✅ **Buscar mutaciones** - Encontrar cambios específicos

**Para tu proyecto:** Comparar variantes de SARS-CoV-2 y crear árboles

---

## 📂 Estructura del Proyecto

```
Tu carpeta/
├── classwork/
│   ├── globalSetup.R              ← Se carga automáticamente
│   ├── assets/
│   │   ├── wuhan.fasta            ← Secuencias de referencia
│   │   └── secuencias_fasta/      ← Tus secuencias
│   └── graficos.R                 ← Aquí va DECIPHER
│
├── assignments/
│   └── analisis1_sp.R             ← Tu código de análisis
│
├── outputs/                       ← Donde guardar resultados
│
├── DECIPHER_CHEATSHEET_ES.md      ← Copiar-pega código aquí
├── DECIPHER_FAQ_ES.md             ← Problemas & soluciones
└── README_ESTUDIANTES.md          ← Este archivo
```

---

## ⚡ Flujo Típico (Paso a Paso)

### Paso 1: Carga tus archivos

```r
library(DECIPHER)
seqs <- readDNAStringSet("mi_archivo.fasta")
```

### Paso 2: Alinea

```r
aligned <- AlignSeqs(seqs, verbose = TRUE)
```

### Paso 3: Calcula distancias

```r
library(ape)
dist <- dist.hamming(aligned)
```

### Paso 4: Crea árbol

```r
tree <- nj(as.dist(dist))
```

### Paso 5: Dibuja

```r
plot(tree)
```

### Paso 6: Guarda

```r
saveRDS(aligned, "mi_alineamiento.rds")
write.csv(dist, "distancias.csv")
```

---

## 🛠️ Instalación (Si Necesitas)

```r
# Instalación rápida
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DECIPHER")

# Verifica que funcione
library(DECIPHER)
packageVersion("DECIPHER")
```

Si falla algo, busca en `decipher_debugging.md` la sección "Instalación".

---

## 📖 Cómo Leer los Documentos

### DECIPHER_CHEATSHEET_ES.md

- **Mejor para:** Copiar código
- **Cómo usarlo:** Busca lo que necesitas, copia-pega
- **Tiempo:** 5 minutos para aprender
- **Tip:** Imprimelo y pega en tu monitor

### DECIPHER_FAQ_ES.md

- **Mejor para:** Cuando algo no funciona
- **Cómo usarlo:** Busca tu error, sigue la solución
- **Tiempo:** Varía según el problema
- **Tip:** Busca con Ctrl+F

### DECIPHER_VS_CODE_GUIDE.md (En inglés, pero completo)

- **Mejor para:** Entender TODO en detalle
- **Cómo usarlo:** Lee sección por sección
- **Tiempo:** 1-2 horas
- **Tip:** Solo si quieres aprender más

---

## 🎯 Tareas Comunes

### "Necesito alinear mis secuencias"

1. Abre `decipher_cheatsheet.md`
2. Ve a "Alineamiento"
3. Copia el código
4. Cambia el nombre del archivo
5. Ejecuta

### "Mi alineamiento tarda años"

1. Ve a `decipher_debugging.md`
2. Busca "El alineamiento tarda DEMASIADO"
3. Sigue los tips

### "Tengo un error y no sé qué significa"

1. Lee el error completo
2. Abre `DECIPHER_FAQ_ES.md`
3. Busca palabras clave del error
4. Sigue la solución

### "Necesito un árbol filogenético"

1. Abre `decipher_debugging.md`
2. Ve a "Construir Árbol Filogenético"
3. Copia-pega
4. Hecho

---

## 🚨 Si Algo Sale Mal

**Paso 1:** Lee el error completo (no tengas miedo)

**Paso 2:** Busca en `decipher_debugging.md`

**Paso 3:** Si no está:

```r
?AlignSeqs  # Para ver ayuda de una función
help(tu_funcion)
```

**Paso 4:** Googlea: `[r] [bioconductor] "tu error"`

**Paso 5:** Pregunta en clase (de verdad, no hay problema)

---

## 📊 Comandos Más Usados

```r
# Cargar
sequences <- readDNAStringSet("file.fasta")

# Alinear
aligned <- AlignSeqs(sequences, verbose = TRUE)

# Distancia
dist <- dist.hamming(aligned)

# Árbol
tree <- nj(as.dist(dist))

# Ver propiedades
width(sequences)      # Longitud
length(sequences)     # Cuántas hay
letterFrequency(sequences, "ACGT")  # Composición

# Guardar
saveRDS(aligned, "file.rds")
write.csv(dist, "file.csv")
```

---

## 👥 Pedir Ayuda

### A TUS COMPAÑEROS

Pregunta directamente, todos estamos aprendiendo

### AL PROFE

Muestra el error y qué intentaste

### ONLINE

- **Bioconductor:** https://support.bioconductor.org/
- **Stack Overflow:** Busca con `[r]` y `[bioconductor]`

---

## Pro Tips

🫵 **Comenta tu código** - Explica qué hace cada línea

🫵 **Guarda frecuente** - Cada 5 minutos, hazlo

🫵 **Usa GitHub** - Para versión control

🫵 **Haz pequeños tests** - Antes de procesos grandes

🫵 **Lee errors completos** - La solución está ahí

🫵 **Usa variables con nombres claros:**

```r
# Malo
d <- AlignSeqs(s)

# BIEN
aligned_sequences <- AlignSeqs(sequences, verbose = TRUE)
```

---

## 📋 Checklist Antes de Presentar

- [ ] Tu código funciona (lo probaste)
- [ ] Hay comentarios explicando qué hace
- [ ] Los gráficos se ven lindos
- [ ] Guardaste los resultados
- [ ] Documentaste parámetros usados
- [ ] Compartiste con tus compañeros
- [ ] Entiendes qué hace cada línea

---

## 🧬 Comandos Útiles (Copy-Paste Ready)

### Carga rápida de 3 archivos

```r
library(DECIPHER)
mexico <- readDNAStringSet("classwork/assets/secuencias_fasta/mexico.fasta")
wuhan <- readDNAStringSet("classwork/assets/wuhan.fasta")
francia <- readDNAStringSet("classwork/assets/secuencias_fasta/francia.fasta")

all_seqs <- c(mexico, wuhan, francia)
names(all_seqs) <- c("Mexico", "Wuhan", "Francia")
```

### Análisis completo en 10 líneas

```r
library(DECIPHER)
library(ape)

seqs <- readDNAStringSet("file.fasta")
aligned <- AlignSeqs(seqs, verbose = TRUE)
dist <- dist.hamming(aligned)
tree <- nj(as.dist(dist))
plot(tree)
title("Mi árbol filogenético")

saveRDS(aligned, "alignment.rds")
```

---

## 🎬 Próximos Pasos

1. **Hoy:** Lee `DECIPHER_CHEATSHEET_ES.md` (15 min)
2. **Mañana:** Copia los ejemplos y prueba (30 min)
3. **Después:** Analiza tus datos (depende)
4. **Cuando falle:** Busca en FAQ (5-10 min)
5. **Al final:** ¡Preséntalo a la clase!

---

## 📞 Preguntas Frecuentes Sobre Este README

**P: ¿Debo leer TODO?**  
R: No. Lee CHEATSHEET, experimenta, y busca en FAQ según necesites.

**P: ¿Puedo usar los documentos en inglés?**  
R: Sí, están disponibles. Pero los españoles son más chill.

**P: ¿Cuánto tarda un análisis típico?**  
R: De 5 minutos (pequeño) a 1 hora (grande). Depende del tamaño.

**P: ¿Puedo ver los gráficos mientras se procesa?**  
R: En VS Code sí. Abre la terminal con R y ejecución paso a paso.

---

_Creado: Abril 2026_  
_Para: Estudiantes de Bio Computacional_  
_Por: Bini Vázquez + Claude Haiku 4.5_
