Autor: Biniza Vázquez  
Realizado con apoyo de Claude Haiku 4.5  
Abril, 2026

---

# Guía de Configuración de VS Code para DECIPHER

Procedimiento paso a paso para configurar tu ambiente de desarrollo de análisis de secuencias.

---

## Requisitos previos

Asegúrate de tener instalado:

- ✅ VS Code (descargar de [aquí](https://code.visualstudio.com/))
- ✅ R en tu Mac (`/usr/local/bin/R`)
- ✅ Herramientas de línea de comandos (XCode)

Si no tienes R, descárgalo de [CRAN](https://cran.r-project.org/bin/macosx/) e instala.

---

## Paso 1: Instala Extensiones en VS Code

Abre VS Code y ve a extensiones (Cmd + Shift + X).

### Extensiones Obligatorias

**1. R** (la oficial)
- Busca: "R"
- Instala: REditorSupport.r
- Con esto escribes y ejecutas código R

**2. Git Lens** (opcional pero útil)
- Busca: "Git Lens"
- Te ayuda a versionar código

**3. Markdown All in One** (opcional pero recomendado)
- Para editar documentación como esta

### Cómo Instalar

```
Cmd + Shift + X  →  Escribe el nombre  →  Click en Install
```

Listo. VS Code te avisa cuando esté instalado.

---

## Paso 2: Configura VS Code

### Abre los Settings

```
Cmd + ,  (acceso rápido)
```

O ve a: `Code > Preferences > Settings`

### Encuentra "R" en los Settings

En la barra de búsqueda, escribe `R`.

### Agrega Estos Valores

Busca o copia esto en JSON (va al final del archivo):

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

**¿Qué significa?**
- `defaultFormatter`: Formatea automáticamente tu código
- `formatOnSave`: Limpia cuando guardas
- `wordWrap`: No se sale del editor
- `sessionWatcher`: Monitorea tu sesión de R
- `useHttpgd`: Los gráficos se ven en VS Code
- `rterm.mac`: Dónde encontrar R en tu Mac

---

## Paso 3: Verifica que R Esté Instalado

Abre Terminal (Cmd + ~) y escribe:

```bash
which R              # Te dice dónde está R
R --version          # Versión de R
Rscript --version    # Versión de Rscript
```

Si dice `/usr/local/bin/R`, ¡perfecto!

Si dice "command not found", instala R desde [CRAN](https://cran.r-project.org/bin/macosx/).

---

## Paso 4: Primera Sesión de R en VS Code

### Crea una Terminal de R

Abre Command Palette:

```
Cmd + Shift + P
```

Escribe:

```
R: Create R Terminal
```

Y presiona Enter.

Se abrirá una ventana donde verás `>`. ¡Eso es R!

### Prueba Escribiendo

```r
2 + 2
```

Presiona Enter. Debe decir `[1] 4`. Si dice eso, R funciona.

---

## Paso 5: Instala DECIPHER y Paquetes

En la terminal de R que acabas de abrir, escribe:

### Primero: BiocManager

```r
install.packages("BiocManager")
```

Espera a que termine.

### Segundo: DECIPHER y Otros Paquetes

```r
# DECIPHER (lo importante)
BiocManager::install("DECIPHER")

# Otros bioinformáticos
BiocManager::install(c("Biostrings", "ape", "ggtree", "ggmsa"))

# Paquetes normales de CRAN
install.packages(c("seqinr", "ggplot2", "dplyr", "tidyr"))
```

### Verifica que Todo Funcione

```r
library(DECIPHER)
library(Biostrings)
library(ape)

packageVersion("DECIPHER")
```

Debe mostrar algo como: `[1] '2.24.0'`

Si todo dice `TRUE` o muestra versiones, ¡ganaste! 🎉

---

## Paso 6: Estructura de Carpetas

Crea esto en tu proyecto. Abre Terminal (Cmd + ~) y escribe:

```bash
cd /Users/binivazquez/UniWorkspace/computational_biology

# Crea las carpetas
mkdir -p classwork/assets/secuencias_fasta
mkdir -p classwork/assets/secuencias_gen
mkdir -p assignments
mkdir -p outputs
mkdir -p cache
mkdir -p functions

# Crea archivo README
touch README.md
```

Ahora tu proyecto se ve así:

```
computational_biology/
├── classwork/
│   ├── assets/
│   │   ├── secuencias_fasta/    (tus FASTA aquí)
│   │   └── secuencias_gen/      (archivos GenBank aquí)
│   └── globalSetup.R
├── assignments/
├── outputs/                      (resultados aquí)
├── cache/                        (archivos temporales)
├── functions/                    (código reutilizable)
└── README.md
```

---

## Paso 7: Crea .gitignore

Para que Git no suba archivos grandes. En VS Code:

```
Cmd + Shift + P  →  "Create File"  →  ".gitignore"
```

Copia esto:

```
# Archivos grandes
*.fasta
*.fa
*.gb
*.gbk

# Resultados de R
*.png
*.pdf
*.html
.Rhistory
.RData

# Caché
cache/
outputs/temp/

# IDE
.Rproj.user/
.vscode/

# Sistema
.DS_Store
Thumbs.db
```

---

## Paso 8: Crea globalSetup.R

Este archivo tiene todo lo que necesitas en cada análisis.

En VS Code:

```
Cmd + N  →  Guarda como "globalSetup.R"  →  Pon en classwork/
```

Copia esto:

```r
# globalSetup.R - Configuración Compartida
# Este archivo se ejecuta al principio de cada análisis

# Función para cargar paquetes seguro
load_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Instala con: BiocManager::install('%s')", pkg))
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# Carga paquetes principales
load_package("seqinr")
load_package("ape")
load_package("ggplot2")
load_package("dplyr")
load_package("Biostrings")

# Intenta cargar opcionales (no falles si no están)
tryCatch(load_package("DECIPHER"), error = function(e) warning(e))
tryCatch(load_package("ggtree"), error = function(e) warning(e))
tryCatch(load_package("ggmsa"), error = function(e) warning(e))

# Muestra información
cat("\n=== Setup de Biología Computacional ===\n")
cat("Directorio de trabajo:", getwd(), "\n")
cat("R versión:", R.version$version.string, "\n")
cat("✅ Configuración lista\n\n")
```

Guarda con Cmd + S.

---

## Paso 9: Crea Tu Primer Script de Prueba

Nuevo archivo:

```
Cmd + N  →  Guarda como "test_decipher.R"
```

Copia esto:

```r
# test_decipher.R - Prueba que todo funcione

# Configura directorio
setwd("/Users/binivazquez/UniWorkspace/computational_biology")

# Carga configuración
source("classwork/globalSetup.R")

cat("\n=== Probando DECIPHER ===\n")

# Carga secuencias
mexico <- readDNAStringSet("classwork/assets/secuencias_fasta/mexico.fasta")
wuhan <- readDNAStringSet("classwork/assets/wuhan.fasta")

cat("✅ Mexico:", length(mexico), "secuencias\n")
cat("✅ Wuhan:", length(wuhan), "secuencias\n")

# Combina
all_seqs <- c(mexico, wuhan)
names(all_seqs) <- c("Mexico", "Wuhan")

# Alinea
cat("\nAlineando... (puede tomar unos segundos)\n")
aligned <- AlignSeqs(all_seqs, verbose = TRUE)

cat("\n✅ ¡Alineamiento completo!\n")
cat("Largo:", width(aligned[1]), "pares de bases\n")

# Guarda
if (!dir.exists("outputs")) dir.create("outputs")
saveRDS(aligned, "outputs/test_alignment.rds")

cat("\n✅ Guardado en: outputs/test_alignment.rds\n")
cat("\n¡Todo funciona! Estás listo para análizar.\n")
```

Guarda con Cmd + S.

---

## Paso 10: Ejecuta Tu Primer Script

Abre `test_decipher.R` si no lo está.

Tienes opciones:

### Opción A: Ejecutar Todo

```
Cmd + Shift + Enter
```

Espera. Verás el progreso en la terminal de R.

### Opción B: Línea por Línea

Coloca el cursor en una línea y presiona:

```
Cmd + Enter
```

Útil para debuggear si algo explota.

---

## Atajos Útiles Para R en VS Code

| Acción | Tecla (Mac) |
|--------|------------|
| Ejecutar una línea | Cmd + Enter |
| Ejecutar selección | Cmd + Enter (con texto seleccionado) |
| Ejecutar script completo | Cmd + Shift + Enter |
| Crear terminal R | Cmd + Shift + P → "R: Create Terminal" |
| Ver documentación | Cmd + Shift + P → "R: Show Help" |
| Nuevo script R | Cmd + Shift + P → "R: Create R Script" |

---

## Si Algo No Funciona

### Problema: "R command not found"

```bash
# Verifica dónde está R
which R

# Si no está en /usr/local/bin/R, instala desde:
# https://cran.r-project.org/bin/macosx/
```

### Problema: "DECIPHER not found"

En R terminal:

```r
# Reinstala
if (!requireNamespace("BiocManager")) install.packages("BiocManager")
BiocManager::install("DECIPHER")

# Verifica
library(DECIPHER)
packageVersion("DECIPHER")
```

### Problema: Permiso denegado en Mac

```bash
# Arregla permisos
sudo chown -R $(whoami) /usr/local/bin/R
```

(Te pedirá tu contraseña)

### Problema: Alineamiento muy lento

```r
# En el script, usa menos iteraciones
aligned <- AlignSeqs(sequences, 
                    iterations = 1,      # En lugar de 2
                    processors = 2)      # Usa 2 núcleos
```

---

## Consejos Para Un Flujo De Trabajo Mejor

### Tip 1: Usa R Markdown

En VS Code, crea archivo `.Rmd`:

```markdown
---
title: "Mi Análisis"
output: html_document
---

## Cargar Datos

```{r setup}
source("classwork/globalSetup.R")
```

## Alineamiento

```{r align}
seqs <- readDNAStringSet("data.fasta")
aligned <- AlignSeqs(seqs)
```
```

Ejecuta con: Cmd + Shift + K

### Tip 2: Guarda Frecuentemente

```r
# Después de cualquier cálculo importante
saveRDS(aligned, "outputs/alignment_20260417.rds")
write.csv(results, "outputs/results.csv")
```

### Tip 3: Mantén Funciones Reutilizables

En `functions/helpers.R`:

```r
# Análisis rápido
quick_align <- function(fasta_file) {
  seqs <- readDNAStringSet(fasta_file)
  AlignSeqs(seqs, verbose = TRUE)
}

quick_tree <- function(aligned) {
  dist <- dist.hamming(aligned)
  nj(as.dist(dist))
}
```

Úsalas después:

```r
source("functions/helpers.R")
result <- quick_align("data.fasta")
```

### Tip 4: Control de Versión

```bash
# En terminal
cd /Users/binivazquez/UniWorkspace/computational_biology

# Inicializa Git
git init

# Agrega archivo
git add classwork/globalSetup.R

# Guarda
git commit -m "Setup inicial con DECIPHER"

# Sube si tienes repositorio remoto
git push origin main
```

---

## Mantenimiento Mensual (No Es Urgente)

```r
# Una vez al mes, actualiza paquetes
update.packages(ask = FALSE)
BiocManager::install()
```

```bash
# Limpia caché
rm -rf cache/*
```

---

## Referencia Rápida De Tu Setup

| Cosa | Valor |
|------|-------|
| R Location | `/usr/local/bin/R` |
| Proyecto | `/Users/binivazquez/UniWorkspace/computational_biology` |
| Comando para DECIPHER | `BiocManager::install("DECIPHER")` |
| R Studio | VS Code con extensión R |

---

## Checklist Final

Marca esto para saber que está todo listo:

- [ ] VS Code instalado
- [ ] Extensión R instalada
- [ ] R funciona en terminal
- [ ] DECIPHER instalado
- [ ] Carpetas creadas
- [ ] globalSetup.R en su lugar
- [ ] test_decipher.R funciona sin errores
- [ ] Puedes ejecutar código con Cmd + Enter

¡Si marqaste todo, estás LISTO! 🚀

---

## Siguiente: Empieza a Analizar

Ahora que todo funciona:

1. Lee [DECIPHER_VS_CODE_GUIDE_ES.md](DECIPHER_VS_CODE_GUIDE_ES.md) para aprender funciones
2. Usa [DECIPHER_CHEATSHEET_ES.md](DECIPHER_CHEATSHEET_ES.md) mientras programas
3. Consulta [DECIPHER_FAQ_ES.md](DECIPHER_FAQ_ES.md) si algo no funciona

---

**¡Felicidades!** 🧬

Tu ambiente de trabajo está configurado. Ahora solo necesitas:
- Curiosidad
- Paciencia
- Café

*Creado: Abril 2026*  
*Para: Estudiantes de Biología Computacional*
