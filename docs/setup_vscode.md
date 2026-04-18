Autoría: Biniza Vázquez y OpenAI Codex  
Abril de 2026

# Configuración de VS Code para este proyecto

Esta guía deja configurado VS Code para trabajar con R y DECIPHER sin confundir los flujos FASTA y GenBank que usa este repositorio.

## Requisitos

- VS Code instalado
- R disponible desde terminal
- Extensión `R` de `REditorSupport.r`
- Paquetes de R necesarios para el flujo que vayas a ejecutar

## Verificación rápida de R en la terminal

```bash
which R
R --version
Rscript --version
```

## Configuración recomendada de VS Code

En `settings.json`:

```json
{
  "[r]": {
    "editor.defaultFormatter": "REditorSupport.r",
    "editor.formatOnSave": true,
    "editor.wordWrap": "on"
  },
  "r.sessionWatcher": true,
  "r.plot.useHttpgd": true
}
```

Si necesitas fijar manualmente la ruta de R en macOS, agrega `r.rterm.mac` con la ruta que devuelva `which R`.

## Paquetes de R usados en el repositorio

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c("DECIPHER", "Biostrings"))
install.packages(c("DBI", "RSQLite", "ape", "seqinr"))
```

## Archivo correcto según el flujo

### Para flujo FASTA

```text
classwork/assets/secuencias_fasta/all_arn_sequences.fasta
```

Se lee así:

```r
seqs <- Biostrings::readDNAStringSet(
  "classwork/assets/secuencias_fasta/all_arn_sequences.fasta"
)
```

### Para flujo GenBank

```text
classwork/assets/secuencias_gb/all_sequences.gen
```

Se procesa así:

```r
dbConn <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")

DECIPHER::Seqs2DB(
  "classwork/assets/secuencias_gb/all_sequences.gen",
  type = "GenBank",
  dbConn = dbConn,
  identifier = "coronavirus"
)
```

## Error de referencia para depuración

Si en la terminal de R aparece esto:

```text
Error in .Call2("fasta_index", filexp_list, nrec, skip, seek.first.rec, :
  reading FASTA file /Users/binivazquez/UniWorkspace/computational_biology/classwork/assets/secuencias_gb/all_sequences.gen: ">" expected at beginning of line 1
```

significa que se intentó leer un GenBank como FASTA.

## Secuencia mínima de diagnóstico

En la terminal integrada de VS Code:

```bash
sed -n '1p' classwork/assets/secuencias_gb/all_sequences.gen
sed -n '1p' classwork/assets/secuencias_fasta/all_arn_sequences.fasta
```

```r
file.exists("classwork/assets/secuencias_gb/all_sequences.gen")
file.exists("classwork/assets/secuencias_fasta/all_arn_sequences.fasta")
```

## Criterio de trabajo

- Si la primera línea es `LOCUS`, usa el flujo GenBank.
- Si la primera línea es `>`, usa el flujo FASTA.
- No cambies de función hasta comprobar el contenido del archivo.

## Documento principal

Para el flujo validado del proyecto, consulta primero `docs/guia_que_si_funciona.md`.
