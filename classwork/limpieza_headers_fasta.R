# ============================================================================
# LIMPIEZA DE HEADERS FASTA PARA ÁRBOLES FILOGENÉTICOS MEJORES
# ============================================================================
# Convierte headers largos en etiquetas cortas y legibles
# Formato: COUNTRY_VARIANT (ej: USA_Omicron, GHA_Beta, THA_Delta)
#
# Autor: Biniza Vázquez
# Elaborado con apoyo de Claude Haiku 4.5
# Abril 2026

library(Biostrings)

# ============================================================================
# CONFIGURACIÓN
# ============================================================================

input_fasta <- "classwork/assets/secuencias_fasta/all_arn_sequences.fasta"
output_fasta <- "classwork/assets/secuencias_fasta/all_arn_sequences_CLEAN.fasta"

# Mapeo de variantes conocidas (por código o palabras clave)
variant_mapping <- list(
    Alpha = c("B.1.1.7", "VOC-20DEC-01", "alpha"),
    Beta = c("B.1.351", "VOC-20MAY-01", "beta", "south african"),
    Gamma = c("P.1", "P1", "gamma", "brazilian"),
    Delta = c("B.1.617.2", "AY", "delta", "indian"),
    Lambda = c("C.37", "lambda", "peruvian"),
    Mu = c("B.1.621", "mu"),
    Omicron = c("B.1.1.529", "BA.", "omicron"),
    RaTG13 = c("RaTG13", "bat"),
    Pangolin = c("pangolin"),
    Wuhan = c("Wuhan-Hu-1", "wuhan", "original", "NC_045512")
)

# Códigos de país
country_codes <- list(
    USA = c("USA", "US"),
    GHA = c("GHA"),
    THA = c("THA"),
    HKG = c("HKG", "hong kong"),
    DEU = c("DEU", "germany"),
    FRA = c("FRA", "france"),
    ITA = c("ITA", "italy"),
    BEL = c("BEL", "belgium"),
    MYS = c("MYS", "malaysia"),
    MEX = c("MEX", "mexico"),
    IND = c("IND", "india"),
    LBY = c("LBY", "libya"),
    RUS = c("RUS", "russia"),
    CO = c("CO", "colombia")
)

# ============================================================================
# FUNCIONES AUXILIARES
# ============================================================================

extract_country <- function(header) {
    # Primero intenta encontrar patrón /human/[CODE]/
    # Este es el patrón más confiable en headers GenBank
    match <- regexpr("/human/([A-Z]{3})/", header, ignore.case = FALSE)
    if (match != -1) {
        # Extraer solo el código de país
        extracted <- regmatches(header, match)
        code <- gsub("/human/|/", "", extracted)
        return(code)
    }

    # Si no encuentra /human/CODE/, busca en la lista de synonyms
    for (country_name in names(country_codes)) {
        patterns <- country_codes[[country_name]]
        for (pattern in patterns) {
            if (grepl(pattern, header, ignore.case = TRUE)) {
                return(country_name)
            }
        }
    }

    # Como último recurso, busca cualquier código de 3 letras entre slashes
    match <- gregexpr("/([A-Z]{2,3})/", header)
    if (match[[1]][1] != -1) {
        extracted <- regmatches(header, match)[[1]]
        code <- gsub("/", "", extracted[1])
        return(code)
    }

    return("UNK")
}

extract_variant <- function(header) {
    # Primero intenta encontrar códigos Pango estrictamente
    pango_patterns <- list(
        Alpha = c("B\\.1\\.1\\.7"),
        Beta = c("B\\.1\\.351"),
        Gamma = c("P\\.1"),
        Delta = c("B\\.1\\.617\\.2", "AY\\."),
        Lambda = c("C\\.37"),
        Mu = c("B\\.1\\.621"),
        Omicron = c("B\\.1\\.1\\.529", "BA\\."),
        RaTG13 = c("RaTG13"),
        Pangolin = c("pangolin"),
        Wuhan = c("Wuhan-Hu-1", "NC_045512")
    )

    # Buscar códigos Pango primero (más específico)
    for (variant_name in names(pango_patterns)) {
        patterns <- pango_patterns[[variant_name]]
        for (pattern in patterns) {
            if (grepl(pattern, header, ignore.case = FALSE)) {
                return(variant_name)
            }
        }
    }

    # Si no encuentra código Pango, busca palabras clave (caso insensible)
    for (variant_name in names(variant_mapping)) {
        patterns <- variant_mapping[[variant_name]]
        for (pattern in patterns) {
            if (grepl(pattern, header, ignore.case = TRUE)) {
                return(variant_name)
            }
        }
    }

    return("WT") # Wild Type si no reconoce variante
}

create_clean_label <- function(header, sequence_number) {
    # Manejar casos especiales que no tienen país de humano
    if (grepl("RaTG13", header, ignore.case = FALSE)) {
        return("RaTG13_RaTG13")
    }
    if (grepl("Pangolin", header, ignore.case = TRUE)) {
        return("Pangolin_Pangolin")
    }
    if (grepl("Wuhan-Hu-1|NC_045512", header, ignore.case = FALSE)) {
        return("Wuhan_Original")
    }

    # Para el resto: COUNTRY_VARIANT
    country <- extract_country(header)
    variant <- extract_variant(header)

    # Formato: COUNTRY_VARIANT (ej: USA_Omicron)
    label <- paste0(country, "_", variant)

    return(label)
}

# ============================================================================
# LECTURA Y PROCESAMIENTO
# ============================================================================

cat("Leyendo archivo FASTA...\n")
seqs <- readDNAStringSet(input_fasta)

# Guardar headers ORIGINALES antes de modificar
original_headers <- names(seqs)

cat("Total de secuencias:", length(seqs), "\n")
cat("Headers originales (primeros 3):\n")
for (i in 1:min(3, length(seqs))) {
    cat("  ", original_headers[i], "\n")
}

# Crear nuevos nombres
cat("\nExtrayendo información y creando nuevos labels...\n")
new_names <- character(length(seqs))

for (i in seq_along(seqs)) {
    original_header <- original_headers[i]
    new_label <- create_clean_label(original_header, i)
    new_names[i] <- new_label

    if (i <= 5 || i > length(seqs) - 2) {
        cat("  ", i, ": ", original_header, " → ", new_label, "\n")
    } else if (i == 6) {
        cat("  ...\n")
    }
}

# Asignar nuevos nombres
names(seqs) <- new_names

# ============================================================================
# ESTADÍSTICAS
# ============================================================================

cat("\n--- RESUMEN DE ASIGNACIÓN ---\n")
cat("Total de secuencias procesadas:", length(seqs), "\n\n")

# Contar por país
countries <- sapply(strsplit(new_names, "_"), "[", 1)
cat("Distribución por país:\n")
country_counts <- table(countries)
for (country in sort(names(country_counts))) {
    cat("  ", country, ": ", country_counts[country], " secuencias\n")
}

# Contar por variante
variants <- sapply(strsplit(new_names, "_"), "[", 2)
cat("\nDistribución por variante:\n")
variant_counts <- table(variants)
for (variant in sort(names(variant_counts))) {
    cat("  ", variant, ": ", variant_counts[variant], " secuencias\n")
}

# ============================================================================
# GUARDAR NUEVO FASTA
# ============================================================================

cat("\nGuardando nuevo archivo FASTA...\n")
writeXStringSet(seqs, file = output_fasta, format = "fasta")

cat("✓ Archivo guardado en: ", output_fasta, "\n")
cat("  Puedes usar este archivo para árboles más legibles.\n\n")

# ============================================================================
# CREAR TABLA DE CORRESPONDENCIA
# ============================================================================

correspondence_table <- data.frame(
    original = original_headers,
    new_label = new_names,
    country = countries,
    variant = variants,
    stringsAsFactors = FALSE
)

output_csv <- "outputs/fasta_headers_correspondencia.csv"
if (!dir.exists("outputs")) {
    dir.create("outputs", recursive = TRUE)
}

write.csv(correspondence_table, file = output_csv, row.names = FALSE)
cat("✓ Tabla de correspondencia guardada: ", output_csv, "\n")

# ============================================================================
# PREVISUALIZACION
# ============================================================================

cat("\nPrimeros 10 headers transformados:\n")
print(head(correspondence_table, 10))

cat("\n", strrep("=", 70), "\n")
cat("PROCESO COMPLETADO\n")
cat(strrep("=", 70), "\n")
