# ============================================================================
# ANÁLISIS COMPARATIVO: 3 MÉTODOS CON FASTA
# ============================================================================
# Compara: Filogenético (NJ) vs Clustering (IdClusters) vs Matriz de Identidad
# Entrada: FASTA (más simple y directo que GenBank)
#
# Autor: Biniza Vázquez
# Elaborado con apoyo de Claude Haiku 4.5
# Abril 2026

library(Biostrings)
library(DECIPHER)
library(ape)

# Configuración de archivos
fasta_file <- "classwork/assets/secuencias_fasta/all_arn_sequences.fasta"
output_dir <- "outputs"

# Crear directorio si no existe
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# ============================================================================
# OPCIÓN 1: ÁRBOL FILOGENÉTICO (Neighbor-Joining + K80)
# ============================================================================
# Mejor para: Entender relaciones evolutivas entre variantes
# Interpreta: Tiempo evolutivo, ramas reflejan distancia genética

cat("\n========== OPCIÓN 1: ÁRBOL FILOGENÉTICO ==========\n")

# Preparar datos desde FASTA
cat("Leyendo FASTA...\n")
seqs1 <- readDNAStringSet(fasta_file)

# Alinear primero
cat("Alineando secuencias...\n")
aligned1 <- AlignSeqs(seqs1, verbose = FALSE)

# Calcular distancias
cat("Calculando distancias (K80)...\n")
d1 <- DistanceMatrix(aligned1, correction = "K80", verbose = FALSE)

# Construir árbol filogenético usando hclust para mejor control
cat("Construyendo árbol (NJ)...\n")
hc1 <- hclust(as.dist(d1), method = "average")
tree_phylo <- as.phylo(hc1)

# Visualizar
cat("Guardando visualización...\n")
png(file.path(output_dir, "arbol_filogenetico.png"), width = 1400, height = 700)
plot.phylo(tree_phylo,
    direction = "rightwards", cex = 0.7,
    main = "Árbol Filogenético (UPGMA + K80)\n28 variantes SARS-CoV-2"
)
dev.off()

cat("✓ Árbol filogenético guardado: ", file.path(output_dir, "arbol_filogenetico.png"), "\n")
cat("  Secuencias cargadas:", length(seqs1), "\n")
cat("  Método: UPGMA (ligamiento promedio)\n")
cat("  Modelo distancia: K80 (Kimura 2-Parameter)\n")

# ============================================================================
# OPCIÓN 2: CLUSTERING POR IDENTIDAD (IdClusters)
# ============================================================================
# Mejor para: Agrupar secuencias prácticamente idénticas
# Interpreta: Identidad de secuencia (%), no tiempo evolutivo

cat("\n========== OPCIÓN 2: CLUSTERING POR IDENTIDAD ==========\n")

# Preparar datos desde FASTA
cat("Leyendo FASTA...\n")
seqs2 <- readDNAStringSet(fasta_file)

# Alinear primero
cat("Alineando secuencias...\n")
aligned2 <- AlignSeqs(seqs2, verbose = FALSE)

# Calcular matriz de distancias
cat("Calculando distancias...\n")
d2 <- DistanceMatrix(aligned2, correction = "K80", verbose = FALSE)

# Clustering jerárquico basado en distancias
cat("Ejecutando clustering jerárquico...\n")
hc <- hclust(as.dist(d2), method = "complete")

# Cortar el árbol para obtener clusters
# cutoff de 0.05 distancia = ~95% identidad
clusters <- cutree(hc, h = 0.05)

cat("✓ Clustering completado\n")
cat("  Secuencias cargadas:", length(seqs2), "\n")
cat("  Umbral de distancia: 0.05\n")
cat("  Número de clusters:", max(clusters), "\n")

# Visualizar clusters
png(file.path(output_dir, "clusters_identidad.png"), width = 1000, height = 600)
barplot(table(clusters),
    main = "Distribución de Clusters por Identidad (95%)",
    xlab = "ID Cluster",
    ylab = "Número de secuencias"
)
dev.off()

cat("✓ Gráfico de clusters guardado: ", file.path(output_dir, "clusters_identidad.png"), "\n")

# Crear tabla de clusters
clusters_df <- data.frame(
    secuencia = names(aligned2),
    cluster = clusters
)

write.csv(clusters_df, file.path(output_dir, "clusters_asignacion.csv"), row.names = FALSE)
cat("✓ Tabla de asignaciones guardada: ", file.path(output_dir, "clusters_asignacion.csv"), "\n")

# ============================================================================
# OPCIÓN 3: MATRIZ DE IDENTIDAD
# ============================================================================
# Mejor para: Visualizar similitud entre secuencias
# Interpreta: Identidad global entre todas las secuencias

cat("\n========== OPCIÓN 3: MATRIZ DE IDENTIDAD ==========\n")

# Preparar datos desde FASTA
cat("Leyendo FASTA...\n")
seqs3 <- readDNAStringSet(fasta_file)

# Alinear primero (requerido para matriz)
cat("Alineando secuencias...\n")
aligned3 <- AlignSeqs(seqs3, verbose = FALSE)

# Crear matriz de identidad (basada en distancias)
cat("Calculando matriz de identidad...\n")
d3 <- DistanceMatrix(aligned3, correction = "K80", verbose = FALSE)
# Convertir distancia a identidad: identidad = (1 - distancia) * 100
identity_matrix <- (1 - d3) * 100

# Visualizar matriz de identidad
cat("Guardando heatmap...\n")
png(file.path(output_dir, "matriz_identidad.png"), width = 1000, height = 900)
heatmap(identity_matrix,
    main = "Matriz de Identidad de Secuencias (%)\n28 variantes SARS-CoV-2",
    col = colorRampPalette(c("white", "yellow", "orange", "red"))(50)
)
dev.off()

cat("✓ Matriz de identidad guardada: ", file.path(output_dir, "matriz_identidad.png"), "\n")

# ============================================================================
# RESUMEN COMPARATIVO
# ============================================================================

cat("\n", strrep("=", 70), "\n")
cat("RESUMEN COMPARATIVO\n")
cat(strrep("=", 70), "\n\n")

cat("┌─ OPCIÓN 1: ÁRBOL FILOGENÉTICO ─────────────────────────────────┐\n")
cat("│ Método:       UPGMA (ligamiento promedio)                      │\n")
cat("│ Modelo:       K80 (Kimura 2-Parameter)                         │\n")
cat("│ Resultado:    Árbol evolutivo (ramas = distancia genética)    │\n")
cat("│ Uso:          Entender relaciones evolutivas y divergencia    │\n")
cat("│ Archivo:      arbol_filogenetico.png                           │\n")
cat("└─────────────────────────────────────────────────────────────────┘\n\n")

cat("┌─ OPCIÓN 2: CLUSTERING JERÁRQUICO ────────────────────────────┐\n")
cat("│ Método:       hclust (complete linkage)                     │\n")
cat("│ Umbral:       Distancia 0.05 (≈95% identidad)             │\n")
cat("│ Resultado:    Grupos de secuencias por similitud          │\n")
cat("│ Uso:          Agrupar variantes cercanas                  │\n")
cat("│ Archivos:     clusters_identidad.png, clusters_asignacion.csv │\n")
cat("└─────────────────────────────────────────────────────────────┘\n\n")

cat("┌─ OPCIÓN 3: MATRIZ DE IDENTIDAD ─────────────────────────────┐\n")
cat("│ Método:       DistanceMatrix K80 → Identidad (%)            │\n")
cat("│ Resultado:    Matriz de similitud entre todas las seqs     │\n")
cat("│ Uso:          Visualizar relaciones de identidad global    │\n")
cat("│ Archivo:      matriz_identidad.png                         │\n")
cat("└──────────────────────────────────────────────────────────────┘\n\n")

cat("ENTRADA DE DATOS:\n")
cat("  Archivo FASTA: ", fasta_file, "\n")
cat("  Secuencias totales: 28 variantes SARS-CoV-2\n\n")

cat("ARCHIVOS GENERADOS EN", output_dir, ":\n")
cat("  ✓ arbol_filogenetico.png\n")
cat("  ✓ clusters_identidad.png\n")
cat("  ✓ clusters_asignacion.csv\n")
cat("  ✓ matriz_identidad.png\n\n")

cat("RECOMENDACIÓN PARA TU PROYECTO:\n")
cat("  • Para análisis filogeográfico → OPCIÓN 1 (filogenético)\n")
cat("  • Para agrupar variantes cercanas → OPCIÓN 2 (clustering)\n")
cat("  • Para comparar identidad global → OPCIÓN 3 (matriz)\n")
cat("  • Para presentación integrada → TODAS LAS 3\n\n")

cat(strrep("=", 70), "\n")
cat("Análisis completado.\n")
cat(strrep("=", 70), "\n")
