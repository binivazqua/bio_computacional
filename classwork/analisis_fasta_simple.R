# ============================================================================
# ANÁLISIS FILOGENÉTICO CON FASTA (SIMPLE Y DIRECTO)
# ============================================================================
# Alternativa simple a GenBank: trabaja directamente con FASTA
# No requiere bases de datos SQLite
#
# Autor: Biniza Vázquez
# Elaborado con apoyo de Claude Haiku 4.5
# Abril 2026

library(Biostrings)
library(DECIPHER)
library(ape)

# ============================================================================
# CONFIGURACIÓN
# ============================================================================

fasta_file <- "classwork/assets/secuencias_fasta/all_arn_sequences_CLEAN.fasta"
output_dir <- "outputs"

if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# ============================================================================
# FLUJO SIMPLE CON FASTA
# ============================================================================

cat("\n", strrep("=", 70), "\n")
cat("ANÁLISIS FILOGENÉTICO CON FASTA\n")
cat(strrep("=", 70), "\n\n")

# Paso 1: Leer FASTA
cat("Paso 1: Leyendo archivo FASTA...\n")
seqs <- readDNAStringSet(fasta_file)
cat("✓ Secuencias cargadas:", length(seqs), "\n")
cat("  Ejemplos de labels:", paste(head(names(seqs), 3), collapse = ", "), "\n\n")

# Paso 2: Alinear secuencias
cat("Paso 2: Alineando secuencias (esto toma tiempo)...\n")
aligned <- AlignSeqs(seqs, verbose = FALSE)
cat("✓ Alineamiento completado\n")
cat("  Largo de alineamiento:", width(aligned)[1], "bp\n\n")

# Paso 3: Calcular matriz de distancias
cat("Paso 3: Calculando distancias (modelo K80)...\n")
d <- DistanceMatrix(aligned, correction = "K80", verbose = FALSE)
cat("✓ Matriz de distancias calculada\n")
cat("  Dimensiones:", nrow(d), "x", ncol(d), "\n\n")

# Construir árbol filogenético usando hclust
cat("Paso 4: Construyendo árbol filogenético (UPGMA)...\n")
hc <- hclust(as.dist(d), method = "average")
tree <- as.phylo(hc)
cat("✓ Árbol construido\n\n")

# Paso 5: Visualizar y guardar
cat("Paso 5: Generando visualizaciones...\n")

# Gráfico 1: Árbol simple horizontal
png(file.path(output_dir, "filogenetico_fasta_simple.png"), width = 1400, height = 700)
plot.phylo(tree,
    direction = "rightwards", cex = 0.7,
    main = "Árbol Filogenético FASTA (UPGMA + K80)\n28 variantes SARS-CoV-2"
)
dev.off()
cat("✓ Árbol simple: ", file.path(output_dir, "filogenetico_fasta_simple.png"), "\n")

# Gráfico 2: Con ggtree (si está disponible)
tryCatch(
    {
        library(ggtree)

        # Convertir a phylo
        tree_phylo <- as.phylo(as.hclust(tree))

        p <- ggtree(tree_phylo) +
            geom_tiplab(size = 3, align = TRUE) +
            geom_treescale(x = 0, y = -2, width = 0.05) +
            theme_tree2() +
            ggtitle("Árbol Filogenético FASTA (ggtree)")

        ggsave(file.path(output_dir, "filogenetico_fasta_ggtree.png"),
            p,
            width = 16, height = 8, dpi = 150
        )
        cat("✓ Árbol mejorado: ", file.path(output_dir, "filogenetico_fasta_ggtree.png"), "\n")
    },
    error = function(e) {
        cat("⚠️  ggtree no disponible, saltando gráfico mejorado\n")
    }
)

# Gráfico 3: Matriz de distancias
png(file.path(output_dir, "matriz_distancias_fasta.png"), width = 1000, height = 900)
heatmap(d,
    main = "Matriz de Distancias K80\n28 variantes SARS-CoV-2",
    col = colorRampPalette(c("white", "yellow", "orange", "red"))(50),
    cexRow = 0.6,
    cexCol = 0.6
)
dev.off()
cat("✓ Matriz de distancias: ", file.path(output_dir, "matriz_distancias_fasta.png"), "\n\n")

# ============================================================================
# ESTADÍSTICAS
# ============================================================================

cat("ESTADÍSTICAS DE DISTANCIA:\n")
cat("  Distancia mínima:", min(d, na.rm = TRUE), "\n")
cat("  Distancia máxima:", max(d, na.rm = TRUE), "\n")
cat("  Distancia media:", mean(d[upper.tri(d)]), "\n\n")

# ============================================================================
# INFORMACIÓN DE SECUENCIAS
# ============================================================================

cat("INFORMACIÓN DE SECUENCIAS:\n")
info_table <- data.frame(
    Label = names(seqs),
    Longitud_bp = width(seqs),
    Tipo = "DNA"
)

cat("\nPrimeras 10 secuencias:\n")
print(head(info_table, 10))

write.csv(info_table, file.path(output_dir, "secuencias_info_fasta.csv"), row.names = FALSE)
cat("\n✓ Tabla de secuencias guardada:", file.path(output_dir, "secuencias_info_fasta.csv"), "\n")

# ============================================================================
# RESUMEN FINAL
# ============================================================================

cat("\n", strrep("=", 70), "\n")
cat("ANÁLISIS COMPLETADO\n")
cat(strrep("=", 70), "\n\n")

cat("ARCHIVOS GENERADOS EN", output_dir, ":\n")
cat("  ✓ filogenetico_fasta_simple.png\n")
cat("  ✓ filogenetico_fasta_ggtree.png (si ggtree disponible)\n")
cat("  ✓ matriz_distancias_fasta.png\n")
cat("  ✓ secuencias_info_fasta.csv\n\n")

cat("PARÁMETROS UTILIZADOS:\n")
cat("  • Modelo de distancia: K80 (Kimura 2-Parameter)\n")
cat("  • Método de árbol: NJ (Neighbor-Joining)\n")
cat("  • Umbral de cutoff: 0.05\n")
cat("  • Formato entrada: FASTA\n\n")

cat("PRÓXIMOS PASOS:\n")
cat("  1. Abre filogenetico_fasta_simple.png para ver el árbol básico\n")
cat("  2. Si ggtree está disponible, usa filogenetico_fasta_ggtree.png para mejor visualización\n")
cat("  3. Consulta matriz_distancias_fasta.png para ver similitud entre secuencias\n")
cat("  4. Usa secuencias_info_fasta.csv para datos de las secuencias\n\n")

cat(strrep("=", 70), "\n")
