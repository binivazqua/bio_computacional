# Análisis de secuencias GenBank con DECIPHER
# Basado en: guia_que_si_funciona.md
# Autor: Biniza Vázquez

library(DBI)
library(RSQLite)
library(DECIPHER)
library(ape)

# Definir rutas
gb_file <- "classwork/assets/secuencias_gb/all_sequences.gen"

# Verificar que el archivo existe
if (!file.exists(gb_file)) {
  stop("Archivo no encontrado: ", gb_file)
}

cat("\n=== Cargar secuencias GenBank ===\n")

# Crear base de datos temporal
db_path <- "outputs/temp_sequences.sqlite"
dbConn <- dbConnect(RSQLite::SQLite(), db_path)

# Cargar las secuencias desde GenBank
Seqs2DB(
  gb_file,
  type = "GenBank",
  dbFile = db_path,
  identifier = "coronavirus"
)

# Recuperar las secuencias
dna <- SearchDB(dbConn)

cat("Secuencias cargadas:", length(dna), "\n")
cat("Longitud de secuencias:", width(dna), "\n\n")

# === ALINEAMIENTO ===
cat("=== Alineando secuencias ===\n")
aligned <- AlignSeqs(dna, verbose = TRUE)

cat("Alineamiento completado\n")
cat("Longitud del alineamiento:", width(aligned[1]), "\n\n")

# === MATRIZ DE DISTANCIAS ===
cat("=== Calculando matriz de distancias ===\n")

d <- DistanceMatrix(aligned, correction = "Jukes-Cantor")

cat("Matriz de distancias calculada\n\n")

# === ÁRBOL FILOGENÉTICO ===
cat("=== Construyendo árbol filogenético ===\n")

tree <- TreeLine(myDistMatrix = d)

cat("Árbol construido\n")
cat("Clase del árbol:", class(tree), "\n\n")

# === VISUALIZACIÓN ===
cat("=== Generando gráfico ===\n")

# Convertir a clase phylo (compatible con ape)
phylo_tree <- as.phylo(tree)

# Crear archivo PNG de alta calidad
png("outputs/arbol_filogenético_gb.png", width = 1200, height = 800, res = 100)
plot(phylo_tree, main = "Árbol Filogenético - SARS-CoV-2 (GenBank)")
axisPhylo()
dev.off()

cat("Árbol guardado en: outputs/arbol_filogenético_gb.png\n\n")

# === RESUMEN ===
cat("=== Resumen del análisis ===\n")
cat("Secuencias procesadas:", length(dna), "\n")
cat("Tipo de árbol:", class(tree), "\n")
cat("Nodo raíz:", tree$root.node, "\n")
cat("Total de nodos:", nrow(tree$edge), "\n\n")

# Guardar resultado
saveRDS(list(
  dna = dna,
  aligned = aligned,
  distance_matrix = d,
  tree = tree,
  phylo_tree = phylo_tree
), "outputs/analisis_gb_completo.rds")

cat("Análisis completo guardado en: outputs/analisis_gb_completo.rds\n")

# Cerrar conexión
dbDisconnect(dbConn)

cat("\n=== Análisis finalizado ===\n")
