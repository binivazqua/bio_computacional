setwd("/Users/binivazquez/UniWorkspace/computational_biology")
source("classwork/globalSetup.R")

fasta_files <- project_fasta_files

sequences <- read_project_fastas()

## LONGITUDES
cat("Sequence lengths:\n")
for (i in seq_along(sequences)) {
    cat(names(sequences)[i], ":", width(sequences)[i], "\n")
}

## GRAFICA
compara_all <- function(seq_set) {
    base_levels <- c("a", "c", "g", "t")
    plot_data <- do.call(
        rbind,
        lapply(seq_along(seq_set), function(i) {
            seq_chars <- tolower(strsplit(as.character(seq_set[[i]]), "")[[1]])
            seq_table <- table(factor(seq_chars, levels = base_levels))

            data.frame(
                sequence = names(seq_set)[i],
                base = c("A", "C", "G", "T"),
                count = as.numeric(seq_table)
            )
        })
    )

    composition_plot <- ggplot(
        plot_data,
        aes(x = base, y = count, fill = base)
    ) +
        geom_col() +
        facet_wrap(~sequence, scales = "free_y") +
        labs(
            title = "Composicion de bases por secuencia",
            x = "Base",
            y = "Conteo"
        ) +
        theme_minimal() +
        theme(legend.position = "none")
    png("base_composition.png", width = 800, height = 600)
    print(composition_plot)
}

# FUNC EN ACTION
quartz()
compara_all(sequences)
png("comparacion_fasta.png")


# porcentajes
gc_table <- data.frame(
    variant = names(sequences),
    gc_percent = rowSums(
        letterFrequency(sequences, letters = c("G", "C"), as.prob = TRUE)
    ) * 100
)

gc_table
png("gc_table.png")

# COMPLEMENTOS

sample_names <- names(sequences)

seq_vectors <- setNames(lapply(seqs, as.vector), sample_names)

reverse_complements <- data.frame(
    Sample = names(seq_vectors),
    ReverseComplement = vapply(
        seq_vectors,
        function(x) paste(reverse_complement(x), collapse = ""),
        character(1)
    )
)

head(reverse_complements)
