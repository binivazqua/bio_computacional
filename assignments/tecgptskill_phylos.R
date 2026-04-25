# Tablero filogenetico interactivo inspirado en ideas de exploracion
# comunes en phylo.io (navegacion/ajuste visual) e iTOL
# (layouts, anotaciones y exportacion).

suppressPackageStartupMessages({
  library(DECIPHER)
  library(Biostrings)
  library(ape)
  library(shiny)
  library(ggplot2)
  library(ggtree)
})

resolve_default_fasta <- function() {
  candidates <- c(
    "classwork/assets/secuencias_fasta/act_gptskill/sequences_gptskill.fasta",
    "../classwork/assets/secuencias_fasta/act_gptskill/sequences_gptskill.fasta"
  )

  for (path in candidates) {
    if (file.exists(path)) {
      return(normalizePath(path, winslash = "/", mustWork = TRUE))
    }
  }

  stop(
    "No pude encontrar el FASTA por defecto. ",
    "Ajusta la ruta en resolve_default_fasta() o carga un archivo desde la app."
  )
}

clean_label <- function(x) {
  gsub("\\s+", " ", trimws(x))
}

infer_host_group <- function(labels) {
  lowered <- tolower(labels)

  ifelse(
    grepl("bat", lowered),
    "Bat",
    ifelse(
      grepl("pangolin", lowered),
      "Pangolin",
      ifelse(
        grepl("sars|severe acute respiratory syndrome|human|wuhan", lowered),
        "SARS-CoV-2 / Human",
        "Other"
      )
    )
  )
}

build_sequence_summary <- function(seqs) {
  gc_matrix <- letterFrequency(seqs, letters = c("G", "C"), as.prob = TRUE)

  data.frame(
    label = clean_label(names(seqs)),
    length_bp = width(seqs),
    gc_content = round(rowSums(gc_matrix) * 100, 2),
    group = infer_host_group(names(seqs)),
    stringsAsFactors = FALSE
  )
}

build_distance_long <- function(distance_matrix) {
  dm <- as.matrix(distance_matrix)
  idx <- which(upper.tri(dm), arr.ind = TRUE)

  data.frame(
    seq_1 = rownames(dm)[idx[, 1]],
    seq_2 = colnames(dm)[idx[, 2]],
    distance = dm[idx],
    stringsAsFactors = FALSE
  )
}

build_phylo_tree <- function(distance_matrix, method = c("NJ", "UPGMA")) {
  method <- match.arg(method)
  dist_object <- as.dist(distance_matrix)

  if (method == "NJ") {
    return(ape::nj(dist_object))
  }

  hc <- hclust(dist_object, method = "average")
  as.phylo(hc)
}

build_phylo_workspace <- function(file_path, correction = "K80", method = "NJ") {
  seqs <- readDNAStringSet(file_path)
  names(seqs) <- clean_label(names(seqs))

  alignment <- AlignSeqs(seqs, verbose = FALSE)
  distance_matrix <- DistanceMatrix(
    alignment,
    correction = correction,
    verbose = FALSE
  )

  tree <- build_phylo_tree(distance_matrix, method = method)
  summary_df <- build_sequence_summary(seqs)

  list(
    file_path = file_path,
    seqs = seqs,
    alignment = alignment,
    distance_matrix = distance_matrix,
    tree = tree,
    summary = summary_df,
    pairwise = build_distance_long(distance_matrix)
  )
}

make_tree_plot <- function(tree, summary_df, layout, branch_width, label_size,
                           tip_size, label_offset, show_labels, ladderize_tree,
                           search_term) {
  tree_to_plot <- if (ladderize_tree) ape::ladderize(tree) else tree

  p <- ggtree(
    tree_to_plot,
    layout = layout,
    linewidth = branch_width,
    color = "#5b6770"
  )

  p$data$group <- summary_df$group[match(p$data$label, summary_df$label)]
  p$data$highlight <- FALSE

  if (nzchar(search_term)) {
    p$data$highlight <- grepl(search_term, p$data$label, ignore.case = TRUE)
  }

  p <- p +
    geom_tippoint(
      data = subset(p$data, isTip),
      aes(x = x, y = y, color = group),
      size = tip_size,
      alpha = 0.9,
      inherit.aes = FALSE
    ) +
    scale_color_manual(
      values = c(
        "Bat" = "#2b8a3e",
        "Pangolin" = "#c77d1f",
        "SARS-CoV-2 / Human" = "#c92a2a",
        "Other" = "#495057"
      ),
      drop = FALSE
    )

  if (show_labels) {
    p <- p +
      geom_tiplab(
        size = label_size,
        color = "#495057",
        align = layout == "rectangular",
        offset = label_offset
      )

    if (any(p$data$highlight, na.rm = TRUE)) {
      p <- p +
        geom_tiplab(
          data = subset(p$data, isTip & highlight),
          size = label_size + 0.5,
          color = "#d9480f",
          fontface = "bold",
          align = layout == "rectangular",
          offset = label_offset
        )
    }
  }

  if (layout == "rectangular") {
    p <- p + theme_tree2() + geom_treescale(fontsize = 3, linesize = 0.6)
  } else {
    p <- p + theme_tree()
  }

  p +
    ggtitle("Tablero Filogenetico Interactivo") +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      legend.title = element_blank(),
      legend.position = "bottom"
    )
}

make_distance_heatmap <- function(distance_matrix) {
  dm <- as.matrix(distance_matrix)
  heat_df <- as.data.frame(as.table(dm), stringsAsFactors = FALSE)
  colnames(heat_df) <- c("seq_x", "seq_y", "distance")

  ggplot(heat_df, aes(seq_x, seq_y, fill = distance)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradient(
      low = "#f1f3f5",
      high = "#c92a2a"
    ) +
    labs(
      title = "Mapa de calor de distancias",
      x = NULL,
      y = NULL,
      fill = "Distancia"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 9),
      panel.grid = element_blank()
    )
}

make_sequence_plot <- function(summary_df) {
  ggplot(summary_df, aes(reorder(label, length_bp), length_bp, fill = group)) +
    geom_col(width = 0.7) +
    geom_text(
      aes(label = paste0(length_bp, " bp")),
      hjust = -0.1,
      size = 3.2
    ) +
    coord_flip() +
    scale_fill_manual(
      values = c(
        "Bat" = "#2b8a3e",
        "Pangolin" = "#c77d1f",
        "SARS-CoV-2 / Human" = "#c92a2a",
        "Other" = "#495057"
      ),
      drop = FALSE
    ) +
    labs(
      title = "Longitud de las secuencias",
      x = NULL,
      y = "Pares de bases",
      fill = "Grupo"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    )
}

make_gc_plot <- function(summary_df) {
  ggplot(summary_df, aes(label, gc_content, color = group)) +
    geom_point(size = 3) +
    geom_segment(
      aes(x = label, xend = label, y = 0, yend = gc_content),
      linewidth = 0.8,
      alpha = 0.5
    ) +
    scale_color_manual(
      values = c(
        "Bat" = "#2b8a3e",
        "Pangolin" = "#c77d1f",
        "SARS-CoV-2 / Human" = "#c92a2a",
        "Other" = "#495057"
      ),
      drop = FALSE
    ) +
    labs(
      title = "Contenido GC por secuencia",
      x = NULL,
      y = "GC (%)",
      color = "Grupo"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 35, hjust = 1),
      legend.position = "bottom"
    )
}

run_phylo_dashboard <- function(default_path = resolve_default_fasta()) {
  ui <- fluidPage(
    titlePanel("Explorador filogenetico estilo phylo.io / iTOL"),
    sidebarLayout(
      sidebarPanel(
        fileInput(
          "fasta_file",
          "Cargar archivo FASTA",
          accept = c(".fasta", ".fa", ".fna")
        ),
        selectInput(
          "correction",
          "Modelo de distancia",
          choices = c("K80", "JC69", "raw"),
          selected = "K80"
        ),
        selectInput(
          "tree_method",
          "Metodo del arbol",
          choices = c(
            "Neighbor-Joining" = "NJ",
            "UPGMA (average linkage)" = "UPGMA"
          ),
          selected = "NJ"
        ),
        selectInput(
          "layout",
          "Layout del arbol",
          choices = c(
            "Rectangular" = "rectangular",
            "Circular" = "circular",
            "Fan" = "fan",
            "Unrooted" = "unrooted"
          ),
          selected = "rectangular"
        ),
        textInput(
          "search_term",
          "Buscar y resaltar etiquetas",
          placeholder = "Ej. bat, pangolin, Wuhan"
        ),
        sliderInput(
          "branch_width",
          "Grosor de ramas",
          min = 0.2,
          max = 2.5,
          value = 0.8,
          step = 0.1
        ),
        sliderInput(
          "label_size",
          "Tamano de etiquetas",
          min = 2,
          max = 8,
          value = 3.5,
          step = 0.2
        ),
        sliderInput(
          "tip_size",
          "Tamano de puntas",
          min = 1,
          max = 6,
          value = 2.8,
          step = 0.2
        ),
        sliderInput(
          "label_offset",
          "Separacion de etiquetas",
          min = 0,
          max = 0.1,
          value = 0.01,
          step = 0.005
        ),
        sliderInput(
          "plot_height",
          "Altura del panel del arbol",
          min = 500,
          max = 1200,
          value = 760,
          step = 20
        ),
        checkboxInput("show_labels", "Mostrar etiquetas", value = TRUE),
        checkboxInput("ladderize_tree", "Ordenar ramas (ladderize)", value = TRUE),
        downloadButton("download_newick", "Exportar Newick"),
        downloadButton("download_png", "Exportar PNG")
      ),
      mainPanel(
        fluidRow(
          column(4, wellPanel(h4("Archivo"), textOutput("active_file"))),
          column(2, wellPanel(h4("Secuencias"), textOutput("n_sequences"))),
          column(3, wellPanel(h4("Longitud media"), textOutput("avg_length"))),
          column(3, wellPanel(h4("Distancia media"), textOutput("avg_distance")))
        ),
        tabsetPanel(
          tabPanel("Arbol", plotOutput("tree_plot")),
          tabPanel("Distancias", plotOutput("distance_heatmap", height = 650)),
          tabPanel("Longitudes", plotOutput("length_plot", height = 650)),
          tabPanel("GC", plotOutput("gc_plot", height = 650)),
          tabPanel("Tabla", tableOutput("summary_table"))
        )
      )
    )
  )

  server <- function(input, output, session) {
    active_path <- reactive({
      if (!is.null(input$fasta_file)) {
        return(input$fasta_file$datapath)
      }
      default_path
    })

    workspace <- reactive({
      build_phylo_workspace(
        file_path = active_path(),
        correction = input$correction,
        method = input$tree_method
      )
    })

    output$active_file <- renderText({
      basename(workspace()$file_path)
    })

    output$n_sequences <- renderText({
      nrow(workspace()$summary)
    })

    output$avg_length <- renderText({
      paste0(round(mean(workspace()$summary$length_bp), 1), " bp")
    })

    output$avg_distance <- renderText({
      round(mean(workspace()$pairwise$distance), 5)
    })

    output$tree_plot <- renderPlot(
      {
        make_tree_plot(
          tree = workspace()$tree,
          summary_df = workspace()$summary,
          layout = input$layout,
          branch_width = input$branch_width,
          label_size = input$label_size,
          tip_size = input$tip_size,
          label_offset = input$label_offset,
          show_labels = input$show_labels,
          ladderize_tree = input$ladderize_tree,
          search_term = input$search_term
        )
      },
      height = function() input$plot_height,
      res = 110
    )

    output$distance_heatmap <- renderPlot({
      make_distance_heatmap(workspace()$distance_matrix)
    }, res = 110)

    output$length_plot <- renderPlot({
      make_sequence_plot(workspace()$summary)
    }, res = 110)

    output$gc_plot <- renderPlot({
      make_gc_plot(workspace()$summary)
    }, res = 110)

    output$summary_table <- renderTable({
      workspace()$summary
    }, striped = TRUE, bordered = TRUE, spacing = "s")

    output$download_newick <- downloadHandler(
      filename = function() {
        paste0("tree_", input$tree_method, "_", input$correction, ".nwk")
      },
      content = function(file) {
        write.tree(workspace()$tree, file = file)
      }
    )

    output$download_png <- downloadHandler(
      filename = function() {
        paste0("tree_", input$tree_method, "_", input$correction, ".png")
      },
      content = function(file) {
        png(file, width = 1800, height = 1200, res = 180)
        print(
          make_tree_plot(
            tree = workspace()$tree,
            summary_df = workspace()$summary,
            layout = input$layout,
            branch_width = input$branch_width,
            label_size = input$label_size,
            tip_size = input$tip_size,
            label_offset = input$label_offset,
            show_labels = input$show_labels,
            ladderize_tree = input$ladderize_tree,
            search_term = input$search_term
          )
        )
        dev.off()
      }
    )
  }

  shinyApp(ui = ui, server = server)
}

if (interactive()) {
  runApp(run_phylo_dashboard(), launch.browser = TRUE)
} else {
  workspace <- build_phylo_workspace(resolve_default_fasta())
  cat(
    "Tablero listo.\n",
    "Archivo: ", basename(workspace$file_path), "\n",
    "Secuencias: ", nrow(workspace$summary), "\n",
    "Metodo por defecto: NJ | Modelo: K80\n",
    "Para abrir la interfaz: source('assignments/tecgptskill_phylos.R') en una sesion interactiva de R.\n",
    sep = ""
  )
}
