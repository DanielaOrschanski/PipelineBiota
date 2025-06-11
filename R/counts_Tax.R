#' @title Counts Taxonomy
#' @import readxl
#' @import data.table
#' @import openxlsx
#' @description kvjdfnkvjdf
#' @param patient_dir is a dataframe fskbfjsf
#' @param source kjsnfn
#' @return cants_original dataframe that contains how many subjects from the same taxonomic level were found.
#' @examples counts_taxonomy <- counts_Tax(otus_table_kraken04)
#' @export
counts_Tax <- function(patients_dir, source, de_host, conEukaryota) {

  if(conEukaryota == TRUE) {
    conE <- "conEUKARYOTA"
  } else {
    conE <- ""
  }

  list_dirs <- list.dirs(patients_dir, recursive = FALSE, full.names = TRUE)
  #list_dirs <- list_dirs[1:50]
  #patient_dir <- list_dirs[1]

  if(length(nchar(list_dirs)) == 0 |  length(nchar(list_dirs[endsWith(list_dirs, "/trimmed")])) !=0 ) {
    list_dirs <- patients_dir
  }

  counts_Tax <- data.frame("ID" = NA,
                               "Domain" = NA,
                               "Kingdom" = NA,
                               "Phylum"  = NA,
                               "Class" = NA,
                               "Order" = NA,
                               "Family" = NA,
                               "Genus" = NA,
                               "Species" = NA,
                               "SubSpecies" = NA,
                               "Unclassified" = NA,
                           "Source" = source)
  i=1
  #patient_dir <- list_dirs[1]
  for (patient_dir in list_dirs) {

    id <- basename(patient_dir)

    if(de_host  == "Bowtie") {
      de_host_file <- "Bo"
      project_id <- "438230795"
      carpeta_bs <- "dehostBowtie"
      nombre_p <- sprintf("%sDHBo", id)

    } else if( de_host == "BWA") {
      de_host_file <- "bwa"
      project_id <- "436873438"
      carpeta_bs <- "dehostBWA"
      nombre_p <- sprintf("%sDHbwa", id)
    } else if(de_host == "RSubread") {
      de_host_file <- "Rs"
      project_id = "436780349"
      carpeta_bs ="dehostRsubread"
      nombre_p <- sprintf("%sDHRs", id)
    } else if(de_host == "") {
      de_host_file <- "sin"
      project_id <- "436799364"
      carpeta_bs <- "MuestrasTrimmed"
      nombre_p <- sprintf("%sT", id)
    } else if(de_host == "crudas") {
      de_host_file <- "cruda"
      project_id <- "438331896"
      carpeta_bs <- "MuestrasCrudas"
      nombre_p <-  paste0(id, "C", sep="")
    } else if(de_host == "Biota") {
      de_host_file <- "trimmed"

    } else if(de_host == "sinDH_PD") {
      de_host_file <- "sinDH_PD"

    } else {
      stop("de_host must be Bowtie, BWA, RSubread, sinDH_PD or empty string")
    }

    print(patient_dir)
    patient_dir <- paste(patient_dir, "/trimmed", sep="")
    print(patient_dir)

    if (source == "KRAKEN") {
      if(de_host == "Biota") {
        ReportSequences <- fread(paste(patient_dir, "/report.sequences", sep=""),  header = FALSE,sep = "\t")
        tabla_otus <-  read_excel(paste(patient_dir, "/TablaOTUS_", id, "_trimmed_", source, ".xlsx", sep=""))
      } else {
        ReportSequences <- fread(paste(patient_dir, "/Resultados_KRAKEN/report_", de_host_file, ".sequences", sep=""),  header = FALSE,sep = "\t")
        tabla_otus <-  read_excel(paste(patient_dir, "/Resultados_KRAKEN/TablaOTUS", conE,"_", id,de_host_file, "_", source, ".xlsx", sep=""))

      }

    } else {
      if(de_host == "Biota") {
        list_files <- list.files(patient_dir, full.names = TRUE)
        path_report <- list_files[which(grepl(sprintf("%s.micriobiome-classification-report.tsv", de_host_file), list_files))]
        ReportSequences <- fread(path_report, header = FALSE, sep = "\t")
        tabla_otus <-  read_excel(paste(patient_dir, "/TablaOTUS_", id, "_trimmed_", source, ".xlsx", sep=""))
      } else {
        list_files <- list.files(paste0(patient_dir, "/Resultados_DRAGEN"), full.names = TRUE)
        path_report <- list_files[which(grepl(sprintf("%s.DRAGEN-report.tsv", de_host_file), list_files))]
        ReportSequences <- fread(path_report, header = FALSE, sep = "\t")
        tabla_otus <-  read_excel(paste(patient_dir, "/Resultados_DRAGEN/TablaOTUS", conE,"_", id,de_host_file, "_", source, ".xlsx", sep=""))
      }

    }


    if(any(tabla_otus[,ncol(tabla_otus)] == 0)) {
      tabla_otus <- tabla_otus[-which(tabla_otus[,ncol(tabla_otus)] == 0),]
    }

    unclassified <- as.numeric(ReportSequences[which(ReportSequences$V6 == "unclassified"), 1])
    cants_original <- data.frame("ID" = id,
                                 "Domain" = length(unique(tabla_otus$Domain)),
                                 "Kingdom" = length(unique(tabla_otus$Kingdom)),
                                 "Phylum"  = length(unique(tabla_otus$Phylum)),
                                 "Class" = length(unique(tabla_otus$Class)),
                                 "Order" = length(unique(tabla_otus$Order)),
                                 "Family" = length(unique(tabla_otus$Family)),
                                 "Genus" = length(unique(tabla_otus$Genus)),
                                 "Species" = length(unique(tabla_otus$Species)),
                                 "SubSpecies" = length(unique(tabla_otus$SubSpecies)),
                                 "Unclassified" = paste(unclassified, "%", sep=" "),
                                 "Source" = source)

    counts_Tax[i,] <- cants_original[1,]
    i <- i+1

  }

  counts_Tax_report <- as.data.frame(t(counts_Tax))
  counts_Tax_report <- data.frame( "Niveles_Taxonomicos" = rownames(counts_Tax_report), counts_Tax_report)
  colnames(counts_Tax_report) <- counts_Tax_report[1,]
  counts_Tax_report <- counts_Tax_report[-c(1, nrow(counts_Tax_report)-1, nrow(counts_Tax_report)) ,]
  colnames(counts_Tax_report)[1] <- "Niveles_Taxonomicos"
  colnames(counts_Tax_report)[-1] <- gsub("_trimmed", "", colnames(counts_Tax_report)[-1])

  library(openxlsx)
  write.xlsx(counts_Tax_report, file = sprintf("%s/trimmed/Resultados_KRAKEN/CountsTaxonomyLevels_%s.xlsx", patients_dir, source))

  return(counts_Tax)
}

#' @title Report Counts Taxonomy
#' @import readxl
#' @import data.table
#' @import openxlsx
#' @description kvjdfnkvjdf
#' @param patient_dir is a dataframe fskbfjsf
#' @param source kjsnfn
#' @return cants_original dataframe that contains how many subjects from the same taxonomic level were found.
#' @examples MetadataB <- as.data.frame(read_excel("~/Daniela/Biota/PipelineBiota/data/Metadata-soloColumnasUsables.xlsx"))
#' reportCountsTax(id = "37", MetadataB = MetadataB, patients_dir = patients_dir)
#' @export

reportCountsTax <- function(id, MetadataB, patients_dir, de_host,  conEukaryota) {

  Biomarcadores_Cantidades <- as.data.frame(read_excel(sprintf("%s/Stats_Cantidades_RangoEtario.xlsx", pipe_data)))
  colnames(Biomarcadores_Cantidades)[2] <- "Rango"
  colnames(Biomarcadores_Cantidades)[5] <- "Q1"
  colnames(Biomarcadores_Cantidades)[6] <- "Q3"

    #pedir id
    #id <- readline(prompt = "Ingrese el ID del paciente que desea analizar: ")
  list_dirs <- list.dirs(patients_dir, recursive = FALSE, full.names = TRUE)
  patient_dir <- path.expand(sprintf("%s/%s", patients_dir, id))

  if(!(patient_dir %in% list_dirs)) {
    return("This id is not on your patients folder")
  }

  #Counts_todos <- counts_Tax(patients_dir = patients_dir, source = "KRAKEN")
  #Counts_todos <- as.data.frame(t(Counts_todos))
  #colnames(Counts_todos) <- Counts_todos[1,]
  #Counts_todos <- Counts_todos[-c(1, nrow(Counts_todos)-1, nrow(Counts_todos) ),]
  #Counts_todos <- cbind(Niveles_Taxonomicos = rownames(Counts_todos), Counts_todos)
  #rownames(Counts_todos) <- NULL

  #intento solucionar:
  Counts_todos <- counts_Tax(patients_dir = patient_dir, source = "KRAKEN", de_host = de_host,  conEukaryota = conEukaryota)
  Counts_todos <- as.data.frame(t(Counts_todos))
  colnames(Counts_todos) <- Counts_todos[1,]
  Counts_todos <- cbind(Niveles_Taxonomicos = rownames(Counts_todos), Counts_todos)
  Counts_todos <- Counts_todos[-c(1, nrow(Counts_todos)-1, nrow(Counts_todos) ),]
  rownames(Counts_todos) <- NULL

  if(any(grepl("_KRAKEN", colnames(Counts_todos)))) {
    colnames(Counts_todos) <- gsub("_KRAKEN", "", colnames(Counts_todos))
  } else if(any(grepl("_DRAGEN", colnames(Counts_todos)))) {
    colnames(Counts_todos) <- gsub("_DRAGEN", "", colnames(Counts_todos))
  } else if(any(grepl("_trimmed", colnames(Counts_todos)))) {
    colnames(Counts_todos) <- gsub("_trimmed", "", colnames(Counts_todos))
  }

  cantidades <- Counts_todos[,which(colnames(Counts_todos) == id)]
  cantidades <- as.data.frame(cantidades)
  cantidades$Taxonomia <- Counts_todos$Niveles_Taxonomicos
  colnames(cantidades)[1] <- id
  cantidades <- cantidades[,c(2,1)]

    #Veo a qué rango etario se aproxima más ----------------
  df_merged <- merge(cantidades, Biomarcadores_Cantidades, by = "Taxonomia")

  edad <- MetadataB$Edad[which(MetadataB$ID == id)]
  rango_real <- MetadataB$`Rango etario`[which(MetadataB$ID == id)]

    # Digo si es alto/Eq/bajo segun su edad
  df_rango_real <- df_merged[which(df_merged$Rango == rango_real),]
  df_rango_real[,2] <- as.numeric(df_rango_real[,2] )
  library(dplyr)

  df_rango_real <- df_rango_real %>%
    mutate(Resultados = ifelse(df_rango_real[,2] > Q1 & df_rango_real[,2] < Q3 | df_rango_real[,2] == Q1 | df_rango_real[,2] == Q3, "Equilibrado",
                                 ifelse(df_rango_real[,2] < Q1, "Bajo", "Alto")))

  colnames(df_rango_real)[2] <- "Valores_Muestra"
  colnames(df_rango_real)[4] <- "Equilibrio"
  df_rango_real$Valores_Muestra <- as.numeric(df_rango_real$Valores_Muestra)
  str(df_rango_real)
  library(ggplot2)
  library(tidyr)

  orden_taxonomia <- c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "SubSpecies")

  df_rango_real$Taxonomia <- factor(df_rango_real$Taxonomia, levels = orden_taxonomia)
  df_rango_real <- df_rango_real[order(df_rango_real$Taxonomia), ]

  # Reorganizar los datos al formato "largo"
  df_largo <- df_rango_real %>%
    pivot_longer(cols = c(Valores_Muestra, Equilibrio),  # Las columnas que se quieren reorganizar
                 names_to = "Tipo",         # Nombre de la nueva columna que indicará si es 'valores' o 'Media'
                 values_to = "Valor")       # Nombre de la columna que contendrá los valores

  # Crear el gráfico de barras
  plot <- ggplot(df_largo, aes(x = Taxonomia, y = Valor, fill = Tipo)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +  # Barras más delgadas
    labs(title = "",
         x = "Nivel taxonómico", y = "Cantidad") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set2") +  # Usar la misma paleta suave
    theme(
      # Personalización del texto del eje x
      axis.text.x = element_text(
        angle = 45,         # Mantener las etiquetas giradas
        hjust = 1,
        size = 10,          # Tamaño similar al segundo gráfico
        face = "italic",    # Texto en cursiva
        color = "gray30"    # Color gris suave
      ),
      # Personalización del texto del eje y
      axis.text.y = element_text(
        size = 10,          # Tamaño igual
        face = "italic",    # Texto en cursiva
        color = "gray30"    # Color gris suave
      ),
      # Leyenda a la derecha
      legend.position = "right",
      legend.text = element_text(size = 10),  # Texto de leyenda más grande
      legend.key.size = unit(0.6, "cm"),  # Tamaño de íconos en la leyenda
      plot.margin = margin(10, 10, 10, 10)  # Ajustar márgenes
    ) +
    # Modificar el nombre de un elemento en la leyenda
    scale_fill_discrete(labels = c("Equilibrio" = "Equilibrio", "Valores_Muestra" = "Valores Muestra")) +
    # Eliminar el título de las leyendas
    guides(fill = guide_legend(title = NULL))

  print(plot)
  colnames(df_rango_real)[6] <- "Min"
  colnames(df_rango_real)[7] <- "Max"
  colnames(df_rango_real)[8] <- "Resultados"

  #Hacer barra ---------------------------------------------------------------------------------
  posiciones <- 1:10 # Posiciones de "Alcohol", "Tabaco" y "Cosméticos"
  niveles <- c("Bajo", "","", "", "", "", "", "", "", "Equilibrado")

  # Definir colores
  color_barra_negra <- "#4D4D4D"  # Gris oscuro para la barra general
  color_barra_roja <- "#69B3E7"   # Rojo para el desbalance

  # Crear el gráfico
  equilibrios <- df_rango_real[which(df_rango_real$Resultados == "Equilibrado"), ]
  n_equilibrios <- nrow(equilibrios) +1
  score_diversidad <- nrow(equilibrios)/nrow(df_rango_real)

  # Crear una nueva columna 'Score_normalizado' entre 0 y 1
  no_equilibrados <- df_rango_real
  no_equilibrados$Score_normalizado <- ifelse(
    no_equilibrados$Valores_Muestra > no_equilibrados$Max,
    (no_equilibrados$Valores_Muestra - no_equilibrados$Max) / (no_equilibrados$Valores_Muestra),  # Escalar sobre Max
    ifelse(
      no_equilibrados$Valores_Muestra < no_equilibrados$Min,
      abs(no_equilibrados$Valores_Muestra - no_equilibrados$Min) / (no_equilibrados$Min),  # Escalar bajo Min
      0  # Si está dentro del rango, score es 0
    )
  )

  no_equilibrados$Score <- 1 - no_equilibrados$Score_normalizado
  score_final <- round(sum(no_equilibrados$Score) / nrow(no_equilibrados),2)
  df_bar <- data.frame(x = c(0, 1), y = c(1, 1))

  bar_plot <- ggplot(df_bar, aes(x = x, y = y)) +
    # Barra de 0 a 1
    geom_segment(aes(x = 0, xend = 1, y = 1, yend = 1), color = "#4D4D4D", size = 3, alpha = 0.8) +
    # Segmento en lugar de triángulo
    geom_segment(aes(x = score_final - 0.02, xend = score_final + 0.02, y = 1, yend = 1),
                 color = "#69B3E7", size = 6, lineend = "round", alpha = 0.7) +
    # Texto estilizado del score
    geom_text(aes(x = score_final, y = 1.01, label = round(score_final, 2)),
              vjust = -1, size = 4, color = "#69B3E7", fontface = "bold") +
    # Etiquetas personalizadas en el eje x
    geom_text(aes(x = 0, y = 0.95, label = "Desequilibrado"), color = "#4D4D4D", size = 4, fontface = "bold") +
    geom_text(aes(x = 1, y = 0.95, label = "Equilibrado"), color = "#4D4D4D", size = 4,  fontface = "bold") +
    scale_x_continuous(breaks = c(0, 1), limits = c(-0.1, 1.1)) +
    #scale_x_continuous(breaks = c(0, 1), labels = c("Desequilibrado", "Equilibrado"), limits = c(-0.1, 1.1)) +
    ylim(0.9, 1.1) +  # Limitar eje y para espacio de texto
    theme_minimal() +  # Tema limpio

    # Estilo estético extraído
    theme(
      # Fondo más suave
      panel.background = element_rect(fill = "#F9F9F9", color = NA),
      plot.background = element_rect(fill = "#F9F9F9", color = NA),
      # Quitar título y etiquetas del eje y
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      # Ajustar el estilo de los ejes
      #axis.line.x = element_line(color = "#A9A9A9", size = 0.8),  # Eje x más marcado
      # Ajustar la posición de los ticks del eje x
      axis.text.x = element_text(color = "#4D4D4D", size = 10, face = "bold"),
      # Texto del título y etiquetas del eje x
      axis.title.x = element_text(size = 1, face = "bold", color = "#4D4D4D")
    )


  print(bar_plot)
  ggsave(sprintf("%s/cantidades_barra.png", patient_dir), plot = bar_plot, width = 6, height = 1.2, dpi = 300)

  if(score_final > 0.74 & score_final < 0.9) {
    dev_cantidades <- "y hemos encontrado levemente menor variabilidad de bacterias de lo habitual"
    mensaje_cantidades <- "La variabilidad de microorganismos detectados en tu piel se encuentra cercana a la esperada según tu rango etario."
    res <- "levemente desequilibrado"

  } else if (score_final < 0.75 ) {
    dev_cantidades <- "y hemos encontrado menor variabilidad de bacterias de lo habitual"
    mensaje_cantidades <- "La variabilidad de microorganismos detectados en tu piel es menor a la esperada según tu rango etario."
    res <- "desequilibrado"
  } else if (score_final >= 0.9 ) {
    dev_cantidades <- "y hemos encontrado una variabilidad de bacterias óptima según tu rango etario"
    mensaje_cantidades <- "¡La variabilidad de microorganismos detectados en tu piel es óptima según tu rango etario!."
    res <- "equilibrado"
  }

  df_rango_real <- df_rango_real[-c(1,2),-c(3,4,5)]
  df_rango_real$Taxonomia <- c("Filos", "Clases", "Órdenes", "Familias", "Géneros", "Especies", "SubEspecies")

  dev_cantidades <- sprintf("Según tus resultados, tu microbioma se encuentra %s con respecto a la variabilidad de microorganismos ", res)

  write.xlsx(df_rango_real, file = sprintf("%s/Tabla_CantidadesTax.xlsx", patient_dir))
  return(list(df_rango_real, bar_plot, mensaje_cantidades, dev_cantidades))

}
