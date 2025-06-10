#' @title Indicador Fototipo
#' @import readxl
#' @import data.table
#' @import openxlsx
#' @import ggplot2
#' @import gridExtra
#' @description It indicates if the skin shows damage because of alcohol drinking by analyzing the relative abundance of determined species.
#' @param id patient's id
#' @param MetadataB dataframe that contains the metadata information. id must be within Metadata$ID.
#' @return list with estados_max (drinking status predicted), mensaje_alcohol (message with the prediction), barra_plot (plot that will be include in the report).
#' @export

indicadorFototipo <- function(id, MetadataB) {

  Biomarcadores_Fototipo <- read_excel("~/Daniela/Biota/PipelineBiota/data/Biomarcadores_Fototipo_sinoutliers.xlsx")

  #pedir id
  #id <- readline(prompt = "Ingrese el ID del paciente que desea analizar: ")
  list_dirs <- list.dirs(patients_dir, recursive = FALSE, full.names = TRUE)
  patient_dir <- path.expand(sprintf("%s/%s", patients_dir, id))

  if(!(patient_dir %in% list_dirs)) {
    return("This id is not on your patients folder")
  }

  otus <-  read_excel(sprintf("%s/trimmed/TablaOTUS_%s_trimmed_KRAKEN.xlsx", patient_dir, id))
  load("~/Daniela/Biota/Muestras/73m/SubEspecies_AUSAR.RData")

  # me quedo con las 138 subespecies y calculo AR -------------------------
  otus_subespecies <- otus[which(otus$SubSpecies %in% subespecies_ausar),]

  AR_subespecies <- otus_subespecies[,9:10]
  AR_biomarcadores <- AR_subespecies
  colnames(AR_biomarcadores)[1] <- "Subespecie"
  AR_biomarcadores[,2:ncol(AR_biomarcadores)] <- prop.table(as.matrix(AR_biomarcadores[,2:ncol(AR_biomarcadores)]), margin = 2) * 100
  sum(AR_biomarcadores[,2:ncol(AR_biomarcadores)])

  df_merged <- merge(AR_biomarcadores, Biomarcadores_Fototipo, by = "Subespecie")
  #Extraigo los biomarcadores ---------------
  # Digo si es alto/Eq/bajo segun su edad

  fototipos <- unique(df_merged$FacilidadBroncearse)
  df_1 <- df_merged[which(df_merged$FacilidadBroncearse == fototipos[1]),]
  df_1 <- df_1 %>%
    mutate(Resultados = ifelse(df_1[,2] > Q1 & df_1[,2] < Q3, "Equilibrado",
                               ifelse(df_1[,2] < Q1, "Bajo", "Alto")))
  t_d <- as.data.frame(table(df_1$Resultados))

  df_2 <- df_merged[which(df_merged$FacilidadBroncearse == fototipos[2]),]
  df_2 <- df_2 %>%
    mutate(Resultados = ifelse(df_2[,2] > Q1 & df_2[,2] < Q3, "Equilibrado",
                               ifelse(df_2[,2] < Q1, "Bajo", "Alto")))
  t_d2 <- as.data.frame(table(df_2$Resultados))

  df_3 <- df_merged[which(df_merged$FacilidadBroncearse == fototipos[3]),]
  df_3 <- df_3 %>%
    mutate(Resultados = ifelse(df_3[,2] > Q1 & df_3[,2] < Q3, "Equilibrado",
                               ifelse(df_3[,2] < Q1, "Bajo", "Alto")))
  t_d3 <- as.data.frame(table(df_3$Resultados))

  #en cual hay mas equilibrados:
  #Cuando no hay equilibrados
  if(length(nchar(t_d$Freq[which(t_d$Var1 == "Equilibrado")])) == 0 ) {
    eq <- data.frame( "Var1" = "Equilibrado", "Freq" = 0)
    t_d <- rbind(t_d, eq)
  }
  if (length(nchar(t_d2$Freq[which(t_d2$Var1 == "Equilibrado")])) == 0) {
    eq <- data.frame( "Var1" = "Equilibrado", "Freq" = 0)
    t_d2 <- rbind(t_d2, eq)
  }
  if (length(nchar(t_d3$Freq[which(t_d3$Var1 == "Equilibrado")])) == 0) {
    eq <- data.frame( "Var1" = "Equilibrado", "Freq" = 0)
    t_d3 <- rbind(t_d3, eq)
  }


  # Crear un dataframe con las frecuencias de "Equilibrado" y los dataframes correspondientes
  equilibrado_data <- data.frame(
    fototipos = fototipos,
    #df_piel = list(df_1, df_2, df_3), # Lista de dataframes
    freq = c(
      t_d$Freq[which(t_d$Var1 == "Equilibrado")],
      t_d2$Freq[which(t_d2$Var1 == "Equilibrado")],
      t_d3$Freq[which(t_d3$Var1 == "Equilibrado")]
    )
  )

  # Ordenar por la frecuencia de mayor a menor
  equilibrado_data <- equilibrado_data[order(equilibrado_data$freq, decreasing = TRUE), ]

  # El dataframe con la mayor frecuencia de "Equilibrado"
  piel_final <- equilibrado_data$fototipos[1]

  # El dataframe con la menor frecuencia de "Equilibrado"
  piel_menor <- equilibrado_data$fototipos[length(equilibrado_data$fototipos)]


  # Creamos el gráfico
  library(ggplot2)

  # Definir colores más suaves y estéticos
  color_barra_negra <- "#4D4D4D"  # Gris oscuro para la barra general
  color_barra_azul <- "#69B3E7"   # Azul pastel para la barra del rango

  # Definir las posiciones correspondientes a cada tipo de piel en el eje x
  fototipos <- c("2", "3", "4")
  posiciones_fototipos <- c(1, 2, 3)  # Asignar una posición para cada piel

  # Crear un dataframe con las frecuencias de "Equilibrado" y las posiciones de las pieles
  equilibrado_data <- data.frame(
    fototipos = fototipos,
    freq = c(
      t_d$Freq[which(t_d$Var1 == "Equilibrado")],
      t_d2$Freq[which(t_d2$Var1 == "Equilibrado")],
      t_d3$Freq[which(t_d3$Var1 == "Equilibrado")]
    ),
    posiciones = posiciones_fototipos
  )

  # Encontrar la frecuencia máxima de "Equilibrados"
  max_freq <- max(equilibrado_data$freq)

  # Verificar si hay un empate entre dos o más tipos de piel
  fototipos_max <- equilibrado_data$fototipos[equilibrado_data$freq == max_freq]
  pos_fototipos_max <- equilibrado_data$posiciones[equilibrado_data$freq == max_freq]


  # Definir la posición del segmento azul:
  # Si hay empate (2 pieles con la misma frecuencia máxima), colocar el segmento entre ellas

  if(length(pos_fototipos_max) == 3) {
    fototipos_max <- "3"
    pos_fototipos_final <- 2    # Si no hay empate, usar la posición normal
    mensaje_fototipos <- sprintf("Según los biomarcadores que utilizamos en BiotaLife para definir el fototipo, tu microbioma indica que tu fototipo es %s", fototipos_max)

  } else if (length(pos_fototipos_max) == 2) {
    pos_fototipos_final <- mean(pos_fototipos_max)  # Posicionar el segmento entre las dos categorías
    mensaje_fototipos <- sprintf("Según los biomarcadores que utilizamos en BiotaLife para definir el fototipo, tu microbioma indica que tu piel se encuentra entre los fototipos %s", paste0(fototipos_max, collapse = " y "))

  } else {
    pos_fototipos_final <- pos_fototipos_max[1]     # Si no hay empate, usar la posición normal
    mensaje_fototipos <- sprintf("Según los biomarcadores que utilizamos en BiotaLife para definir el fototipo, tu microbioma indica que tu fototipo es %s", fototipos_max)
  }

  # Crear el gráfico
  barra_plot <- ggplot() +
    # Barra negra general que representa todas las pieles
    geom_segment(aes(x = 1, xend = 3, y = 1, yend = 1), color = color_barra_negra, size = 3, alpha = 0.8) +
    # Barra azul que marca la piel indicada por piel_final o el empate
    geom_segment(aes(x = pos_fototipos_final - 0.2, xend = pos_fototipos_final + 0.2, y = 1, yend = 1),
                 color = color_barra_azul, size = 8, lineend = "round", alpha = 0.7) +
    # Ajustes estéticos
    scale_x_continuous(limits = c(0.8, 3.2), breaks = posiciones_fototipos, labels = fototipos) +
    scale_y_continuous(limits = c(0.9999, 1.0001)) +
    theme_minimal(base_size = 15) +
    theme(
      # Fondo más suave
      panel.background = element_rect(fill = "#F9F9F9", color = NA),
      plot.background = element_rect(fill = "#F9F9F9", color = NA),
      # Quitar título y etiquetas del eje y
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      # Ajustar el estilo de los ejes
      axis.line.x = element_line(color = "#A9A9A9", size = 0.8),  # Eje x más marcado
      # Ajustar la posición de los ticks del eje x
      axis.text.x = element_text(color = "#4D4D4D", size = 10, face = "bold"),
      # Texto del título y etiquetas del eje x
      axis.title.x = element_text(size = 1, face = "bold", color = "#4D4D4D")
    ) +
    # Etiqueta del eje x
    labs(x = NULL)

  # Mostrar el gráfico
  print(barra_plot)

  ggsave(sprintf("%s/fototipo_barra.png", patient_dir), plot = barra_plot, width = 6, height = 1, dpi = 300)

  return(list(fototipos_max, mensaje_fototipos, barra_plot))

}
