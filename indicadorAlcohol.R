#' @title Indicador Alcohol
#' @import readxl
#' @import data.table
#' @import openxlsx
#' @import ggplot2
#' @import gridExtra
#' @description
#' @description It indicates if the skin shows damage because of alcohol drinking by analyzing the relative abundance of determined species.
#' @param id patient's id
#' @param MetadataB dataframe that contains the metadata information. id must be within Metadata$ID.
#' @return list with estados_max (drinking status predicted), mensaje_alcohol (message with the prediction), barra_plot (plot that will be include in the report).
#' @export

indicadorAlcohol <- function(id, MetadataB) {

  Biomarcadores_Alcohol <- read_excel("~/Daniela/Biota/PipelineBiota/data/Biomarcadores_Alcohol_sinoutliers.xlsx")

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

  df_merged <- merge(AR_biomarcadores, Biomarcadores_Alcohol, by = "Subespecie")
  #Extraigo los biomarcadores ---------------
  # Digo si es alto/Eq/bajo segun su edad

  estados <- unique(df_merged$Alcohol)
  df_1 <- df_merged[which(df_merged$Alcohol == estados[1]),]
  df_1 <- df_1 %>%
    mutate(Resultados = ifelse(df_1[,2] > Q1 & df_1[,2] < Q3, "Equilibrado",
                               ifelse(df_1[,2] < Q1, "Bajo", "Alto")))
  t_d <- as.data.frame(table(df_1$Resultados))

  df_2 <- df_merged[which(df_merged$Alcohol == estados[2]),]
  df_2 <- df_2 %>%
    mutate(Resultados = ifelse(df_2[,2] > Q1 & df_2[,2] < Q3, "Equilibrado",
                               ifelse(df_2[,2] < Q1, "Bajo", "Alto")))
  t_d2 <- as.data.frame(table(df_2$Resultados))


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


  # Creamos el gráfico
  library(ggplot2)

  # Definir colores más suaves y estéticos
  color_barra_negra <- "#4D4D4D"  # Gris oscuro para la barra general
  color_barra_azul <- "#69B3E7"   # Azul pastel para la barra del rango

  # Definir las posiciones correspondientes a cada tipo de piel en el eje x
  estados <- c("SI", "NO")
  posiciones_estados <- c(1, 2)  # Asignar una posición para cada piel

  # Crear un dataframe con las frecuencias de "Equilibrado" y las posiciones de las pieles
  equilibrado_data <- data.frame(
    estados = estados,
    freq = c(
      t_d$Freq[which(t_d$Var1 == "Equilibrado")],
      t_d2$Freq[which(t_d2$Var1 == "Equilibrado")]
    ),
    posiciones = posiciones_estados
  )

  # Encontrar la frecuencia máxima de "Equilibrados"
  max_freq <- max(equilibrado_data$freq)

  # Verificar si hay un empate entre dos o más tipos de piel
  estados_max <- equilibrado_data$estados[equilibrado_data$freq == max_freq]
  pos_estados_max <- equilibrado_data$posiciones[equilibrado_data$freq == max_freq]


  # Definir la posición del segmento azul:
  # Si hay empate (2 pieles con la misma frecuencia máxima), colocar el segmento entre ellas

  if (length(pos_estados_max) == 2) {
    pos_estados_final <- mean(pos_estados_max)  # Posicionar el segmento entre las dos categorías
    mensaje_alcohol <- "Según los biomarcadores que utilizamos en BiotaLife para definir el afecto por ingesta de alcohol, tu microbioma indica que tu piel se encuentra en un punto intermedio"

  } else {
    pos_estados_final <- pos_estados_max[1]     # Si no hay empate, usar la posición normal
    if(pos_estados_final == 2) {
      res <- "balanceado"
    } else {
      res = "desbalanceado"
    }
    mensaje_alcohol <- sprintf("Según los biomarcadores que utilizamos en BiotaLife para definir el efecto por ingesta de alcohol, tu microbioma se encuentra %s", res)
  }

  # Crear el gráfico
  estados_balanceados <- c("Desbalanceado", "Balanceado")
  barra_plot <- ggplot() +
    # Barra negra general que representa todas las pieles
    geom_segment(aes(x = 1, xend = 2, y = 1, yend = 1), color = color_barra_negra, size = 3, alpha = 0.8) +
    # Barra azul que marca la piel indicada por piel_final o el empate
    geom_segment(aes(x = pos_estados_final - 0.2, xend = pos_estados_final + 0.2, y = 1, yend = 1),
                 color = color_barra_azul, size = 6, lineend = "round", alpha = 0.7) +
    # Ajustes estéticos
    scale_x_continuous(limits = c(0.8, 2.2), breaks = posiciones_estados, labels = estados_balanceados) +
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

  ggsave(sprintf("%s/alcohol_barra.png", patient_dir), plot = barra_plot, width = 6, height = 1, dpi = 300)

  return(list(estados_max, mensaje_alcohol, barra_plot))

}
