#' @title Indicador Fagos
#' @import readxl
#' @import data.table
#' @import openxlsx
#' @import ggplot2
#' @import gridExtra
#' @description It indicates if the skin shows presence of phages.
#' @param id patient's id
#' @param MetadataB dataframe that contains the metadata information. id must be within Metadata$ID.
#' @return mensaje_fagos (message with the prediction).
#' @export
indicadorFagos <- function(patient_dir, MetadataB) {

  nivel <- "Species"
  Biomarcadores_Especies_abundantes <- as.data.frame(read_excel(sprintf("%s/%sAbundantes_RangoEtario_sinoutliers.xlsx", pipe_data, nivel)))

  #pedir id
  #id <- readline(prompt = "Ingrese el ID del paciente que desea analizar: ")
  #list_dirs <- list.dirs(patients_dir, recursive = FALSE, full.names = TRUE)
  #patient_dir <- path.expand(sprintf("%s/%s", patients_dir, id))

  id <- basename(patient_dir)
  #if(!(patient_dir %in% list_dirs)) {
  #  return("This id is not on your patients folder")
  #}

  otus <-  read_excel(sprintf("%s/trimmed/Resultados_KRAKEN/TablaOTUS_%sBo_KRAKEN.xlsx", patient_dir, id))
  colnames(Biomarcadores_Especies_abundantes)[1] <- nivel
  Biomarcadores_Especies_abundantes[,1] <- gsub("\\.", " ", Biomarcadores_Especies_abundantes[,1])
  especies_ausar <- unique(Biomarcadores_Especies_abundantes$Species)

  # me quedo con las especies y calculo AR
  otus_especies <- otus[which(otus$Species %in% especies_ausar),]

  AR_especies <- otus_especies[,c(8,10)]
  AR_biomarcadores <- AR_especies
  colnames(AR_biomarcadores)[1] <- "Especie"

  #Sumar los que son iguales:
  AR_biomarcadores <- AR_biomarcadores %>%
    # Agrupar por nivel
    group_by(Especie) %>%
    summarise(across(everything(), sum, na.rm = TRUE)) %>%
    ungroup()

  AR_biomarcadores[,2:ncol(AR_biomarcadores)] <- prop.table(as.matrix(AR_biomarcadores[,2:ncol(AR_biomarcadores)]), margin = 2) * 100
  sum(AR_biomarcadores[,2:ncol(AR_biomarcadores)])

  AR_virus_rojos <- AR_biomarcadores[which(grepl("Lactococcus phage", AR_biomarcadores$Especie) | grepl("Propionibacterium phage", AR_biomarcadores$Especie) | grepl("Staphylococcus phage", AR_biomarcadores$Especie) | grepl("Streptococcus phage", AR_biomarcadores$Especie) ),]

  #Agrupo por fago:
  library(stringr)
  library(dplyr)

  AR_fagos <- AR_virus_rojos %>%
    # Crear una nueva columna con las dos primeras palabras de la especie
    mutate(Grupo_Especie = word(Especie, 1, 2)) %>%
    group_by(Grupo_Especie) %>%
    summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%
    ungroup()

  AR_fagos$Indicador <- ifelse(AR_fagos[,2] > 0, "Presente", "")


  # Crear el gráfico
  library(ggplot2)

  # Datos ficticios
  posiciones <- 1:4
  fagos <- c("Lactococcus fago", "Propionibacterium fago", "Staphylococcus fago", "Streptococcus fago")


  # Definir colores
  color_barra_negra <- "#4D4D4D"  # Gris oscuro para la barra general
  color_barra_roja <- "#27AE60"   # Rojo para el desbalance

  # Crear el gráfico
  presentes <- AR_fagos[which(AR_fagos$Indicador == "Presente"), ]
  #presentes <- AR_fagos
  presentes$x <- which(AR_fagos$Indicador == "Presente") - 0.2
  presentes$xend <- which(AR_fagos$Indicador == "Presente") + 0.2

  barra_plot <- ggplot() +
    # Barra negra general que representa todas las categorías
    geom_segment(aes(x = 1, xend = 4, y = 1, yend = 1), color = color_barra_negra, size = 2, alpha = 0.8) +
    # Barra roja para los que están en "Riesgo"
    geom_segment(data = presentes,
                 aes(x = x,
                     xend = xend,
                     y = 1, yend = 1),
                 color = color_barra_roja, size = 6, lineend = "round", alpha = 0.7) +

    # Ajustes estéticos
    scale_x_continuous(limits = c(0.8, 4.2), breaks = posiciones, labels = fagos) +
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
      axis.text.x = element_text(color = "#4D4D4D", size = 7, face = "bold"),
      # Texto del título y etiquetas del eje x
      axis.title.x = element_text(size = 1, face = "bold", color = "#4D4D4D")
    )


  print(barra_plot)

  ggsave(sprintf("%s/fagos_barra.png", patient_dir), plot = barra_plot, width = 7, height = 1, dpi = 300)

  if(nrow(presentes) > 0) {
    mensaje_fagos <- "El color verde destaca la presencia de fagos que ayudan a mantener el equilibrio y salud de tu piel."
  } else {
    mensaje_fagos <- "La ausencia de estos fagos no implica amenazas directas hacia tu microbioma de piel. Pero se ofrecen recomendaciones para fomentar su crecimiento."
  }

  return(mensaje_fagos)

}
