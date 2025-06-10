#' @title Indicador Virus
#' @import readxl
#' @import data.table
#' @import openxlsx
#' @import ggplot2
#' @import gridExtra
#' @description It indicates if the skin present relative abundance of determined viruses above the normal values.
#' @param id patient's id
#' @param MetadataB dataframe that contains the metadata information. id must be within Metadata$ID.
#' @return mensaje_virus (message with the conclusion).
#' @export

indicadorVirus <- function(id, MetadataB) {

  nivel <- "Species"
  Biomarcadores_Especies_abundantes <- as.data.frame(read_excel(sprintf("/media/4tb2/Daniela/Biota/PipelineBiota-master/data/%sAbundantes_RangoEtario_sinoutliers.xlsx", nivel)))
  colnames(Biomarcadores_Especies_abundantes)[1] <- nivel
  Biomarcadores_Especies_abundantes[,1] <- gsub("\\.", " ", Biomarcadores_Especies_abundantes[,1])

  #pedir id
  #id <- readline(prompt = "Ingrese el ID del paciente que desea analizar: ")
  list_dirs <- list.dirs(patients_dir, recursive = FALSE, full.names = TRUE)
  patient_dir <- path.expand(sprintf("%s/%s", patients_dir, id))

  if(!(patient_dir %in% list_dirs)) {
    return("This id is not on your patients folder")
  }

  otus <-  read_excel(sprintf("%s/trimmed/Resultados_KRAKEN/TablaOTUS_%sBo_KRAKEN.xlsx", patient_dir, id))
  #otus <-  read_excel(sprintf("%s/OTUs_84Pacientes_KRAKEN.xlsx", patients_dir))

  especies_ausar <- unique(Biomarcadores_Especies_abundantes$Species)

  # me quedo con las especies y calculo AR
  otus_especies <- otus[which(otus$Species %in% especies_ausar),]

  AR_especies <- otus_especies[,c(8,10)]
  #AR_especies <- otus_especies[,-c(1,2,3,4,5,6,7,9)]
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
  colnames(AR_biomarcadores)[1] <- nivel
  df_merged <- merge(AR_biomarcadores, Biomarcadores_Especies_abundantes, by = nivel)

  edad <- MetadataB$Edad[which(MetadataB$ID == id)]
  rango_real <- MetadataB$`Rango etario`[which(MetadataB$ID == id)]

  # Digo si es alto/Eq/bajo segun su edad
  df_rango_real <- df_merged[which(df_merged$Rango == rango_real),]

  virus_rojos <- c("Betapapillomavirus 1", "Betapapillomavirus 2", "Betapapillomavirus 3", "Betapapillomavirus 4", "Betapapillomavirus 5", "Human polyomavirus 5", "Merkel cell polyomavirus")

  AR_virus_biom <- df_rango_real[which(df_rango_real$Species %in% virus_rojos),]

  AR_virus_biom$Indicador <- ifelse(AR_virus_biom[,2] > AR_virus_biom$Media*5 | AR_virus_biom[,2] > 1 , "Riesgo", "")


  AR_virus_biom[, c(2,7)] <- round(AR_virus_biom[, c(2,7)], 4)
  AR_virus_biom$Q3 <- AR_virus_biom$Media

  #Por cada virus cuyo indicador diga Riesgo? se la agrega su funcion y su recomendacion: -----------------------

  Recomendacion_Viruses <- as.data.frame(read_excel("/media/4tb2/Daniela/Biota/PipelineBiota-master/data/Viruses_Recomendaciones.xlsx"))
  paraCompleto <- Recomendacion_Viruses[which(Recomendacion_Viruses$Species %in% AR_virus_biom$Species), c("Species", "Que hace", "Recomendaciones")]
  mensaje_paraCompleto <- c()

  virus_riesgo <- AR_virus_biom$Species[which(AR_virus_biom$Indicador == "Riesgo")]

  #Una misma funcion y recomendacion para cualquiera de los betapapilomavirus:
  if(any(grepl("Betapapillomavirus", virus_riesgo))) {
    funcion <-  paraCompleto$`Que hace`[which(paraCompleto$Species == "Betapapillomavirus 1")]
    recomendacion <- paraCompleto$Recomendaciones[which(paraCompleto$Species == "Betapapillomavirus 1")]

    if(length(virus_riesgo)>1) {
      viruses <- paste0(virus_riesgo[which(grepl("Betapapillomavirus",virus_riesgo))], collapse = ", ")
      mensaje_paraCompleto <- c(mensaje_paraCompleto, sprintf("Los virus %s cumplen con las siguientes funciones: %s Por lo tanto, se sugieren las siguientes recomendaciones %s", viruses, funcion, recomendacion))
    } else {
      virus <- virus_riesgo
      mensaje_paraCompleto <- c(mensaje_paraCompleto, sprintf("El virus %s cumple con las siguientes funciones: %s Por lo tanto, se sugieren las siguientes recomendaciones %s", virus, funcion, recomendacion))
    }
  }

  #Si esta el otro virus:
  virus_riesgo_sinbeta <- virus_riesgo[-which(grepl("Betapapillomavirus", virus_riesgo))]
  if(length(virus_riesgo_sinbeta) >0) {
    funcion <-  paraCompleto$`Que hace`[which(paraCompleto$Species == "Human polyomavirus 5")]
    recomendacion <- paraCompleto$Recomendaciones[which(paraCompleto$Species == "Human polyomavirus 5")]
    mensaje_paraCompleto <- c(mensaje_paraCompleto, sprintf("El virus %s cumple con las siguientes funciones: %s Por lo tanto, se sugieren las siguientes recomendaciones %s", virus_riesgo_sinbeta, funcion, recomendacion))
  }


  #--------------------------------------------------------------------

  # Crear el gráfico
  library(ggplot2)

  # Datos ficticios
  posiciones <- 1:6 # Posiciones de "Alcohol", "Tabaco" y "Cosméticos"
  virus <- virus_rojos[1:6]

  # Definir colores
  color_barra_negra <- "#4D4D4D"  # Gris oscuro para la barra general
  color_barra_roja <- "#E74C3C"   # Rojo para el desbalance

  # Crear el gráfico
  presentes <- AR_virus_biom[which(AR_virus_biom$Indicador == "Riesgo"), ]
  #presentes <- AR_virus_biom
  presentes$x <- which(AR_virus_biom$Indicador == "Riesgo") - 0.2
  presentes$xend <- which(AR_virus_biom$Indicador == "Riesgo") + 0.2


  barra_plot <- ggplot() +
    # Barra negra general que representa todas las categorías
    geom_segment(aes(x = 1, xend = 6, y = 1, yend = 1), color = color_barra_negra, size = 2, alpha = 0.8) +
    # Barra roja para los que están en "Riesgo"
    geom_segment(data = presentes,
                 aes(x = x,
                     xend = xend,
                     y = 1, yend = 1),
                 color = color_barra_roja, size = 6, lineend = "round", alpha = 0.7) +

    # Ajustes estéticos
    scale_x_continuous(limits = c(0.8, 6.2), breaks = posiciones, labels = virus) +
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
      axis.text.x = element_text(color = "#4D4D4D", size = 6, face = "bold"),
      # Texto del título y etiquetas del eje x
      axis.title.x = element_text(size = 1, face = "bold", color = "#4D4D4D")
    )


  print(barra_plot)

  ggsave(sprintf("%s/virus_barra.png", patient_dir), plot = barra_plot, width = 7, height = 1, dpi = 300)
  dev.off()

  if(nrow(presentes)>0) {
    mensaje_virus <- "El color rojo destaca la presencia de los virus que podrían afectar la salud de tu piel. Se recomienda consultar con un especialista."
    dev_virus <- " y presenta virus que pueden amenazar la salud de tu piel."
  } else {
    mensaje_virus <- "La ausencia de estos virus implica ausencia de amenazas hacia tu microbioma de piel."
    dev_virus <- "."
  }

  AR_virus_biom <- AR_virus_biom[, -c(3,4,5)]
  colnames(AR_virus_biom)[c(2,3,4)] <- c("AR", "Min", "Max" )
  AR_virus_biom[, c(2,3,4)] <- round(AR_virus_biom[, c(2,3,4)], 4)
  rownames(AR_virus_biom) <- NULL

  virus_presentes <- AR_virus_biom$Species[which(AR_virus_biom$Indicador == "Riesgo")]
  if(length(virus_presentes)>1) {
    virus_presentes <- paste0(virus_presentes, collapse=",")
  }


  return(list(mensaje_virus, dev_virus, mensaje_paraCompleto, AR_virus_biom))

  #return(virus_presentes)

}

