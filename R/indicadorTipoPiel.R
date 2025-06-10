#' @title Indicador Tipo de Piel
#' @import readxl
#' @import data.table
#' @import openxlsx
#' @import ggplot2
#' @import gridExtra
#' @description It indicates if the id's skin is oil, dry or mixed by analyzing the relative abundance of determined species.
#' @param id patient's id
#' @param MetadataB dataframe that contains the metadata information. id must be within Metadata$ID.
#' @return list with pieles_max (type of skin predicted), mensaje_piel (message with the prediction), barra_plot (plot that will be include in the report).
#' @export

indicadorTipoPiel <- function(patient_dir, MetadataB) {

  id <- basename(patient_dir)
  patients_dir <- dirname(patient_dir)

  #Biomarcadores_PielGrasa <- read_excel("~/Daniela/Biota/PipelineBiota/data/Biomarcadores_TipodePiel_sinoutliers.xlsx")
  Biomarcadores_PielGrasa <- read_excel("/media/4tb2/Daniela/Biota/PipelineBiota-master/data/NEW-Biomarcadores_TipodePiel_sinoutliers.xlsx")
  #Biomarcadores_PielGrasa <- read_excel("~/Daniela/Biota/PipelineBiota/data/NEW-sinfagos-Biomarcadores_TipodePiel_sinoutliers.xlsx")


  #pedir id
  #id <- readline(prompt = "Ingrese el ID del paciente que desea analizar: ")
  #list_dirs <- list.dirs(patients_dir, recursive = FALSE, full.names = TRUE)
  #patient_dir <- path.expand(sprintf("%s/%s", patients_dir, id))

  #if(!(patient_dir %in% list_dirs)) {
  #  return("This id is not on your patients folder")
  #}

  otus <-  read_excel(sprintf("%s/trimmed/Resultados_KRAKEN/TablaOTUS_%sBo_KRAKEN.xlsx", patient_dir, id))
  #load("~/Daniela/Biota/Muestras/73m/SubEspecies_AUSAR.RData")
  load("/media/4tb2/Daniela/Biota/PipelineBiota-master/data/SubSpecies_AUSAR.RData")
  

  # me quedo con las 138 subespecies y calculo AR -------------------------
  subespecies_ausar <- unique(Biomarcadores_PielGrasa$Subespecie)
  otus_subespecies <- otus[which(otus$SubSpecies %in% subespecies_ausar),]

  AR_subespecies <- otus_subespecies[,9:10]
  AR_biomarcadores <- AR_subespecies
  colnames(AR_biomarcadores)[1] <- "Subespecie"
  AR_biomarcadores[,2:ncol(AR_biomarcadores)] <- prop.table(as.matrix(AR_biomarcadores[,2:ncol(AR_biomarcadores)]), margin = 2) * 100

  sum(AR_biomarcadores[,2:ncol(AR_biomarcadores)])

  df_merged <- merge(AR_biomarcadores, Biomarcadores_PielGrasa, by = "Subespecie")
  #Extraigo los biomarcadores ---------------
  # Digo si es alto/Eq/bajo segun su edad

      tiposdepiel <- c("Piel Grasa", "Piel Mixta", "Piel Seca")

      #Para piel grasa:
      df_1 <- df_merged[which(df_merged$Tipodepiel == tiposdepiel[1]),]
      df_1 <- df_1 %>%
        mutate(Resultados = ifelse(df_1[,2] >= Q1 & df_1[,2] <= Q3, "Equilibrado",
                                   ifelse(df_1[,2] < Q1, "Bajo", "Alto")))

      #Se hace un ajuste para los casos en los cuales es elevado para pieles grasa:
      elevados_grasa <- c("Cutibacterium acnes subsp. acnes", "Cutibacterium acnes TypeIA2 P.acn17", "Staphylococcus epidermidis RP62A")
      filas_a_cambiar <- which(df_1$Resultados == "Alto" & df_1$Subespecie %in% elevados_grasa)
      # Cambiar el valor de 'Resultados' en esas filas a "Equilibrado"
      df_1$Resultados[filas_a_cambiar] <- "Equilibrado"
      t_d <- as.data.frame(table(df_1$Resultados))

      #Para piel mixta:
      df_2 <- df_merged[which(df_merged$Tipodepiel == tiposdepiel[2]),]
      df_2 <- df_2 %>%
        mutate(Resultados = ifelse(df_2[,2] >=  Q1 & df_2[,2] <= Q3, "Equilibrado",
                                   ifelse(df_2[,2] < Q1, "Bajo", "Alto")))
      t_d2 <- as.data.frame(table(df_2$Resultados))

      df_3 <- df_merged[which(df_merged$Tipodepiel == tiposdepiel[3]),]
      df_3 <- df_3 %>%
        mutate(Resultados = ifelse(df_3[,2] >= Q1 & df_3[,2] <= Q3, "Equilibrado",
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



  # Creamos el gráfico
  library(ggplot2)

      # Definir colores más suaves y estéticos
  color_barra_negra <- "#4D4D4D"  # Gris oscuro para la barra general
  color_barra_azul <- "#69B3E7"   # Azul pastel para la barra del rango

      # Definir las posiciones correspondientes a cada tipo de piel en el eje x
  pieles <- tiposdepiel
  posiciones_piel <- c(1,2,3)  # Asignar una posición para cada piel

      # Crear un dataframe con las frecuencias de "Equilibrado" y las posiciones de las pieles
  equilibrado_data <- data.frame(
    tiposdepiel = pieles,
    freq = c(
      t_d$Freq[which(t_d$Var1 == "Equilibrado")],
      t_d2$Freq[which(t_d2$Var1 == "Equilibrado")],
      t_d3$Freq[which(t_d3$Var1 == "Equilibrado")]
    ),
    posiciones = posiciones_piel
    )


  # Encontrar la frecuencia máxima de "Equilibrados"
  max_freq <- max(equilibrado_data$freq)

  # Verificar si hay un empate entre dos o más tipos de piel
  pieles_max <- equilibrado_data$tiposdepiel[equilibrado_data$freq == max_freq]
  pos_piel_max <- equilibrado_data$posiciones[equilibrado_data$freq == max_freq]

  #Hacer una salvedad para mixta-seca:
  if(pieles_max == "Piel Seca" & equilibrado_data$freq[which(equilibrado_data$tiposdepiel == "Piel Mixta")]>5) {
    pieles_max <- "Piel Mixta"
    pos_piel_max <- 2
  }


  # Definir la posición del segmento azul:
  # Si hay empate (2 pieles con la misma frecuencia máxima), colocar el segmento entre ellas

  if(length(pos_piel_max) == 3) {
    pieles_max <- "Piel Mixta"
    pos_piel_final <- 2    # Si no hay empate, usar la posición normal
    tipo <- strsplit(pieles_max, split = " ")[[1]][2]
    mensaje_piel <- sprintf("El color celeste indica que tu piel es de tipo %s, según los biomarcadores que utilizamos en Biotalife.", tipo)
    dev_piel <- ", presenta características de Piel Mixta"

  } else if (length(pos_piel_max) == 2) {
    pos_piel_final <- mean(pos_piel_max)  # Posicionar el segmento entre las dos categorías
    mensaje_piel <- sprintf("El color celeste indica que tu piel se encuentra entre %s, según los biomarcadores que utilizamos en Biotalife.", paste0(pieles_max, collapse = " y "))
    dev_piel <- sprintf(", presenta características de %s", paste0(pieles_max, collapse = " y "))

  } else {
    pos_piel_final <- pos_piel_max[1]     # Si no hay empate, usar la posición normal
    tipo <- strsplit(pieles_max, split = " ")[[1]][2]
    mensaje_piel <- sprintf("El color celeste indica que tu piel es de tipo %s, según los biomarcadores que utilizamos en Biotalife.", tipo)
    dev_piel <- sprintf(", presenta características de Piel %s", tipo)
  }

      # Crear el gráfico
  barra_plot <- ggplot() +
        # Barra negra general que representa todas las pieles
    geom_segment(aes(x = 1, xend = 3, y = 1, yend = 1), color = color_barra_negra, size = 3, alpha = 0.8) +
        # Barra azul que marca la piel indicada por piel_final o el empate
    geom_segment(aes(x = pos_piel_final - 0.2, xend = pos_piel_final + 0.2, y = 1, yend = 1),
               color = color_barra_azul, size = 8, lineend = "round", alpha = 0.7) +
        # Ajustes estéticos
    scale_x_continuous(limits = c(0.8, 3.2), breaks = posiciones_piel, labels = pieles) +
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

  ggsave(sprintf("%s/tipodepiel_barra.png", patient_dir), plot = barra_plot, width = 6, height = 1, dpi = 300)


  return(list(pieles_max, mensaje_piel, barra_plot, dev_piel))

}
