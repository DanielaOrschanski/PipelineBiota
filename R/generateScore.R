#' @title Generate Balance Score
#' @description It generates the report in pdf format that register all the predictions and the analizes of the microbiome composition of a determined patient.
#' @param patients_dir indicated the patient directory that will be analyzed.
#' @export
#' @import rmarkdown
#' @import ggplot2
#' @import kableExtra
#' @import gridExtra
#' @import readxl
#' @import data.table
#' @import openxlsx
generateScore <- function(patient_dir) {

  #subesp_masAbundantes <- read_excel(sprintf("%s/Tabla_SubSpecies_masAbundantes.xlsx", patient_dir))
  #subesp_masAbundantes <- subesp_masAbundantes[1:5,]
  esp_masAbundantes <- read_excel(sprintf("%s/Tabla_Species_masAbundantes.xlsx", patient_dir))
  #esp_masAbundantes <- esp_masAbundantes[1:5,]
  Filos_masAbundantes <- read_excel(sprintf("%s/Tabla_Phylum_masAbundantes.xlsx", patient_dir))
  Genus_masAbundantes <- read_excel(sprintf("%s/Tabla_Genus_masAbundantes.xlsx", patient_dir))

  colnames(Genus_masAbundantes)[1] <- "X"
  colnames(Filos_masAbundantes)[1] <- "X"
  colnames(esp_masAbundantes)[1] <- "X"
  #colnames(subesp_masAbundantes)[1] <- "X"


  no_equilibrados <- rbind(Genus_masAbundantes,
                           Filos_masAbundantes,
                           esp_masAbundantes
  )
  no_equilibrados <- na.omit(no_equilibrados)
  # Crear una nueva columna 'Score_normalizado' entre 0 y 1
  no_equilibrados$Score_normalizado <- ifelse(
    no_equilibrados$AR > no_equilibrados$Max,
    (no_equilibrados$AR - no_equilibrados$Max) / (no_equilibrados$AR),  # Escalar sobre Max
    ifelse(
      no_equilibrados$AR < no_equilibrados$Min,
      abs(no_equilibrados$AR - no_equilibrados$Min) / (no_equilibrados$Min),  # Escalar bajo Min
      0  # Si está dentro del rango, score es 0
    )
  )

  no_equilibrados$Score <- 1 - no_equilibrados$Score_normalizado
  score_final <- round(sum(no_equilibrados$Score) / nrow(no_equilibrados),2)

  library(ggplot2)

  # Crear una base de datos ficticia para la barra de 0 a 1
  df_bar <- data.frame(x = c(0, 1), y = c(1, 1))


  bar_score_plot <- ggplot(df_bar, aes(x = x, y = y)) +
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



  # Mostrar el gráfico
  print(bar_score_plot)
  ggsave(sprintf("%s/score_barra.png", patient_dir), plot = bar_score_plot, width = 6, height = 1.2, dpi = 300)
  dev.off()

  extra <- ifelse(score_final < 1, "En la sección de recomendaciones se sugieren cambios en los hábitos que pueden colaborar a la mejora en el equilibrio de su piel." , "")


  if(score_final < 0.9 & score_final > 0.74) {
    res <- "levemente menor a la esperada"
    dev_score_gral <- sprintf("y hemos encontrado una cantidad de bacterias %s.", res)
  } else if(score_final < 0.75 ) {
    res <- "menor a la esperada"
    dev_score_gral <- sprintf("y hemos encontrado una cantidad de bacterias %s.", res)
  } else if( score_final >= 0.9) {
    res <- "la esperada"
    dev_score_gral <- sprintf("y hemos encontrado una cantidad de bacterias óptima con respecto a la esperada.")
  }

  mensaje_score <- sprintf("La cantidad de microorganismos detectados en tu piel es %s según tu rango etario.", res)


  return(list(mensaje_score, dev_score_gral))

}
