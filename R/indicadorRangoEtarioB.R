#' @title Indicador Rango Etario
#' @import readxl
#' @import data.table
#' @import openxlsx
#' @import ggplot2
#' @import gridExtra
#' @import dplyr
#' @description kvjdfnkvjdf
#' @param id patient's id
#' @param MetadataB dataframe that contains the metadata information. id must be within Metadata$ID.
#' @return list with df_rango_real (dataframe with comparison of relative abundances within the biological age),
#' rango_final (prediction of age), plot1 (pie chart), combined_plot (combination of pie charts, one for each age range),
#' mensaje (message with the resolution), barra_plot (plot that will be included in the report).
#' @export

indicadorRangoEtario <- function(patient_dir, MetadataB) {

  Biomarcadores_RangoEtario <- read_excel(sprintf("%s/Biomarcadores_RangoEtario_sinoutliers.xlsx", pipe_data))

  #pedir id
  #id <- readline(prompt = "Ingrese el ID del paciente que desea analizar: ")
  #list_dirs <- list.dirs(patients_dir, recursive = FALSE, full.names = TRUE)
  #patient_dir <- path.expand(sprintf("%s/%s", patients_dir, id))
  id <- basename(patient_dir)
  print(id)
  
  #if(!(patient_dir %in% list_dirs)) {
  #  return("This id is not on your patients folder")
  #}

  otus <-  read_excel(sprintf("%s/trimmed/Resultados_KRAKEN/TablaOTUS_%sBo_KRAKEN.xlsx", patient_dir, id))
  load(sprintf("%s/SubSpecies_AUSAR.RData", pipe_data))

  # me quedo con las 138 subespecies y calculo AR -------------------------
  otus_subespecies <- otus[which(otus$SubSpecies %in% subespecies_ausar),]
  sub_out <- subespecies_ausar[which(!(subespecies_ausar %in% otus_subespecies$SubSpecies))]
  #AR_SubSpecies_KRAKEN$'1_KRAKEN'[which(AR_SubSpecies_KRAKEN$SubSpecies %in% sub_out)]

  AR_subespecies <- otus_subespecies[,9:10]
  AR_biomarcadores <- AR_subespecies[which(AR_subespecies$SubSpecies %in% Biomarcadores_RangoEtario$Subespecie),]
  colnames(AR_biomarcadores)[1] <- "Subespecie"

  AR_biomarcadores[,2:ncol(AR_biomarcadores)] <- prop.table(as.matrix(AR_biomarcadores[,2:ncol(AR_biomarcadores)]), margin = 2) * 100
  sum(AR_biomarcadores[,2:ncol(AR_biomarcadores)])

  #Extraigo los biomarcadores ---------------

  #Veo a qué rango etario se aproxima más ----------------

  df_merged <- merge(AR_biomarcadores, Biomarcadores_RangoEtario, by = "Subespecie")

  #Probar con solo algunas subespecies:
  #df_merged <- df_merged[which(df_merged$Subespecie %in% c("Corynebacterium kroppenstedtii DSM 44385", "Cutibacterium acnes subsp. defendens", "Cutibacterium acnes subsp. defendens ATCC 11828")),]
  #df_merged <- df_merged[which(df_merged$Subespecie %in% rownames(Informe_sign)),]
  #df_merged <- df_merged[which(df_merged$Subespecie %in% subespecies_correlacion),]

  edad <- MetadataB$Edad[which(MetadataB$ID == id)]
  rango_real <- MetadataB$`Rango etario`[which(MetadataB$ID == id)]

  # Digo si es alto/Eq/bajo segun su edad
  df_rango_real <- df_merged[which(df_merged$Rango == rango_real),]

  library(dplyr)
  df_rango_real[,2] <- as.numeric(df_rango_real[,2] )
  df_rango_real <- df_rango_real %>%
    mutate(Resultados = ifelse(df_rango_real[,2] > Q1 & df_rango_real[,2] < Q3 | df_rango_real[,2] == Q1 | df_rango_real[,2] == Q3, "Equilibrado",
                               ifelse(df_rango_real[,2] < Q1, "Bajo", "Alto")))


  t <- as.data.frame(table(df_rango_real$Resultados))
  t_max <- t$Var1[which(t$Freq == max(t$Freq))]

  if(any(t_max == "Equilibrado")) {
    rango_final <- rango_real

  } else {

    df_cory <- df_merged[which(df_merged$Subespecie == "Corynebacterium kroppenstedtii DSM 44385"),]
    if(df_cory[1,2] > df_cory$Q3[which(df_cory$Rango == ">55")]) {
      rango_final = ">55"

    #} else if (df_cory[1,2] < df_cory$Q1[which(df_cory$Rango == "18-35")]) {
    #  rango_final = "18-35"

    } else {
      rangos_no_reales <- unique(df_merged$Rango)[-which(unique(df_merged$Rango) == rango_real)]
      df_1 <- df_merged[which(df_merged$Rango == rangos_no_reales[1]),]
      df_1 <- df_1 %>%
        mutate(Resultados = ifelse(df_1[,2] >= Q1 & df_1[,2] <= Q3, "Equilibrado",
                                   ifelse(df_1[,2] < Q1, "Bajo", "Alto")))
      t_d <- as.data.frame(table(df_1$Resultados))

      df_2 <- df_merged[which(df_merged$Rango == rangos_no_reales[2]),]
      df_2 <- df_2 %>%
        mutate(Resultados = ifelse(df_2[,2] >= Q1 & df_2[,2] <= Q3, "Equilibrado",
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

      if(t_d$Freq[which(t_d$Var1 == "Equilibrado")] > t_d2$Freq[which(t_d2$Var1 == "Equilibrado")] ) {
        rango_final <- rangos_no_reales[1]
      } else if (t_d$Freq[which(t_d$Var1 == "Equilibrado")] < t_d2$Freq[which(t_d2$Var1 == "Equilibrado")]) {
        rango_final <- rangos_no_reales[2]
      } else if (t_d$Freq[which(t_d$Var1 == "Equilibrado")] == 0 & t_d2$Freq[which(t_d2$Var1 == "Equilibrado")] == 0) {
        rango_final <- rango_real

      } else { #cuando tienen la misma cantidad de equilibrados:
        if(rango_real == ">55") {
          rango_final <- "35-55"
        } else if(t_d$Freq[which(t_d$Var1 == "Equilibrado")] < t$Freq[which(t$Var1 == "Equilibrado")]) {
          rango_final <- rango_real
        } else if(t_d$Freq[which(t_d$Var1 == "Equilibrado")] > t$Freq[which(t$Var1 == "Equilibrado")]) {
          rango_final <- rango_real

        }
      }

    }

  }

  colnames(df_rango_real)[2] <- gsub("KRAKEN","AR", colnames(df_rango_real)[2])
  colnames(df_rango_real)[6] <- "Min"
  colnames(df_rango_real)[7] <- "Max"
  colnames(df_rango_real)[2] <- "AR"
  #df_rango_real$Subespecie <- c("Biomarcador 1", "Biomarcador 2", "Biomarcador 3", "Biomarcador 4", "Biomarcador 5", "Biomarcador 6",
  #                              "Biomarcador 7", "Biomarcador 8", "Biomarcador 9", "Biomarcador 10", "Biomarcador 11", "Biomarcador 12")

  df_rango_real$Subespecie <- factor(df_rango_real$Subespecie,
                                     levels = (unique(df_rango_real$Subespecie)))

  plot1 <- ggplot(df_rango_real, aes(x = "", y = AR, fill = Subespecie)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y") +  # Para hacer el gráfico circular
    theme_void() +  # Elimina los ejes para un gráfico de pastel
    labs(title = sprintf("Biomarcadores Rango Etario - Muestra %s", id)) +
    theme(
      legend.title = element_blank(),  # Elimina el título de la leyenda
      plot.title = element_text(size = 10, margin = margin(b = 20)),  # Añade un margen inferior al título
      legend.text = element_text(size = 6),  # Cambia el tamaño de la letra de la leyenda
      legend.key.size = unit(0.5, "cm"),
      plot.margin = margin(t = 1, r = 1, b = 1, l = 1)  # Ajusta los márgenes alrededor del gráfico
    )


  df_1 <- df_merged[which(df_merged$Rango == "18-35"),]
  df_1 <- df_1 %>%
    mutate(Resultados = ifelse(df_1[,2] > Q1 & df_1[,2] < Q3, "Equilibrado",
                               ifelse(df_1[,2] < Q1, "Bajo", "Alto")))
  #df_1$Subespecie <- c("Biomarcador 1", "Biomarcador 2", "Biomarcador 3", "Biomarcador 4", "Biomarcador 5", "Biomarcador 6",
  #                              "Biomarcador 7", "Biomarcador 8", "Biomarcador 9", "Biomarcador 10", "Biomarcador 11", "Biomarcador 12")

  plot_rango1 <- ggplot(df_1, aes(x = "", y = Media, fill = Subespecie)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y") +  # Para hacer el gráfico circular
    theme_void() +  # Elimina los ejes para un gráfico de pastel
    labs(title = "Biomarcadores Equilibrados Rango 18-35") +
    theme(
      legend.title = element_blank(),  # Elimina el título de la leyenda
      plot.title = element_text(size = 7, margin = margin(b = 15)),  # Añade un margen inferior al título
      legend.text = element_text(size = 6),  # Cambia el tamaño de la letra de la leyenda
      legend.key.size = unit(0.3, "cm"),
      plot.margin = margin(t = 1, r = 1, b = 1, l = 1)  # Ajusta los márgenes alrededor del gráfico
    )

  df_1 <- df_merged[which(df_merged$Rango == "35-55"),]
  df_1 <- df_1 %>%
    mutate(Resultados = ifelse(df_1[,2] > Q1 & df_1[,2] < Q3, "Equilibrado",
                               ifelse(df_1[,2] < Q1, "Bajo", "Alto")))
  #df_1$Subespecie <- c("Biomarcador 1", "Biomarcador 2", "Biomarcador 3", "Biomarcador 4", "Biomarcador 5", "Biomarcador 6",
  #                     "Biomarcador 7", "Biomarcador 8", "Biomarcador 9", "Biomarcador 10", "Biomarcador 11", "Biomarcador 12")

  plot_rango2 <- ggplot(df_1, aes(x = "", y = Media, fill = Subespecie)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y") +  # Para hacer el gráfico circular
    theme_void() +  # Elimina los ejes para un gráfico de pastel
    labs(title = "Biomarcadores Equilibrados Rango 35-55") +
    theme(
      legend.title = element_blank(),  # Elimina el título de la leyenda
      plot.title = element_text(size = 7, margin = margin(b = 15)),  # Añade un margen inferior al título
      legend.text = element_text(size = 6),  # Cambia el tamaño de la letra de la leyenda
      legend.key.size = unit(0.3, "cm"),
      plot.margin = margin(t = 1, r = 1, b =1, l = 1)  # Ajusta los márgenes alrededor del gráfico
    )

  df_1 <- df_merged[which(df_merged$Rango == ">55"),]
  df_1 <- df_1 %>%
    mutate(Resultados = ifelse(df_1[,2] > Q1 & df_1[,2] < Q3, "Equilibrado",
                               ifelse(df_1[,2] < Q1, "Bajo", "Alto")))
  #df_1$Subespecie <- c("Biomarcador 1", "Biomarcador 2", "Biomarcador 3", "Biomarcador 4", "Biomarcador 5", "Biomarcador 6",
  #                     "Biomarcador 7", "Biomarcador 8", "Biomarcador 9", "Biomarcador 10", "Biomarcador 11", "Biomarcador 12")

  plot_rango3 <- ggplot(df_1, aes(x = "", y = Media, fill = Subespecie)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y") +  # Para hacer el gráfico circular
    theme_void() +  # Elimina los ejes para un gráfico de pastel
    labs(title = "Biomarcadores Equilibrados Rango >55") +
    theme(
      legend.title = element_blank(),  # Elimina el título de la leyenda
      plot.title = element_text(size = 7, margin = margin(b = 15)),  # Añade un margen inferior al título
      legend.text = element_text(size = 6),  # Cambia el tamaño de la letra de la leyenda
      legend.key.size = unit(0.3, "cm"),
      plot.margin = margin(t = 1, r = 1, b = 1, l = 1)  # Ajusta los márgenes alrededor del gráfico
    )

  library(gridExtra)
  combined_plot <- grid.arrange(
    plot1,                   # Primera fila, con un gráfico
    plot_rango1, plot_rango2, plot_rango3,  # Segunda fila, con tres gráficos
    ncol = 3,                 # Para la segunda fila, especificar tres columnas
    layout_matrix = rbind(c(1, 1, 1),      # La primera fila tiene un gráfico que ocupa 3 columnas
                          c(2, 3, 4)),      # La segunda fila tiene tres gráficos, cada uno en su propia columna
    heights = c(2, 2)
  )


  df_rango_real <- df_rango_real[,-c(3,4,5)]

  write.xlsx(df_rango_real, file = sprintf("%s/Tabla_SubEspecies_RangoEtario.xlsx", patient_dir))


  #Dentro del rango etario ver si es para arriba o para abajo ----------------------------------
  df_predicho <- df_merged[which(df_merged$Rango == rango_final),]
  str(df_predicho)
  df_predicho$DifQ1 <- abs(df_predicho[,2] - df_predicho$Q1)
  df_predicho$DifQ3 <- abs(df_predicho[,2] - df_predicho$Q3)

  df_predicho$SubRango <- "Viejo"
  f=1
  for(f in 1:nrow(df_predicho)) {
    print(f)
    if(df_predicho$DifQ1[f] < df_predicho$DifQ3[f] ) { #Si la AR del cory esta mas cerca del Q1: joven
      if( df_predicho$Subespecie[f] == "Corynebacterium kroppenstedtii DSM 44385") {
        df_predicho$SubRango[f] <- "Joven"
      }
    } else { #Si està mas cerca del Q3 y no es cory: joven
      if( df_predicho$Subespecie[f] != "Corynebacterium kroppenstedtii DSM 44385") {
        df_predicho$SubRango[f] <- "Joven"
      }

    }
  }

  t_predicho <- as.data.frame(table(df_predicho$SubRango))
  if(length(nchar(t_predicho$Freq[which(t_predicho$Var1 == "Joven")])) == 0 ) {
    eq <- data.frame( "Var1" = "Joven", "Freq" = 0)
    t_predicho <- rbind(t_predicho, eq)
  }
  if(length(nchar(t_predicho$Freq[which(t_predicho$Var1 == "Viejo")])) == 0 ) {
    eq <- data.frame( "Var1" = "Viejo", "Freq" = 0)
    t_predicho <- rbind(t_predicho, eq)
  }

  if(t_predicho$Freq[which(t_predicho$Var1 == "Joven")] > t_predicho$Freq[which(t_predicho$Var1 == "Viejo")]){
    subrango <- "Dentro del rango etario dado, su piel se encuentra dentro de la población más joven."
    sub_rango <- "inferior"
  } else {
    subrango <- "Dentro del rango etario dado, su piel se encuentra dentro de la población más adulta."
    sub_rango <- "superior"
  }

  if(rango_final == ">55") {
    if(sub_rango == "superior") {
      resultado <- "mayor a 70"
    } else {
      resultado <- "entre 55 y 70"
    }
  } else if (rango_final == "35-55") {
    if(sub_rango == "superior") {
      resultado <- "entre 45 y 55"
    } else {
      resultado <- "entre 35 y 45"
    }
  } else if(rango_final == "18-35") {
    if(sub_rango == "superior") {
      resultado <- "entre 26 y 35"
    } else {
      resultado <- "entre 18 y 25"
    }
  }

  #AJUSTE!!!!!

  if(edad < 35 ) {
    if( edad > 26) {
      if(resultado == "entre 45 y 55" | resultado ==  "entre 55 y 70" | resultado == "mayor a 70" ){
        resultado <- "entre 35 y 45"
      }
    } else {
      if(resultado == "entre 35 y 45" | resultado == "entre 45 y 55" | resultado ==  "entre 55 y 70" | resultado == "mayor a 70" ){
        resultado <- "entre 26 y 35"
      }
    }
  }

  if(edad >= 35 & edad < 56) {
    if( edad > 45) {
      if( resultado == "mayor a 70" ) {
        resultado <- "entre 55 y 70"
      } else if(resultado == "entre 26 y 35" |  resultado == "entre 18 y 25") {
        resultado <- "entre 35 y 45"
      }
    } else {
      if(resultado ==  "entre 55 y 70" | resultado == "mayor a 70" ) {
        resultado <- "entre 45 y 55"
      } else if(resultado == "entre 18 y 25") {
        resultado <- "entre 26 y 35"
      }
    }
  }

  if(edad >56) {
    if( edad > 70) {
      if (resultado == "entre 45 y 55" | resultado == "entre 35 y 45" | resultado == "entre 26 y 35" |  resultado == "entre 18 y 25") {
        resultado <- "entre 55 y 70"
      }
    } else {
      if(resultado ==  "entre 26 y 35" | resultado == "entre 18 y 25" | resultado == "entre 35 y 45" ) {
        resultado <- "entre 45 y 55"
      }
    }
  }
  ################################

  #Mensaje segun como predice:
  mensaje <- sprintf("El color celeste indica que tu microbioma de piel corresponde al grupo de edad %s años.", resultado)

  numeros <- as.numeric(str_extract_all(resultado, "\\d+")[[1]])
  rango_min <- numeros[1]
  rango_max <- numeros[2]
  if(resultado == "mayor a 70") {
    rango_min <- 70
    rango_max <- 85
  }


  # Creamos el gráfico
  library(ggplot2)

  # Definir colores más suaves y estéticos
  color_barra_negra <- "#4D4D4D"  # Gris oscuro para la barra general
  color_barra_rosa <- "#69B3E7"  # Azul pastel para la barra del rango
  color_barra_azul <- "#FF6F91"  # Rosa suave para el subrango

  dev.new(width = 8, height = 0.5)
  barra_plot <- ggplot() +
    # Barra negra general que va de 18 a 90
    geom_segment(aes(x = 18, xend = 90, y = 1, yend = 1),
                 color = color_barra_negra, size = 3, alpha = 0.8) +
    geom_segment(aes(x = rango_min, xend = rango_max, y = 1, yend = 1),
                 color = color_barra_rosa, size = 6, lineend = "round", alpha = 0.9) +
    # Ajustes estéticos
    scale_x_continuous(limits = c(18, 90), breaks = seq(20, 90, by = 5)) +
    scale_y_continuous(limits = c(0.9999, 1.0001)) +
    theme_minimal(base_size = 15) +  # Ajustar el tamaño base para que las fuentes se vean mejor
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
      axis.title.x = element_text(size = 10, face = "bold", color = "#4D4D4D")
      # Márgenes de la gráfica
      #plot.margin = margin(-10, -50,-10, -50)
    ) +
    # Etiqueta del eje x
    labs(x = "Edad")
    #+
    # Agregar una barra de título superior
    #annotate("text", x = 54, y = 1.05, label = "Rango etario indicado según microbioma", size = 6, color = "#4D4D4D", fontface = "bold")

  print(barra_plot)
  ggsave(sprintf("%s/rango_etario_barra.png", patient_dir), plot = barra_plot, width = 6, height = 1.5, dpi = 300)
  dev.off()

  if(rango_min > edad) {
    res <- "un poco más envejecida"
  } else if(rango_max < edad) {
    res <- "rejuvenecida"
  } else {
    res <- "en óptimas condiciones según su rango etario"
  }

  dev_edad <- sprintf("Además, el análisis de las bacterias indica que tu piel se encuentra %s", res)

  return(list(df_rango_real, resultado, plot1, combined_plot, mensaje, barra_plot, dev_edad))

}
