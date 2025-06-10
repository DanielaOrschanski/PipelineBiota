#' @title Indicador Habitos
#' @import readxl
#' @import data.table
#' @import openxlsx
#' @import ggplot2
#' @import gridExtra
#' @import dplyr
#' @description It generates a prediction on smoking, alcohol drinking, making excercise and using cosmetics by analizing the relative abundance of determined subspecies.
#' @param id patient's id
#' @param MetadataB dataframe that contains the metadata information. id must be within Metadata$ID.
#' @return list with one message for each habits resolution.
#' @export

indicadorHabitos <- function(id, MetadataB) {

  #COSMETICOS ------------------------------------------------------------------------

  Biomarcadores_Cosmeticos <- read_excel("~/Daniela/Biota/PipelineBiota/data/Biomarcadores_CuandoMaquilla_sinoutliers.xlsx")

  #pedir id
  #id <- readline(prompt = "Ingrese el ID del paciente que desea analizar: ")
  list_dirs <- list.dirs(patients_dir, recursive = FALSE, full.names = TRUE)
  patient_dir <- path.expand(sprintf("%s/%s", patients_dir, id))

  if(!(patient_dir %in% list_dirs)) {
    return("This id is not on your patients folder")
  }

  otus <-  read_excel(sprintf("%s/trimmed/TablaOTUS_%s_trimmed_KRAKEN.xlsx", patient_dir, id))
  load("~/Daniela/Biota/Muestras/73m/SubEspecies_AUSAR.RData")

  # me quedo con las 138 subespecies y calculo AR
  otus_subespecies <- otus[which(otus$SubSpecies %in% subespecies_ausar),]

  AR_subespecies <- otus_subespecies[,9:10]
  AR_biomarcadores <- AR_subespecies
  colnames(AR_biomarcadores)[1] <- "Subespecie"
  AR_biomarcadores[,2:ncol(AR_biomarcadores)] <- prop.table(as.matrix(AR_biomarcadores[,2:ncol(AR_biomarcadores)]), margin = 2) * 100
  sum(AR_biomarcadores[,2:ncol(AR_biomarcadores)])

  subespecies_maquillaje <-c("Anaerococcus prevotii DSM 20548",
                             "Corynebacterium jeikeium K411",
                             "Finegoldia magna ATCC 53516",
                             "Fusobacterium nucleatum subsp. animalis",
                             "Fusobacterium nucleatum subsp. polymorphum",
                             "Fusobacterium nucleatum subsp. vincentii",
                             "Haemophilus parainfluenzae T3T1",
                             "Lactobacillus crispatus ST1",
                             "Lactococcus piscium MKFS47",
                             "Leptotrichia buccalis C-1013-b",
                             "Prevotella denticola F0289",
                             "Propionibacterium phage Attacne",
                             "Propionibacterium phage Ouroboros",
                             "Rothia dentocariosa ATCC 17931",
                             "Rothia mucilaginosa DY-18",
                             "Selenomonas sputigena ATCC 35185",
                             "Streptococcus equi subsp. zooepidemicus",
                             "Streptococcus oralis subsp. tigurinus",
                             "Streptococcus parasanguinis FW213",
                             "Streptococcus sanguinis SK36",
                             "Streptococcus thermophilus TH1435",
                             "Xanthomonas campestris pv. raphani")

  AR_biomarcadores <- AR_biomarcadores[which(AR_biomarcadores$Subespecie %in% subespecies_maquillaje),]
  AR_biomarcadores[,2:ncol(AR_biomarcadores)] <- prop.table(as.matrix(AR_biomarcadores[,2:ncol(AR_biomarcadores)]), margin = 2) * 100
  sum(AR_biomarcadores[,2:ncol(AR_biomarcadores)])

  df_merged <- merge(AR_biomarcadores, Biomarcadores_Cosmeticos, by = "Subespecie")

  estados <- unique(df_merged$CuándoMaquillaje)
  df_1 <- df_merged[which(df_merged$CuándoMaquillaje == estados[1]),]
  df_1 <- df_1 %>%
    mutate(Resultados = ifelse(df_1[,2] > Q1 & df_1[,2] < Q3, "Equilibrado",
                               ifelse(df_1[,2] < Q1, "Bajo", "Alto")))
  t_d <- as.data.frame(table(df_1$Resultados))

  df_2 <- df_merged[which(df_merged$CuándoMaquillaje == estados[2]),]
  df_2 <- df_2 %>%
    mutate(Resultados = ifelse(df_2[,2] > Q1 & df_2[,2] < Q3, "Equilibrado",
                               ifelse(df_2[,2] < Q1, "Bajo", "Alto")))
  t_d2 <- as.data.frame(table(df_2$Resultados))

  df_3 <- df_merged[which(df_merged$CuándoMaquillaje == estados[3]),]
  df_3 <- df_3 %>%
    mutate(Resultados = ifelse(df_2[,2] > Q1 & df_2[,2] < Q3, "Equilibrado",
                               ifelse(df_2[,2] < Q1, "Bajo", "Alto")))
  t_d3 <- as.data.frame(table(df_3$Resultados))


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

  # Definir las posiciones correspondientes a cada tipo de piel en el eje x
  maquillaje <- estados
  posiciones_maquillaje <- c(1, 2, 3)  # Asignar una posición para cada piel

  # Crear un dataframe con las frecuencias de "Equilibrado" y las posiciones de las pieles
  equilibrado_data <- data.frame(
    tipos_maquillaje = maquillaje,
    freq = c(
      t_d$Freq[which(t_d$Var1 == "Equilibrado")],
      t_d2$Freq[which(t_d2$Var1 == "Equilibrado")],
      t_d3$Freq[which(t_d3$Var1 == "Equilibrado")]
    ),
    posiciones = posiciones_maquillaje
  )

  # Encontrar la frecuencia máxima de "Equilibrados"
  max_freq <- max(equilibrado_data$freq)

  # Verificar si hay un empate entre dos o más tipos de piel
  maquillaje_max <- equilibrado_data$tipos_maquillaje[equilibrado_data$freq == max_freq]
  pos_maquillaje_max <- equilibrado_data$posiciones[equilibrado_data$freq == max_freq]

  sexo <- MetadataB$Sexo[which(MetadataB$ID == id)]
  if("Diariamente" %in% maquillaje_max & sexo == "Femenino") {
    res_maquillaje <- "Desbalanceado"
    mensaje_maquillaje <- "Según los biomarcadores que utilizamos en BiotaLife para definir el efecto del uso de cosméticos, tu microbioma indica que el estado de tu piel se encuentra desbalanceado."

  } else {
    res_maquillaje <- "Balanceado"
    mensaje_maquillaje <- ""
  }

  ################################################################################################
  # ALCOHOL --------------------------------------------------------------------------------

  Biomarcadores_Alcohol <- read_excel("~/Daniela/Biota/PipelineBiota/data/Biomarcadores_Alcohol_sinoutliers.xlsx")

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

  # Definir las posiciones correspondientes a cada tipo de piel en el eje x
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

  if("SI" %in% estados_max) {
    res_alcohol <- "Desbalanceado"
    mensaje_alcohol <- "Según los biomarcadores que utilizamos en BiotaLife para definir el efecto de la ingesta de alcohol, tu microbioma indica que el estado de tu piel se encuentra desbalanceado."

  } else {
    res_alcohol <- "Balanceado"
    mensaje_alcohol <- ""
  }

  ##################################################################################################

  # TABACO --------------------------------------------------------------------------------

  Biomarcadores_Tabaco <- read_excel("~/Daniela/Biota/PipelineBiota/data/Biomarcadores_Tabaco_sinoutliers.xlsx")

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

  df_merged <- merge(AR_biomarcadores, Biomarcadores_Tabaco, by = "Subespecie")
  #Extraigo los biomarcadores ---------------
  # Digo si es alto/Eq/bajo segun su edad

  estados <- unique(df_merged$Tabaco)
  df_1 <- df_merged[which(df_merged$Tabaco == estados[1]),]
  df_1 <- df_1 %>%
    mutate(Resultados = ifelse(df_1[,2] > Q1 & df_1[,2] < Q3, "Equilibrado",
                               ifelse(df_1[,2] < Q1, "Bajo", "Alto")))
  t_d <- as.data.frame(table(df_1$Resultados))

  df_2 <- df_merged[which(df_merged$Tabaco == estados[2]),]
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

  # Definir las posiciones correspondientes a cada tipo de piel en el eje x
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

  if("SI" %in% estados_max) {
    res_tabaco <- "Desbalanceado"
    mensaje_tabaco <- "Según los biomarcadores que utilizamos en BiotaLife para definir el efecto de la ingesta de alcohol, tu microbioma indica que el estado de tu piel se encuentra desbalanceado."

  } else {
    res_tabaco <- "Balanceado"
    mensaje_tabaco <- ""
  }

  ##################################################################################################
  #ACTIVIDAD FÍSICA  ------------------------------------------------------------------------

  Biomarcadores_ActividadFisica <- read_excel("~/Daniela/Biota/PipelineBiota/data/Biomarcadores_ActividadFisica_sinoutliers.xlsx")
  list_dirs <- list.dirs(patients_dir, recursive = FALSE, full.names = TRUE)
  patient_dir <- path.expand(sprintf("%s/%s", patients_dir, id))

  if(!(patient_dir %in% list_dirs)) {
    return("This id is not on your patients folder")
  }

  otus <-  read_excel(sprintf("%s/trimmed/TablaOTUS_%s_trimmed_KRAKEN.xlsx", patient_dir, id))
  load("~/Daniela/Biota/Muestras/73m/Especies_AUSAR.RData")

  otus_especies <- otus[which(otus$Species %in% especies_ausar),]

  AR_especies <- otus_especies[,c(8,10)]
  AR_especies <- AR_especies %>%
    # Agrupar por nivel
    group_by(Species) %>%
    # Sumar los valores de las columnas para cada especie
    summarise(across(everything(), sum, na.rm = TRUE)) %>%
    ungroup()
  length(unique(AR_especies$Species))

  AR_biomarcadores <- AR_especies
  colnames(AR_biomarcadores)[1] <- "Subespecie"
  AR_biomarcadores[,2:ncol(AR_biomarcadores)] <- prop.table(as.matrix(AR_biomarcadores[,2:ncol(AR_biomarcadores)]), margin = 2) * 100
  sum(AR_biomarcadores[,2:ncol(AR_biomarcadores)])

  df_merged <- merge(AR_biomarcadores, Biomarcadores_ActividadFisica, by = "Subespecie")

  estados <- c("Nunca", "1 o 2 veces por semana", "Todos los días")
  df_1 <- df_merged[which(df_merged$ActividadFísica == estados[1]),]
  df_1 <- df_1 %>%
    mutate(Resultados = ifelse(df_1[,2] > Q1 & df_1[,2] < Q3, "Equilibrado",
                               ifelse(df_1[,2] < Q1, "Bajo", "Alto")))
  t_d <- as.data.frame(table(df_1$Resultados))

  df_2 <- df_merged[which(df_merged$ActividadFísica == estados[2]),]
  df_2 <- df_2 %>%
    mutate(Resultados = ifelse(df_2[,2] > Q1 & df_2[,2] < Q3, "Equilibrado",
                               ifelse(df_2[,2] < Q1, "Bajo", "Alto")))
  t_d2 <- as.data.frame(table(df_2$Resultados))

  df_3 <- df_merged[which(df_merged$ActividadFísica == estados[3]),]
  df_3 <- df_3 %>%
    mutate(Resultados = ifelse(df_2[,2] > Q1 & df_2[,2] < Q3, "Equilibrado",
                               ifelse(df_2[,2] < Q1, "Bajo", "Alto")))
  t_d3 <- as.data.frame(table(df_3$Resultados))


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

  # Definir las posiciones correspondientes a cada tipo de piel en el eje x
  actividadfisica <- estados
  posiciones_actividadfisica<- c(1, 2, 3)  # Asignar una posición para cada piel

  # Crear un dataframe con las frecuencias de "Equilibrado" y las posiciones de las pieles
  equilibrado_data <- data.frame(
    actividadfisica = actividadfisica,
    freq = c(
      t_d$Freq[which(t_d$Var1 == "Equilibrado")],
      t_d2$Freq[which(t_d2$Var1 == "Equilibrado")],
      t_d3$Freq[which(t_d3$Var1 == "Equilibrado")]
    ),
    posiciones = posiciones_actividadfisica
  )

  # Encontrar la frecuencia máxima de "Equilibrados"
  max_freq <- max(equilibrado_data$freq)

  # Verificar si hay un empate entre dos o más tipos de piel
  actividadfisica_max <- equilibrado_data$actividadfisica[equilibrado_data$freq == max_freq]
  pos_actividadfisica_max <- equilibrado_data$posiciones[equilibrado_data$freq == max_freq]

  if("Nunca" %in% actividadfisica_max) {
    res_actividadfisica <- "Desbalanceado"
    mensaje_actividadfisica <- "Según los biomarcadores que utilizamos en BiotaLife para definir el efecto de la realización de actividad física, tu microbioma indica que el estado de tu piel se encuentra desbalanceado."

  } else {
    res_actividadfisica <- "Balanceado"
    mensaje_actividadfisica <- ""
  }
  ##############################################################################################################


  # Crear el gráfico
  library(ggplot2)

  # Datos ficticios
  posiciones_piel <- c(1, 2, 3, 4)  # Posiciones de "Alcohol", "Tabaco" y "Cosméticos"
  pieles <- c("Alcohol", "Tabaco", "Cosméticos", "Actividad Física")

  # Definir colores
  color_barra_negra <- "#4D4D4D"  # Gris oscuro para la barra general
  color_barra_roja <- "#E74C3C"   # Rojo para el desbalance

  # Crear el gráfico
  barra_plot <- ggplot() +
    # Barra negra general que representa todas las categorías
    geom_segment(aes(x = 1, xend = 4, y = 1, yend = 1), color = color_barra_negra, size = 3, alpha = 0.8) +

    # Ajustes estéticos
    scale_x_continuous(limits = c(0.8, 4.2), breaks = posiciones_piel, labels = pieles) +
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
    )

  # Barra roja para "Alcohol" si está desbalanceado
  if (res_alcohol == "Desbalanceado") {
    barra_plot <- barra_plot + geom_segment(aes(x = 1 - 0.2, xend = 1 + 0.2, y = 1, yend = 1),
                 color = color_barra_roja, size = 8, lineend = "round", alpha = 0.7)
  }

  # Barra roja para "Tabaco" si está desbalanceado
  if (res_tabaco == "Desbalanceado") {
    barra_plot <- barra_plot + geom_segment(aes(x = 2 - 0.2, xend = 2 + 0.2, y = 1, yend = 1),
                   color = color_barra_roja, size = 8, lineend = "round", alpha = 0.7)
  }

    # Barra roja para "Cosméticos" si está desbalanceado
  if (res_maquillaje == "Desbalanceado") {
    barra_plot <- barra_plot + geom_segment(aes(x = 3 - 0.2, xend = 3 + 0.2, y = 1, yend = 1),
                 color = color_barra_roja, size = 8, lineend = "round", alpha = 0.7)
  }

  if (res_actividadfisica == "Desbalanceado") {
    barra_plot <- barra_plot + geom_segment(aes(x = 4 - 0.2, xend = 4 + 0.2, y = 1, yend = 1),
                                            color = color_barra_roja, size = 8, lineend = "round", alpha = 0.7)
  }

  print(barra_plot)

  ggsave(sprintf("%s/habitos_barra.png", patient_dir), plot = barra_plot, width = 6, height = 1, dpi = 300)

  mensaje_habitos <- "El color rojo destaca los factores que alteran la salud de tu piel. Según los biomarcadores de Biotalife, tu microbioma sugiere que tu piel está "
  mensaje_habitos <- ifelse(any(c(mensaje_alcohol, mensaje_actividadfisica, mensaje_tabaco, mensaje_maquillaje) != ""),
                            paste0(mensaje_habitos, "desequilibrada debido al "),
                            paste0(mensaje_habitos, "equilibrada con respecto a los hábitos y factores externos.")
                    )

  aporte_alcohol <- ifelse(mensaje_alcohol != "",  "consumo de alcohol ", "")
  aporte_tabaco <- ifelse(mensaje_tabaco != "", "consumo de tabaco", "")
  aporte_maquillaje <- ifelse(mensaje_maquillaje != "",  "uso de maquillaje", "")
  aporte_activ <- ifelse(mensaje_actividadfisica != "",   "falta de actividad física", "")

  list_aportes <- list(aporte_alcohol, aporte_tabaco, aporte_maquillaje, aporte_activ)
  aportes <- c()
  for(aporte in list_aportes) {
    if(aporte != "") {
      aportes <- c(aportes, aporte)
    }
  }

  if(length(aportes)>1) {
    mensaje_habitos <- paste0(mensaje_habitos, paste0(aportes[1:(length(aportes)-1)], collapse =", "))
    mensaje_habitos <- paste0(mensaje_habitos, " y ", aportes[length(aportes)], ".")
  } else {
    mensaje_habitos <- paste0(mensaje_habitos, aportes, ".")
  }


  dev_habitos <- sprintf(", está afectada por %s", paste0(aportes, collapse =", "))

  return(list(mensaje_habitos, dev_habitos))

}
