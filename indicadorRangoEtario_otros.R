library(readxl)
save(df_rangos, file = "~/Daniela/Biota/df_rangos.RData")
load( "~/Daniela/Biota/df_rangos.RData")

patients_dir = "~/Daniela/Biota/Muestras/73m"
MetadataB <- read_excel("~/Daniela/Biota/Metadata-soloColumnasUsables.xlsx")

MetadataB$RangoEtario2 <- ifelse(MetadataB$Edad>45, ">45", "<45")
MetadataB$RangoEtario3 <- ifelse(MetadataB$Edad>50, ">50", "<50")
MetadataB$RangoEtario4 <- ifelse(MetadataB$Edad>55, ">55", "<55")
MetadataB$RangoJovenes <- ifelse(MetadataB$Edad>35, ">35", "<35")

indicadorRangoEtario(id = "125")

# Analizar performance:
list_ids <- list.dirs("~/Daniela/Biota/Muestras/73m", full.names = F, recursive = F)
MetadataB$ID[which(MetadataB$ID == "118")] <- "118-1"
MetadataB$ID[which(MetadataB$ID == "184")] <- "184-1"

resultados <- data.frame("ID" = c(), "Real" =c(), "Predicho" = c(), "Edad" = c())
i=1
id = "27"
id= "156"
id = "140"
id = "1"
id = "173"
for (i in 1:length(list_ids)) {
  id <- list_ids[i]


  #if(rango_predicho != ">55") {
  #  rango_predicho <- indicadorRangoEtario_logFC(id=id)
  #  if(rango_predicho != "18-35") {
  #    rango_predicho <- "35-55"
  #  }
  #}

  out <- indicadorRangoEtario(id = id)
  rango_predicho <- out[[2]]
  df_rango_predicho <- out[[1]]
  #rango_predicho <- indicadorRangoEtario_55(id = id)
  #rango_predicho <- indicadorRangoEtario_35(id = id)
  #if(length(rango_predicho)>1) { rango_predicho <- rango_predicho[1]}

  print(id)
  print(rango_predicho)

  rango_real <- MetadataB$`Rango etario`[which(MetadataB$ID == id)]
  #rango_real <- MetadataB$RangoEtario4[which(MetadataB$ID == id)]
  #rango_real <- MetadataB$RangoJovenes[which(MetadataB$ID == id)]

  edad <- MetadataB$Edad[which(MetadataB$ID == id)]
  resultados[i, "ID"] <- id
  resultados[i, "Real"] <- rango_real
  resultados[i, "Predicho"] <- paste(rango_predicho)
  resultados[i, "Edad"] <- edad
}

# Evaluar la matriz de confusión
resultados$Real <- factor(resultados$Real, levels = c("18-35", "35-55", ">55"))
resultados$Predicho <- factor(resultados$Predicho, levels = c("18-35", "35-55", ">55"))

resultados$Real <- factor(resultados$Real)
resultados$Predicho <- factor(resultados$Predicho)

library(caret)
confusion_matrix <- confusionMatrix(resultados$Real, resultados$Predicho, mode = "everything")
print(confusion_matrix)

res_mal <- resultados[which(resultados$Predicho != resultados$Real),]

res_mal$Predicho <- as.factor(res_mal$Predicho)
res_mal$Edad <- as.numeric(res_mal$Edad)
boxplot(Edad ~ Predicho,
        data = res_mal,
        xlab = "Rango Etario Predicho",
        ylab = "Edad REAL",
        main = "Distribución de MAL clasificados",
        col = "lightblue")


id= "108"
id= "123"
id = "29"
indicadorRangoEtario <- function(id) {

  Biomarcadores_RangoEtario <- read_excel("~/Daniela/Biota/Muestras/73m/Biomarcadores_RangoEtario_sinoutliers.xlsx")

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

    } else if (df_cory[1,2] < df_cory$Q1[which(df_cory$Rango == "18-35")]) {
      rango_final = "18-35"

    } else {
      rangos_no_reales <- unique(df_merged$Rango)[-which(unique(df_merged$Rango) == rango_real)]
      df_1 <- df_merged[which(df_merged$Rango == rangos_no_reales[1]),]
      df_1 <- df_1 %>%
        mutate(Resultados = ifelse(df_1[,2] > Q1 & df_1[,2] < Q3, "Equilibrado",
                                   ifelse(df_1[,2] < Q1, "Bajo", "Alto")))
      t_d <- as.data.frame(table(df_1$Resultados))

      df_2 <- df_merged[which(df_merged$Rango == rangos_no_reales[2]),]
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

      if(t_d$Freq[which(t_d$Var1 == "Equilibrado")] > t_d2$Freq[which(t_d2$Var1 == "Equilibrado")] ) {
        rango_final <- rangos_no_reales[1]
      } else if (t_d$Freq[which(t_d$Var1 == "Equilibrado")] < t_d2$Freq[which(t_d2$Var1 == "Equilibrado")]) {
        rango_final <- rangos_no_reales[2]
      } else if (t_d$Freq[which(t_d$Var1 == "Equilibrado")] == 0 & t_d2$Freq[which(t_d2$Var1 == "Equilibrado")] == 0) {
        rango_final <- rango_real
      } else { #cuando tienen la misma cantidad de equilibrados:
        if(rango_real == ">55") {
          rango_final <- "35-55"
        }
      }

    }

  }

  colnames(df_rango_real)[2] <- gsub("KRAKEN","AR", colnames(df_rango_real)[2])
  colnames(df_rango_real)[6] <- "Min"
  colnames(df_rango_real)[7] <- "Max"
  colnames(df_rango_real)[2] <- "AR"
  df_rango_real$Subespecie <- c("Biomarcador 1", "Biomarcador 2", "Biomarcador 3", "Biomarcador 4", "Biomarcador 5", "Biomarcador 6",
                                "Biomarcador 7", "Biomarcador 8", "Biomarcador 9", "Biomarcador 10", "Biomarcador 11", "Biomarcador 12")

  df_rango_real$Subespecie <- factor(df_rango_real$Subespecie,
                                     levels = (unique(df_rango_real$Subespecie)))

  plot1 <- ggplot(df_rango_real, aes(x = "", y = AR, fill = Subespecie)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y") +  # Para hacer el gráfico circular
    theme_void() +  # Elimina los ejes para un gráfico de pastel
    labs(title = sprintf("Biomarcadores Rango Etario - Muestra %s", id)) +
    theme(
      legend.title = element_blank(),  # Elimina el título de la leyenda
      plot.title = element_text(size = 10, margin = margin(b = 10)),  # Añade un margen inferior al título
      legend.text = element_text(size = 6),  # Cambia el tamaño de la letra de la leyenda
      legend.key.size = unit(0.5, "cm"),
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0)  # Ajusta los márgenes alrededor del gráfico
    )


  plot2 <- ggplot(df_rango_real, aes(x = "", y = Media, fill = Subespecie)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y") +  # Para hacer el gráfico circular
    theme_void() +  # Elimina los ejes para un gráfico de pastel
    labs(title = sprintf("Biomarcadores Equilibrados Rango %s", rango_real)) +
    theme(
      legend.title = element_blank(),  # Elimina el título de la leyenda
      plot.title = element_text(size = 10, margin = margin(b = 10)),  # Añade un margen inferior al título
      legend.text = element_text(size = 6),  # Cambia el tamaño de la letra de la leyenda
      legend.key.size = unit(0.5, "cm"),
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0)  # Ajusta los márgenes alrededor del gráfico
    )


  library(gridExtra)
  combined_plot <- grid.arrange(plot1, plot2, ncol = 2)
  png(sprintf("%s/PieChart_BiomarcadoresRangoEtario.png", patient_dir), width = 900, height = 500)  # Define the output file and size
  grid.arrange(plot1, plot2, ncol = 2)  # Arrange the plots side by side
  dev.off()

  df_rango_real <- df_rango_real[,-c(3,4,5)]


  write.xlsx(df_rango_real, file = sprintf("%s/Tabla_SubEspecies_RangoEtario.xlsx", patient_dir))
  return(list(df_rango_real, rango_final, combined_plot))

}


# SOLO <55 -----------------------------
i = 77
id <- list_ids[i]
indicadorRangoEtario_55_35(id)

indicadorRangoEtario_55 <- function(id) {

  Biomarcadores_RangoEtario_55 <- read_excel("~/Daniela/Biota/Resumen_Rangos><55.xlsx")

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
  sub_out <- subespecies_ausar[which(!(subespecies_ausar %in% otus_subespecies$SubSpecies))]
  #AR_SubSpecies_KRAKEN$'1_KRAKEN'[which(AR_SubSpecies_KRAKEN$SubSpecies %in% sub_out)]

  AR_subespecies <- otus_subespecies[,9:10]
  AR_subespecies[,2:ncol(AR_subespecies)] <- prop.table(as.matrix(AR_subespecies[,2:ncol(AR_subespecies)]), margin = 2) * 100
  sum(AR_subespecies[,2:ncol(AR_subespecies)])


  #Extraigo los biomarcadores para ><55 ---------------
  AR_biomarcadores <- AR_subespecies[which(AR_subespecies$SubSpecies %in% Biomarcadores_RangoEtario_55$Subespecie),]
  colnames(AR_biomarcadores)[1] <- "Subespecie"

  df_merged <- merge(AR_biomarcadores, Biomarcadores_RangoEtario_55, by = "Subespecie")

  df_merged$diferencia_media <- abs(df_merged[,2] - df_merged$Media)
  df_merged$diferencia_mediana <- abs(df_merged[,2] - df_merged$Mediana)

  df_resultado <- df_merged %>%
    group_by(Subespecie) %>%
    slice_min(diferencia_mediana) %>%
    ungroup()

  t <- as.data.frame(table(df_resultado$Rango))
  rango_final <- t$Var1[which(t$Freq == max(t$Freq))]

  if(length(rango_final)>1) {
    df_resultado <- df_merged %>%
      group_by(Subespecie) %>%
      slice_min(diferencia_media) %>%
      ungroup()

    t <- as.data.frame(table(df_resultado$Rango))
    rango_final <- t$Var1[which(t$Freq == max(t$Freq))]
  }

  return(rango_final)

}

# solo <>35 -----------------------------------------------------------------------
indicadorRangoEtario_35 <- function(id) {

  Biomarcadores_RangoEtario_35 <- read_excel("~/Daniela/Biota/Resumen_Rangos><35.xlsx")

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
  sub_out <- subespecies_ausar[which(!(subespecies_ausar %in% otus_subespecies$SubSpecies))]
  #AR_SubSpecies_KRAKEN$'1_KRAKEN'[which(AR_SubSpecies_KRAKEN$SubSpecies %in% sub_out)]

  AR_subespecies <- otus_subespecies[,9:10]
  AR_subespecies[,2:ncol(AR_subespecies)] <- prop.table(as.matrix(AR_subespecies[,2:ncol(AR_subespecies)]), margin = 2) * 100
  sum(AR_subespecies[,2:ncol(AR_subespecies)])


  AR_biomarcadores <- AR_subespecies[which(AR_subespecies$SubSpecies %in% Biomarcadores_RangoEtario_35$Subespecie),]
  colnames(AR_biomarcadores)[1] <- "Subespecie"

  df_merged <- merge(AR_biomarcadores, Biomarcadores_RangoEtario_35, by = "Subespecie")

  df_merged$diferencia_media <- abs(df_merged[,2] - df_merged$Media)
  df_merged$diferencia_mediana <- abs(df_merged[,2] - df_merged$Mediana)

  df_resultado <- df_merged %>%
    group_by(Subespecie) %>%
    slice_min(diferencia_mediana) %>%
    ungroup()

  t <- as.data.frame(table(df_resultado$Rango))
  rango_final <- t$Var1[which(t$Freq == max(t$Freq))]

  if(length(rango_final)>1) {
    df_resultado <- df_merged %>%
      group_by(Subespecie) %>%
      slice_min(diferencia_media) %>%
      ungroup()

    t <- as.data.frame(table(df_resultado$Rango))
    rango_final <- t$Var1[which(t$Freq == max(t$Freq))]

  }

  #rango_final <- ifelse(rango_final == ">35", "35-55", "18-35")
  #if(length(rango_final)>1){ rango_final <- rango_final[1]}

  #return(paste("La edad del microbioma del paciente es de: ", rango_final, ". Su edad real es: ", edad))
  return(rango_final)

}


############################################################
# Log FC ----------------------------------------------
#############################################################

library(limma)
df_completo <- as.data.frame(read_excel("~/Daniela/Biota/df_completo_91_138SubSpecies.xlsx"))

# Seleccionar solo las columnas de abundancias
abundancias <- df_rangos[, 1:9]
abundancias <- df_completo[, 2:139]
# Convertir a matriz (necesario para limma)
abundancias_matrix <- as.matrix(abundancias)
rownames(abundancias_matrix) <- df_completo$ID

# Crear un diseño basado en los rangos etarios
design <- model.matrix(~1 + df_completo$`Rango etario`)
#~1 indica que solo estás incluyendo la intersección
colnames(design) <- gsub("df_rangos$Rango","", colnames(design))

# Ajustar el modelo lineal a los datos
abundancias_matrix <- t(abundancias_matrix)

#2. Transform RNA-Seq Data Ready for Linear Modelling
v <- voom(abundancias_matrix, design, plot = F)

#3. Fit linear model for each gene given a series of arrays:
vfit <- lmFit(v, design)
#Empirical Bayes Statistics for Differential Expression:
lfc = 0.5
efit <- treat(vfit, lfc = lfc)

#Para el informe que quieren:
Informe <- as.data.frame(abundancias_matrix)
Informe$p_value <- efit$p.value[,2]
Informe$p_adj <- p.adjust(efit$p.value[,2], method="fdr")
Informe$stat <- efit$t[,2]
Informe$log2FC <- efit$coefficients[,2]
Informe$stdev <- efit$stdev.unscaled[,2]

Informe_sign <- Informe[which(Informe$p_value<0.05),]
rownames(Informe_sign)

M <- Informe_sign[, -((ncol(Informe_sign)-4) : ncol(Informe_sign))]

library(ComplexHeatmap)
#remotes::install_version("rjson", version = "0.2.20", repos = "https://cloud.r-project.org/")
#install.packages('GetoptLong', dependencies = TRUE, force =TRUE)
#BiocManager::install("ComplexHeatmap", force =TRUE)

Ms <- t(scale(t(M)))

#En vez de los nombres de las muestras se pone el nivel del group:
#all(colnames(Ms) == colnames(Informe_sign))


hm_annot <- ComplexHeatmap::Heatmap(
  Ms,
  column_title = "",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  border =1,
  column_names_gp = gpar(fontsize = 10),# Ajustar el tamaño de la fuente de los nombres de las columnas
  row_names_gp = gpar(fontsize = 8),
  show_column_names = TRUE,
  top_annotation = HeatmapAnnotation(
    Rango = df_completo$`Rango etario`
  )
)


hm_split <- ComplexHeatmap::Heatmap(
  Ms,
  column_title = "Expresión solo Genes Diferenciales por subtipo",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  border =1,
  column_names_gp = gpar(fontsize = 10),# Ajustar el tamaño de la fuente de los nombres de las columnas
  row_names_gp = gpar(fontsize = 8),
  column_split =  df_completo$`Rango etario`,
  show_column_names = TRUE,
  top_annotation = HeatmapAnnotation(
    Rango = df_completo$`Rango etario`
  )
)


# corynebacterium vs las demas --------------------------------------------------------
abundancias <- df_rangos[, 1:9]
colnames(abundancias) <- gsub(" ", "_", colnames(abundancias))   # Nombres de las subespecies

subespecie_interes <- colnames(abundancias)[1]

abundancias_interes <- abundancias[, subespecie_interes]
abundancias_otros <- abundancias[, !colnames(abundancias) %in% subespecie_interes]


# Calcular logFC entre el coryne y otras subesp
logFC <- function(abundancia_interes, abundancia_otros) {
  mean_interes <- mean(abundancia_interes, na.rm = TRUE)
  mean_otros <- mean(abundancia_otros, na.rm = TRUE)

  if (mean_otros == 0) {
    return(NA)
  }

  log2(mean_interes / mean_otros)
}


logFC_values <- apply(abundancias_otros, 2, function(col) logFC(abundancias_interes, col))
logFC_df <- data.frame(Subespecie = colnames(abundancias_otros), LogFC = logFC_values)


#Por rango ---------------------------------------
library(dplyr)

logFC_results_list <- list()
rangos <- unique(df_rangos$Rango)

rango <- rangos[1]
for (rango in rangos) {
  # Filtrar los datos por rango etario
  subset_df <- df_rangos %>% filter(Rango == rango)

  abundancias <- subset_df[, 1:9]
  colnames(abundancias) <- gsub(" ", "_", colnames(abundancias))
  abundancias_matrix <- as.matrix(abundancias)
  rownames(abundancias_matrix) <- subset_df$ID

  abundancias_interes <- abundancias_matrix[,"Corynebacterium_kroppenstedtii_DSM_44385" ]

  # Calcular el logFC para cada subespecie comparada con la coryneb
  abundancias_otros <- abundancias_matrix[,colnames(abundancias_matrix) != "Corynebacterium_kroppenstedtii_DSM_44385" ]
  logFC_values <- apply(abundancias_otros, 2, function(col) {
    log2(mean(col) / mean(abundancias_interes))
  })

  logFC_df <- data.frame(Subespecie = colnames(abundancias_otros), LogFC = logFC_values, Rango = rango)

  logFC_results_list[[rango]] <- logFC_df
}

# Combinar todos los resultados en un solo dataframe
logFC_final_df <- do.call(rbind, logFC_results_list)
write.xlsx(logFC_final_df, file = "~/Daniela/Biota/Biomarcadores_Rango_logFC.xlsx")

# Imprimir el dataframe final
print(logFC_final_df)

id="168"
indicadorRangoEtario_logFC <- function(id) {

  logFC_final_df <- as.data.frame(read_excel("~/Daniela/Biota/Biomarcadores_Rango_logFC.xlsx"))
  Biomarcadores_RangoEtario <- read_excel("~/Daniela/Biota/Biomarcadores_RangoEtario.xlsx")

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
  sub_out <- subespecies_ausar[which(!(subespecies_ausar %in% otus_subespecies$SubSpecies))]
  AR_SubSpecies_KRAKEN$'1_KRAKEN'[which(AR_SubSpecies_KRAKEN$SubSpecies %in% sub_out)]

  AR_subespecies <- otus_subespecies[,9:10]
  AR_subespecies[,2:ncol(AR_subespecies)] <- prop.table(as.matrix(AR_subespecies[,2:ncol(AR_subespecies)]), margin = 2) * 100
  sum(AR_subespecies[,2:ncol(AR_subespecies)])

  #Extraigo los biomarcadores ---------------
  colnames(AR_subespecies)[1] <- "Subespecie"
  #AR_subespecies$Subespecie <- gsub(" ", "_", AR_subespecies$Subespecie)
  AR_biomarcadores <- as.data.frame(AR_subespecies[which(AR_subespecies$Subespecie %in% Biomarcadores_RangoEtario$Subespecie),])

  #Calculo el logFC
  AR_biomarcadores[,2] <- as.numeric(AR_biomarcadores[,2])
  valor_referencia <- AR_biomarcadores[,2][AR_biomarcadores$Subespecie == "Corynebacterium kroppenstedtii DSM 44385"]
  AR_biomarcadores$logFC_Paciente <- log2(AR_biomarcadores[,2] / valor_referencia)


  #Veo a qué rango etario se aproxima más ----------------
  AR_biomarcadores$Subespecie <- gsub(" ", "_", AR_biomarcadores$Subespecie)
  df_merged <- merge(AR_biomarcadores, logFC_final_df, by = "Subespecie", all = T)


  # Calculamos la diferencia absoluta entre el valor del paciente y la media de cada rango etario
  df_merged$diferencia <- abs(df_merged$logFC_Paciente - df_merged$LogFC)

  df_resultado <- df_merged %>%
    group_by(Subespecie) %>%
    slice_min(diferencia) %>%
    ungroup()

  t <- as.data.frame(table(df_resultado$Rango))
  rango_final <- t$Var1[which(t$Freq == max(t$Freq))]


  if(length(rango_final) >1) {
   rango_final <- "35-55"
  }

  edad <- MetadataB$Edad[which(MetadataB$ID == id)]


  # Ver los resultados

  #return(paste("La edad del microbioma del paciente es de: ", rango_final, ". Su edad real es: ", edad))
  return(rango_final)

}


