OTUs_82Pacientes_KRAKEN <- read_excel("~/Daniela/Biota/Muestras/73m/OTUs_83Pacientes_KRAKEN.xlsx")

virus_todos <- OTUs_82Pacientes_KRAKEN[which(OTUs_82Pacientes_KRAKEN$Domain == "Viruses"),]

# me quedo con los virus cuyas AR >1 en al menos 1 muestra
virus_todos[,10:ncol(virus_todos)] <- prop.table(as.matrix(virus_todos[,10:ncol(virus_todos)]), margin = 2) * 100
#virus_todos <- virus_todos[rowSums(virus_todos[, -c(1:9)] >= 1) > 0, ]
#virus_todos[,10:ncol(virus_todos)] <- prop.table(as.matrix(virus_todos[,10:ncol(virus_todos)]), margin = 2) * 100
colSums(virus_todos[,10:ncol(virus_todos)])

colnames(virus_todos) <- gsub("_KRAKEN", "", colnames(virus_todos))
write.xlsx(virus_todos, file = "~/Daniela/Biota/Tabla_viruses_AR>1_82p.xlsx")

out <- group_TaxonomicLevels(patients_dir = "/home/daniela/Daniela/Biota/Muestras/73m", virus_todos, source= "KRAKEN", nombre_extra = "Viruses")
list_AR <- out[[1]]
AR_species <- list_AR[[8]]


AR_esp_imp <- out[[2]]
AR_sub_esp_imp <- out[[3]]

#Me quedo con los 6 peligrosos:
AR_species <- as.data.frame(read_excel("~/Daniela/Biota/Muestras/73m/AR__Species_KRAKEN.xlsx"))
nivel = "Species"
Biomarcadores_Especies_abundantes <- as.data.frame(read_excel(sprintf("~/Daniela/Biota/PipelineBiota/data/%sAbundantes_RangoEtario_sinoutliers.xlsx", nivel)))
colnames(Biomarcadores_Especies_abundantes)[1] <- nivel
Biomarcadores_Especies_abundantes[,1] <- gsub("\\.", " ", Biomarcadores_Especies_abundantes[,1])
especies_ausar <- unique(Biomarcadores_Especies_abundantes$Species)
AR_species <- AR_species[which(AR_species$Species %in% especies_ausar),]
AR_species[,10:ncol(AR_species)] <- prop.table(as.matrix(AR_species[,10:ncol(AR_species)]), margin = 2) * 100
colSums(AR_species[,10:ncol(AR_species)])

virus_rojos <- c("Betapapillomavirus 1", "Betapapillomavirus 2", "Betapapillomavirus 3", "Betapapillomavirus 4", "Betapapillomavirus 5", "Human polyomavirus 5", "Merkel cell polyomavirus")
AR_virus_rojos <- AR_species[which(AR_species$Species %in% virus_rojos),]
AR_virus_rojosT <- as.data.frame(t(AR_virus_rojos))
colnames(AR_virus_rojosT) <- AR_virus_rojosT[1,]
AR_virus_rojosT <- AR_virus_rojosT[-1,]
str(AR_virus_rojosT)
AR_virus_rojosT <- as.data.frame(sapply(AR_virus_rojosT, as.numeric))

resumen <- AR_virus_rojosT %>%
  summarise(
    across(1:(ncol(AR_virus_rojosT)),  # Seleccionar las primeras 8 columnas (las subespecies)
           list(Media = ~mean(.),
                Mediana = ~median(.),
                Q1 = ~quantile(., 0.25),
                Q3 = ~quantile(., 0.75)),
           .names = "{col}_{fn}")) %>%  # Generar nombres de columna automáticamente
  pivot_longer(cols = everything(),  # Convertir todas las columnas excepto 'Rango' en filas
               names_to = c("Subespecie", ".value"),  # Separar el nombre de la subespecie y la métrica
               names_sep = "_")

write.xlsx(resumen, file = "~/Daniela/Biota/PipelineBiota/data/Biomarcadores_Virus_sinoutliers.xlsx", overwrite = TRUE)


#me quedo con fagos:
AR_virus_rojos <- AR_species[which(grepl("Lactococcus phage", AR_species$Species) | grepl("Propionibacterium phage", AR_species$Species) | grepl("Staphylococcus phage", AR_species$Species) | grepl("Streptococcus phage", AR_species$Species) ),]
AR_virus_rojos <- AR_virus_rojos[rowSums(AR_virus_rojos[, -1] != 0) > 0, ]
library(tidyr)
library(ggplot2)


colnames(AR_virus_rojos) <- gsub("_KRAKEN", "", colnames(AR_virus_rojos))
# Reformatear los datos a formato largo
AR_virus_long <- AR_virus_rojos %>%
  pivot_longer(cols = -Species, names_to = "Muestra", values_to = "AR")

# Crear el boxplot con etiquetas de outliers y color en muestras específicas
ggplot(AR_virus_long, aes(x = Species, y = AR)) +
  geom_boxplot(outlier.shape = NA) +  # Evitar que ggplot pinte los outliers por defecto
  #geom_jitter(aes(color = "black"), width = 0.2, height = 0) +  # Añadir puntos para todas las muestras
  geom_text(
    data = AR_virus_long %>%
      group_by(Species) %>%
      filter(AR > 1),
    aes(label = Muestra),
    vjust = -0.5, size = 3, color = "black"
  ) +  # Etiquetar los outliers
  geom_point(
    data = AR_virus_long %>% filter(Muestra %in% c( "37", "138", "165", "175", "35")),
    aes(x = Species, y = AR),
    color = "red", size = 3
  ) +  # Colorear las muestras "36" y "138" en rojo
  labs(title = "Distribución de AR por especie con outliers etiquetados", x = "Especies", y = "AR") +
  theme_minimal() +
  theme(legend.position = "none")  # Opcional: Eliminar la leyenda

#Buscar diferenciales -----------------------------------------------------------------------

AR_SubSpecies_KRAKEN <- AR_species

colSums(AR_SubSpecies_KRAKEN[,-1])
AR_subespeciesT <- as.data.frame(t(AR_SubSpecies_KRAKEN))
colnames(AR_subespeciesT) <- AR_subespeciesT[1,]
AR_subespeciesT <- AR_subespeciesT[-1,]
AR_subespeciesT <- cbind("ID" = rownames(AR_subespeciesT), AR_subespeciesT)
colnames(AR_subespeciesT)[1] <- "ID"
AR_subespeciesT$ID <- sub("_KRAKEN", "", AR_subespeciesT$ID)

df_completo <- merge(AR_subespeciesT, MetadataB, by = "ID")
library(readxl)
df_completo$FacilidadBronceado[which(df_completo$FacilidadBronceado == "1")] <- "2"
df_completo$CuándoMaquillaje[which(df_completo$CuándoMaquillaje == "1 o 2 veces por semana")] <- "Diariamente"
df_completo$AplicacionProtectorSolar[which(df_completo$AplicacionProtectorSolar == "1 o 2 veces por semana")] <- "Solo ante exposición en verano"

n <- ncol(df_completo) - ncol(MetadataB)
colnames(df_completo)[n+2]
str(df_completo)
#Hacer factor las categorias y numericas las AR:
rowSums(df_completo[,2:n+1])
df_completo[,(n+2):ncol(df_completo)] <- lapply(df_completo[,(n+2):ncol(df_completo)], as.factor)
df_completo[,2: (n+1)] <- lapply(df_completo[,2: (n+1)], as.numeric)


p_por_categoria <- data.frame("Categoria" = c(), "SubEspecie" = c(), "T-test" = c(), "Wilcoxon"= c(), "Welch" = c(), "P-adj" =c())
c= n+2
j = 2
i=1
for (c in (n+2):ncol(df_completo)) {
  categoria <- colnames(df_completo)[c]
  print(categoria)
  for (j in 2:(n+1)) {
    genero <- colnames(df_completo)[j]
    print(genero)
    df_sin_na <- df_completo[complete.cases(df_completo[, categoria]), ]

    if(categoria == "RangoJovenes") {
      df_sin_na <- df_sin_na[-which(df_sin_na$`Rango etario`== ">55"),]
    }

    if (any(df_sin_na[,categoria] == "Desconocido")) {
      df_sin_na <- df_sin_na[-which(df_sin_na[,categoria] == "Desconocido"),]
    }

    if( categoria == "Maquillaje_Base") {
      df_sin_na <- df_sin_na[-which(df_sin_na$Sexo == "Masculino"),]
    }

    if( length(levels(df_sin_na[, categoria])) > 2 ) {
      k_test <- kruskal.test(df_sin_na[, genero] ~ df_sin_na[, categoria], data = df_sin_na)
      p_v <- k_test$p.value
      p_valor <- p_v
      p_val <- p_v
      p_adj <- p_v
    } else if (length(levels(df_sin_na[, categoria])) == 2 ) {
      w_test <- wilcox.test(df_sin_na[,genero] ~ df_sin_na[,categoria], df_sin_na)
      p_valor <- w_test$p.value
      #t_test <- t.test(df_sin_na[,genero] ~ df_sin_na[,categoria], df_sin_na)
      #p_val <- t_test$p.value
      #welch <- t.test(df_sin_na[,genero] ~ df_sin_na[,categoria], df_sin_na, var.equal = FALSE)
      #p_v <- welch$p.value
      #p_adj <- p.adjust(p_val, method="BH")
    } else {
      p_valor <- 1
    }

    # Calcular la media para cada grupo en la categoría
    mean_values <- aggregate(df_sin_na[, genero], by = list(df_sin_na[, categoria]), FUN = mean)
    colnames(mean_values) <- c("Grupo", "Media")
    grupo_dominante <- mean_values$Grupo[which.max(mean_values$Media)]

    p_por_categoria[i,"Categoria"] <- categoria
    p_por_categoria[i,"SubEspecie"] <- genero
    #p_por_categoria[i,"T-test"] <- p_val
    p_por_categoria[i,"Wilcoxon"] <- p_valor
    #p_por_categoria[i,"Welch"] <- p_v
    #p_por_categoria[i,"P-adj"] <- p_adj
    p_por_categoria[i, "Grupo Dominante"] <- grupo_dominante
    i <- i+1
  }
}

p_significativos <- p_por_categoria[which(p_por_categoria$Wilcoxon < 0.05),]


#Grafico los significativos
#p_significativos <- p_significativos[which(p_significativos$Categoria == "Rango etario"),]
i=1
library(ggplot2)
library(gridExtra)

unique(p_significativos$Categoria)
k=3
k = k+1

#Hace una grilla de graficos por cada categoria ---------------------------
for(cat in unique(p_significativos$Categoria)) {
  cat <- unique(p_significativos$Categoria)[k]
  print(cat)
  p <- p_significativos[which(p_significativos$Categoria == cat),]

  #p <- p[which(p$SubEspecie %in% rownames(Informe_sign)),]
  lista_graficos <- list()
  i=1
  for (i in 1:nrow(p)) {
    categoria <- p[i, "Categoria"]
    filo <- p[i, 2]

    df_sin_na <- df_completo[complete.cases(df_completo[,categoria]), ]
    df_sin_na <- df_sin_na[, order(colnames(df_sin_na))]
    df_sin_na <- df_sin_na[, unique(colnames(df_sin_na))]
    if (any(df_sin_na[,categoria] == "Desconocido")) {
      df_sin_na <- df_sin_na[-which(df_sin_na[,categoria] == "Desconocido"),]
    }

    # Generar gráfico
    if( categoria == "Rango etario") {
      df_sin_na[,categoria] <- factor(df_sin_na[,categoria], levels = c("18-35", "35-55", ">55"))

    }
    if( categoria == "RangoJovenes") {
      df_sin_na[,categoria] <- factor(df_sin_na[,categoria], levels = c("<35", ">35"))

    }
    if( categoria == "AplicacionProtectorSolar") {
      df_sin_na[,categoria] <- factor(df_sin_na[,categoria], levels = c("Nunca", "Solo ante exposición en verano", "1 o 2 veces por semana","Todos los días"))

    }
    if( categoria == "CuándoMaquillaje") {
      df_sin_na[,categoria] <- factor(df_sin_na[,categoria], levels = c("Nunca", "1 o 2 veces por mes", "1 o 2 veces por semana", "Diariamente"))
      df_sin_na <- df_sin_na[-which(df_sin_na$Sexo == "Masculino"),]
    }
    if( categoria == "Maquillaje_Base") {
      df_sin_na <- df_sin_na[-which(df_sin_na$Sexo == "Masculino"),]
    }
    if( categoria == "ActividadFísica") {
      df_sin_na[,categoria] <- factor(df_sin_na[,categoria], levels = c("Nunca", "1 o 2 veces por semana", "Todos los días"))
    }


    g <- ggplot(df_sin_na, aes_string(x = df_sin_na[,categoria], y = df_sin_na[,filo])) +
      geom_boxplot() +
      theme(axis.text.x = element_text(size = 7),
            axis.text.y = element_text(size = 8)) +
      theme(axis.title.x = element_text(size = 7),
            axis.title.y = element_text(size = 8)) +
      labs(y = paste(filo), x = paste(categoria))


    lista_graficos[[i]] <- g
    i <- i+1
  }

  #selec <- c(1,5,8,9,17, 21,24, 26, 35, 36, 34, 38, 41, 49, 50, 53, 54, 57, 58, 60, 63,64, 65, 66, 67, 74, 76)
  #grid.arrange(grobs = lista_graficos[c(29:55)], ncol = 7)
  grid.arrange(grobs = lista_graficos, ncol = 1)  # Ajusta ncol según la cantidad de gráficos
  lista_graficos <- list()
  k<-k+1
}




# MDS
library(MASS)
AR_genus <- AR_virus_rojos

df_numeric <- as.data.frame(sapply(AR_genus[,-1],as.numeric))
df_numeric <- AR_genus[,-1]
rownames(df_numeric) <- AR_genus$Species
df_numeric <- as.data.frame(t(df_numeric))
dim(df_numeric)
str(df_numeric)

# UMAP
library(umap)
umap_result <- umap(as.matrix(df_numeric), n_neighbors = 15, min_dist = 0.1)
umap_coords <- as.data.frame(umap_result$layout)
colnames(umap_coords) <- c("UMAP1", "UMAP2")
umap_coords$ID <-  MetadataB$ID
umap_coords$Sexo <-  MetadataB$Sexo
umap_coords$Edad <-  MetadataB$Edad
umap_coords$Rango <-  MetadataB$`Rango etario`
umap_coords$AntecedenteEnfermedad <-  MetadataB$AntecedenteEnfermedad
umap_coords$FacilidadBroncearse <-  MetadataB$FacilidadBroncearse
umap_coords$ColorCabello <-  MetadataB$ColorCabello
umap_coords$Peso <-  MetadataB$Peso
umap_coords$TratamientoMédico <-  MetadataB$TratamientoMédico
umap_coords$Tabaco <-  MetadataB$Tabaco
umap_coords$Alcohol <-  MetadataB$Alcohol
umap_coords$ActividadFísica <-  MetadataB$ActividadFísica
umap_coords$AplicacionProtectorSolar <-  MetadataB$AplicacionProtectorSolar

ggplot(umap_coords, aes(x = UMAP1, y = UMAP2, color = Rango)) +
  geom_point(size = 3) +
  geom_text(aes(label = ID), vjust = -1, hjust = 0.5) +
  labs(title = "UMAP de Muestras",
       x = "UMAP1",
       y = "UMAP2") +
  theme_minimal()



# PCA
pca_result <- prcomp(df_numeric, scale. = TRUE)
pca_coords <- as.data.frame(pca_result$x[, 1:2])  # Usar las dos primeras componentes principales
colnames(pca_coords) <- c("PC1", "PC2")
pca_coords$ID <-  MetadataB$ID
pca_coords$Sexo <-  MetadataB$Sexo
pca_coords$Edad <-  MetadataB$Edad
pca_coords$Rango <-  MetadataB$`Rango etario`
pca_coords$AntecedenteEnfermedad <-  MetadataB$AntecedenteEnfermedad
pca_coords$FacilidadBroncearse <-  MetadataB$FacilidadBroncearse
pca_coords$ColorCabello <-  MetadataB$ColorCabello
pca_coords$Peso <-  MetadataB$Peso
pca_coords$TratamientoMédico <-  MetadataB$TratamientoMédico
pca_coords$Tabaco <-  MetadataB$Tabaco
pca_coords$Alcohol <-  MetadataB$Alcohol
pca_coords$ActividadFísica <-  MetadataB$ActividadFísica
pca_coords$AplicacionProtectorSolar <-  MetadataB$AplicacionProtectorSolar

ggplot(pca_coords, aes(x = PC1, y = PC2, color = Rango)) +
  geom_point(size = 3) +
  geom_text(aes(label = ID), vjust = -1, hjust = 0.5) +
  labs(title = "PCA",
       x = "PC1",
       y = "PC2") +
  theme_minimal()

ggplot(pca_coords, aes(x = PC1, y = PC2, color = cantidad_alcohol)) +
  geom_point(size = 3) +
  geom_text(aes(label = ID), vjust = -1, hjust = 0.5) +
  labs(title = "PCA",
       x = "PC1",
       y = "PC2") +
  theme_minimal()


# MDS
library(MASS)
mds_result <- isoMDS(dist(df_numeric))
mds_coords <- as.data.frame(mds_result$points)
colnames(mds_coords) <- c("MDS1", "MDS2")
mds_coords$ID <- df_enf$ID


ggplot(mds_coords, aes(x = MDS1, y = MDS2)) +
  geom_text(aes(label = ""), vjust = -1, hjust = 0.5) +
  geom_point(size = 3) +
  labs(title = "MDS de Muestras",
       x = "MDS1",
       y = "MDS2") +
  theme_minimal()

# t-SNE
library(Rtsne)
pca_coords <- Rtsne(as.matrix(df_numeric), perplexity = 8, theta = 0.5)
pca_coords <- as.data.frame(tsne_result$Y)
colnames(pca_coords) <- c("tSNE1", "tSNE2")
pca_coords$ID <-  MetadataB$ID
pca_coords$ID <-  MetadataB$ID
pca_coords$Sexo <-  MetadataB$Sexo
pca_coords$Edad <-  MetadataB$Edad
pca_coords$Rango <-  MetadataB$`Rango etario`
pca_coords$AntecedenteEnfermedad <-  MetadataB$AntecedenteEnfermedad
pca_coords$FacilidadBroncearse <-  MetadataB$FacilidadBroncearse
pca_coords$ColorCabello <-  MetadataB$ColorCabello
pca_coords$Peso <-  MetadataB$Peso
pca_coords$TratamientoMédico <-  MetadataB$TratamientoMédico
pca_coords$Tabaco <-  MetadataB$Tabaco
pca_coords$Alcohol <-  MetadataB$Alcohol
pca_coords$ActividadFísica <-  MetadataB$ActividadFísica
pca_coords$AplicacionProtectorSolar <-  MetadataB$AplicacionProtectorSolar
pca_coords$ActividadFísica <-  MetadataB$ActividadFísica
pca_coords$ActividadFísica <-  MetadataB$ActividadFísica

ggplot(pca_coords, aes(x = tSNE1, y = tSNE2, color = AplicacionProtectorSolar)) +
  geom_point(size = 3) +
  geom_text(aes(label = ID), vjust = -1, hjust = 0.5) +
  labs(title = "t-SNE de Muestras",
       x = "tSNE1",
       y = "tSNE2") +
  theme_minimal()
