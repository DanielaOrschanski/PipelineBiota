#Alcohol SubEspecie ------------------------------------------------------
load("~/Daniela/Biota/p_significativos-138subspecies_83p.RData")
load("~/Daniela/Biota/p_significativos-99species_83p.RData")
load("~/Daniela/Biota/p_significativos-138subspecies_85p.RData")

str(df_completo)

df_completo <- as.data.frame(read_excel("~/Daniela/Biota/df_completo_85_138SubSpecies.xlsx"))
df_completo <- as.data.frame(read_excel("~/Daniela/Biota/df_completo_83_138SubSpecies.xlsx"))
df_completo <- as.data.frame(read_excel("~/Daniela/Biota/df_completo_83_99Species.xlsx"))

BIOTALIFE_SKIN_Respuestas_Actualizado <- as.data.frame(read_excel("~/Daniela/Biota/BIOTALIFE SKIN  (Respuestas).xlsx", sheet = 4))

#BIOTALIFE_SKIN_Respuestas_Actualizado <- as.data.frame(read_excel("~/Daniela/Biota/BIOTALIFE SKIN (Respuestas)- Actualizado.xlsx", sheet=4))
ids_seq <- list.dirs(path =patients_dir, full.names = F, recursive = F)
BIOTALIFE_SKIN_Respuestas_Actualizado$SacadoOutlier <- "-"
BIOTALIFE_SKIN_Respuestas_Actualizado$SacadoOutlier[which(BIOTALIFE_SKIN_Respuestas_Actualizado$ID %in% ids_seq)] <- "NO"
BIOTALIFE_SKIN_Respuestas_Actualizado$EnLaPlataforma <- "NO"
BIOTALIFE_SKIN_Respuestas_Actualizado$EnLaPlataforma[which(BIOTALIFE_SKIN_Respuestas_Actualizado$ID %in% ids_seq)] <- "SI"

ids_seq_out_acne <- c( "217", "25")
ids_seq_out_norm <- c("169", "170", "174", "CM")
ids_seq_out_rango <- c("62", "114", "123", "126", "128", "129", "137", "153")
BIOTALIFE_SKIN_Respuestas_Actualizado$EnLaPlataforma[which(BIOTALIFE_SKIN_Respuestas_Actualizado$ID %in% c(ids_seq, ids_seq_out_acne, ids_seq_out_norm, ids_seq_out_rango))] <- "SI"

MetadataB <- as.data.frame(read_excel("~/Daniela/Biota/PipelineBiota/data/Metadata-soloColumnasUsables.xlsx"))


BIOTALIFE_SKIN_Respuestas_Actualizado$SacadoOutlier[which(BIOTALIFE_SKIN_Respuestas_Actualizado$ID %in% ids_seq_out_acne)] <- "Si, por acné"
BIOTALIFE_SKIN_Respuestas_Actualizado$SacadoOutlier[which(BIOTALIFE_SKIN_Respuestas_Actualizado$ID %in% ids_seq_out_norm)] <- "Si, por normalización"
BIOTALIFE_SKIN_Respuestas_Actualizado$SacadoOutlier[which(BIOTALIFE_SKIN_Respuestas_Actualizado$ID %in% ids_seq_out_rango)] <- "Si, para rango"

write.xlsx(BIOTALIFE_SKIN_Respuestas_Actualizado, file = "~/Daniela/Biota/BIOTALIFE SKIN (Respuestas)- Actualizado.xlsx")

#MetadataExtra <- read_excel("~/Daniela/Biota/BIOTALIFE SKIN- ActualizadoSept(3).xlsx")
#MetadataExtra$ID <- gsub("\\.0$", "", as.character(MetadataExtra$ID))
#MetadataExtra <- MetadataExtra[which(MetadataExtra$ID %in% MetadataB$ID),]
#all(MetadataB$ID == MetadataExtra$ID)
#MetadataExtra <- MetadataExtra[order(match(MetadataExtra$ID, MetadataB$ID)),]
#MetadataB$AlcoholporSemana <- MetadataExtra$AlcoholporSemana

unique(p_significativos$Categoria)
categoria = "Tabaco"
categoria = "Alcohol"
categoria = "Tipodepiel"
categoria = "FacilidadBroncearse"
categoria = "ActividadFísica"
categoria = "CuándoMaquillaje"
categoria = "AfeccionesPiel"

print(categoria)
p <- p_significativos[which(p_significativos$Categoria == categoria),]
p$SubEspecie
i=1
lista_graficos <- list()
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

  #ELIMINO LOS OUTLIERS!!!!!!!!!!!!!!!!
  if( categoria == "FacilidadBroncearse") {
    df_sin_na <- df_sin_na[-which(df_sin_na$ID %in% ids_out),]
  }

  if(categoria == "Tipodepiel") {
    df_sin_na[,categoria] <- factor(df_sin_na[,categoria], levels = c("Piel Seca", "Piel Mixta", "Piel Grasa"))
  }

  if(categoria == "AfeccionesPiel") {
    df_sin_na[,categoria] <- factor(df_sin_na[,categoria], levels = c("Ninguna", "Acné"))
    df_sin_na <- df_sin_na[-which(df_sin_na$`Rango etario` != "18-35"),]
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

grid.arrange(grobs = lista_graficos, ncol = 5)

# SOLO PARA ACNE: --------------------------------------------------
df_completo <- df_completo[-which(df_completo$`Rango etario` != "18-35"),]
table(df_completo$AfeccionesPiel)

categoria = "AfeccionesPiel"
cepas_acne <- c("Cutibacterium acnes subsp. acnes",
                "Cutibacterium acnes subsp. defendens", #este es bueno, protege contra acné
                "Cutibacterium acnes subsp. defendens ATCC 11828", #este es bueno, protege contra acné
                "Cutibacterium acnes SK137",
                "Cutibacterium acnes TypeIA2 P.acn33",
                "Cutibacterium acnes hdn-1",
                "Cutibacterium acnes TypeIA2 P.acn31",
                "Cutibacterium acnes 266",
                "Cutibacterium acnes C1",
                "Cutibacterium acnes TypeIA2 P.acn17",
                "Cutibacterium acnes KPA171202" )

i=1
lista_graficos <- list()
for (i in 1:length(cepas_acne)) {
  filo <- cepas_acne[i]

  df_sin_na <- df_completo[complete.cases(df_completo[,categoria]), ]
  df_sin_na <- df_sin_na[, order(colnames(df_sin_na))]
  df_sin_na <- df_sin_na[, unique(colnames(df_sin_na))]

  df_sin_na[,categoria] <- factor(df_sin_na[,categoria], levels = c("Ninguna", "Acné"))


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

grid.arrange(grobs = lista_graficos, ncol = 6)

#-----------------------------


#Selecciono las que se separan mejor:
#alcohol:
subesp_enf <- p$SubEspecie[c(1,4,5,6,8,9,12,16, 17,18)]#Estas andan bien
subesp_enf <- p$SubEspecie[c(6,9,12,13,14)]
subesp_enf <- p$SubEspecie[c(1,2,3,4,5,7,10,13,14,15)]

#tipo de piel:
  #especies:
subesp_enf <- p$SubEspecie[c(6,8,9,11,12,13,15,16,17,19,21,22)]
  #subespecies:
subesp_enf <- p$SubEspecie[c(4,20, 21,22,23,24)]
subesp_enf <- p$SubEspecie[c(12,13,14,15,16,17,18,19)]# fagos
  #Las que usé finalmente:
#subesp_enf <- p$SubEspecie[c(4,8, 12,14, 17, 19, 20,21,22,23)]
subesp_enf <- p$SubEspecie[c(4,8, 12,14, 15, 17, 19, 20,21,23)] #esto es una prueba: elimine epidermidis PM221
subesp_enf <- p$SubEspecie[c(4,8, 20,21,23)] #esto es una prueba: elimine fagos y epidermidis PM221
subesp_enf <- p$SubEspecie[c(4,8, 12,14, 15, 17, 19, 20,21,22,23)]

#Fototipo:
subesp_enf <- p$SubEspecie[c(4,5,6,10,11)]
  #Fototipo sin outliers:
subesp_enf <- p$SubEspecie[c(1,7,16,18,21,22,23,24,28,29,30,31,32)]

#CuandoMaquillaje y Tabaco:
subesp_enf <- p$SubEspecie

#Actividad Fisica
subesp_enf <- p$SubEspecie

#ACNE:
subesp_enf <- p$SubEspecie[c(1,2,5,6,7,8,9)]

p <- p[which(p$SubEspecie %in% subesp_enf),]
df_enf <- df_completo[, which(colnames(df_completo) %in% subesp_enf)]

df_enf$ID <- df_completo$ID
df_enf$categoria <- as.factor(df_completo[, categoria])
colnames(df_enf)


df_enf <- df_enf[order(match(df_enf$ID, MetadataB$ID)),]
all(df_enf$ID == MetadataB$ID)
#df_enf$cantidad_alcohol <- MetadataB$AlcoholporSemana
#df_enf$cantidad_alcohol[which(is.na(df_enf$cantidad_alcohol))] <- 0

df_enf$RangoEtario <- MetadataB$`Rango etario`[which(MetadataB$ID %in% df_enf$ID)]
df_enf$Sexo <- MetadataB$Sexo[which(MetadataB$ID %in% df_enf$ID)]

table(df_enf$categoria)

#Elimino outliers: Alcohol ----------------------------------
ids_out <- c("131", "246", "36", "180",
             "13", "115", "133", "175", "144")
#Elimino outliers: fototipo
ids_out <- c("185", "45", "36", "33","29", "207", "138", "179")

#Elimino outliers: acné
ids_out <- df_enf$ID[which(df_enf$RangoEtario != "18-35")]


df_enf <- df_enf[-which(df_enf$ID %in% ids_out),]

#cuando maquilla
if( categoria == "CuándoMaquillaje") {
  df_enf$categoria <- factor(df_enf$categoria, levels = c("Nunca", "1 o 2 veces por mes", "Diariamente"))
  df_enf <- df_enf[-which(df_enf$Sexo == "Masculino"),]
}
#---------------------------------------------------------------

if(any(is.na(df_enf$categoria))) {
  df_enf <- df_enf[-which(is.na(df_enf$categoria)),]
}

str(df_enf)
#df_enf[,1:3] <- lapply(df_enf[,1:3], as.numeric)

resumen <- df_enf %>%
  group_by(categoria) %>%  # Agrupar por variable
  summarise(
    across(1:(ncol(df_enf)-4),  # Seleccionar las primeras 8 columnas (las subespecies)
           list(Media = ~mean(.),
                Mediana = ~median(.),
                Q1 = ~quantile(., 0.25),
                Q3 = ~quantile(., 0.75)),
           .names = "{col}_{fn}")) %>%  # Generar nombres de columna automáticamente
  pivot_longer(cols = -categoria,  # Convertir todas las columnas excepto 'Rango' en filas
               names_to = c("Subespecie", ".value"),  # Separar el nombre de la subespecie y la métrica
               names_sep = "_")
#resumen[,3:6] <- round(resumen[,3:6], 2)
colnames(resumen)[which(colnames(resumen) == "categoria")] <- categoria
colnames(df_enf)[which(colnames(df_enf) == "categoria")] <- categoria
table(df_enf[, categoria])

#all(resumen$Subespecie == Biomarcadores_PielGrasa$Subespecie)

#resumen <- resumen[which(resumen$Subespecie == "Cutibacterium granulosum"),]

write.xlsx(resumen, file = "~/Daniela/Biota/PipelineBiota/data/Biomarcadores_TipodePiel_sinoutliers.xlsx", overwrite = TRUE)
write.xlsx(resumen, file = "~/Daniela/Biota/PipelineBiota/data/Biomarcadores_Fototipo_sinoutliers.xlsx", overwrite = TRUE)
write.xlsx(resumen, file = "~/Daniela/Biota/PipelineBiota/data/Biomarcadores_Alcohol_sinoutliers.xlsx", overwrite = TRUE)
write.xlsx(resumen, file = "~/Daniela/Biota/PipelineBiota/data/Biomarcadores_Tabaco_sinoutliers.xlsx", overwrite = TRUE)
write.xlsx(resumen, file = "~/Daniela/Biota/PipelineBiota/data/Biomarcadores_ActividadFisica_sinoutliers.xlsx", overwrite = TRUE)

#Resumen pero con AR repartidas entre las subespecies seleccionadas:
df_enf_AR <- df_enf
str(df_enf_AR)
df_enf_AR[,-((ncol(df_enf)-3):ncol(df_enf))] <- as.data.frame(prop.table(as.matrix(df_enf_AR[,-((ncol(df_enf_AR)-3):ncol(df_enf_AR))]), margin = 1) * 100)
colnames(df_enf_AR)[which(colnames(df_enf_AR) == categoria)] <- "categoria"
rowSums(df_enf_AR[,-((ncol(df_enf)-3):ncol(df_enf))])

if(any(is.na(df_enf_AR[,-((ncol(df_enf)-3):ncol(df_enf))] ))) {
  df_enf_AR[which(is.na(df_enf_AR[,-((ncol(df_enf)-3):ncol(df_enf))] )),-((ncol(df_enf)-3):ncol(df_enf))]  <- 0
}

#Esto solo para cuando maquillaje: ----------------
subesp_enf <- p$SubEspecie[c(2,4,6,7,11)]
p <- p[which(p$SubEspecie %in% subesp_enf),]
df_enf_AR <- df_enf_AR[,which(colnames(df_enf_AR) %in% subesp_enf)]
df_enf_AR$categoria <- as.factor(df_enf[, categoria])
df_enf_AR$Sexo <- df_enf$Sexo
# -----------------------------------------------


resumen <- df_enf_AR %>%
  group_by(categoria) %>%  # Agrupar por variable
  summarise(
    across(1:(ncol(df_enf_AR)-4),  # Seleccionar las primeras 8 columnas (las subespecies)
           list(Media = ~mean(.),
                Mediana = ~median(.),
                Q1 = ~quantile(., 0.25),
                Q3 = ~quantile(., 0.75)),
           .names = "{col}_{fn}")) %>%  # Generar nombres de columna automáticamente
  pivot_longer(cols = -categoria,  # Convertir todas las columnas excepto 'Rango' en filas
               names_to = c("Subespecie", ".value"),  # Separar el nombre de la subespecie y la métrica
               names_sep = "_")
#resumen[,3:6] <- round(resumen[,3:6], 2)

colnames(resumen)[which(colnames(resumen) == "categoria")] <- categoria
colnames(df_enf_AR)[which(colnames(df_enf_AR) == "categoria")] <- categoria
table(df_enf[, categoria])

write.xlsx(resumen, file = "~/Daniela/Biota/PipelineBiota/data/NEW-sinfagos-Biomarcadores_TipodePiel_sinoutliers.xlsx", overwrite = TRUE)
write.xlsx(resumen, file = "~/Daniela/Biota/PipelineBiota/data/Biomarcadores_CuandoMaquilla_sinoutliers.xlsx", overwrite = TRUE)

#boxplots pero con AR de esas 8 subespecies
p <- p[which(p$SubEspecie %in% subesp_enf),]

str(df_enf_AR)
rowSums(df_enf_AR[,1:(ncol(df_enf_AR)-4)])

lista_graficos <- list()
i=1
for (i in 1:nrow(p)) {
  categoria <- p[i, "Categoria"]
  filo <- p[i, 2]

  df_sin_na <- df_enf_AR

  if(categoria == "Tipodepiel") {
    df_sin_na[,categoria] <- factor(df_sin_na[,categoria], levels = c("Piel Seca", "Piel Mixta", "Piel Grasa"))
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
grid.arrange(grobs = lista_graficos, ncol = 4)


#Visualizaciones ------------------------------------------------------
library(umap)
df_numeric <- df_enf[,-c(ncol(df_enf)-5, ncol(df_enf))]
df_numeric <- as.data.frame(prop.table(as.matrix(df_enf[,-((ncol(df_enf)-5):ncol(df_enf))]), margin = 1) * 100)

# UMAP
umap_result <- umap(as.matrix(df_numeric), n_neighbors = 12, min_dist = 0.1)
umap_coords <- as.data.frame(umap_result$layout)
colnames(umap_coords) <- c("UMAP1", "UMAP2")
umap_coords$ID <-  df_enf$ID
umap_coords$Categoria <-  df_enf[, categoria]

#df_enf$cantidad_alcohol <- MetadataB$AlcoholporSemana
#umap_coords$cantidad_alcohol <- df_enf$cantidad_alcohol

ggplot(umap_coords, aes(x = UMAP1, y = UMAP2, color = Categoria)) +
  geom_point(size = 3) +
  geom_text(aes(label = ID), vjust = -1, hjust = 0.5) +
  labs(title = "UMAP de Muestras",
       x = "UMAP1",
       y = "UMAP2") +
  theme_minimal()

ggplot(umap_coords, aes(x = UMAP1, y = UMAP2, color = cantidad_alcohol)) +
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
pca_coords$ID <- df_enf$ID
pca_coords$Categoria <- df_enf[, categoria]
pca_coords$cantidad_alcohol <- df_enf$cantidad_alcohol

ggplot(pca_coords, aes(x = PC1, y = PC2, color = Categoria)) +
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
mds_coords$Categoria <- df_enf[, categoria]
mds_coords$cantidad_alcohol <- df_enf$cantidad_alcohol
mds_coords$rango <- df_enf$RangoEtario
mds_coords$sexo <- df_enf$Sexo


ggplot(mds_coords, aes(x = MDS1, y = MDS2, color = Categoria)) +
  geom_text(aes(label = ID), vjust = -1, hjust = 0.5) +
  geom_point(size = 3) +
  labs(title = "MDS de Muestras",
       x = "MDS1",
       y = "MDS2") +
  theme_minimal()

# t-SNE
library(Rtsne)
tsne_result <- Rtsne(as.matrix(df_numeric), perplexity = 15, theta = 0.5)
tsne_coords <- as.data.frame(tsne_result$Y)
colnames(tsne_coords) <- c("tSNE1", "tSNE2")
tsne_coords$ID <-  df_enf$ID
tsne_coords$Categoria <-  df_enf[, categoria]

ggplot(tsne_coords, aes(x = tSNE1, y = tSNE2, color = Categoria)) +
  geom_point(size = 3) +
  geom_text(aes(label = ID), vjust = -1, hjust = 0.5) +
  labs(title = "t-SNE de Muestras",
       x = "tSNE1",
       y = "tSNE2") +
  theme_minimal()

#Grafico de barras
df_rangosT <- as.data.frame(t(df_enf[,-c((ncol(df_enf)-3): ncol(df_enf))]))
df_sub_rango <- as.data.frame(prop.table(as.matrix(df_rangosT), margin = 2) * 100)
colSums(df_sub_rango)
colnames(df_sub_rango) <- df_enf$ID

df_sub_rango$SubSpecies <- rownames(df_sub_rango)
df_sub_rango$SubSpecies <- factor(df_sub_rango$SubSpecies, levels = df_sub_rango$SubSpecies)
df_sub_rango <- df_sub_rango[order(df_sub_rango[,2], decreasing = TRUE),]

si <- MetadataB$ID[which(MetadataB$Alcohol == "SI")]
no <- MetadataB$ID[which(MetadataB$Alcohol == "NO")]

acne <- MetadataB$ID[which(MetadataB$AfeccionesPiel == "Acné")]
ninguna <- MetadataB$ID[which(MetadataB$AfeccionesPiel == "Ninguna")]

fototipo2 <- MetadataB$ID[which(MetadataB$FacilidadBroncearse == "2")]
fototipo3 <- MetadataB$ID[which(MetadataB$FacilidadBroncearse == "3")]
fototipo4 <- MetadataB$ID[which(MetadataB$FacilidadBroncearse == "4")]

library(tidyr)
AR_patient_long <- df_sub_rango %>%
  pivot_longer(cols = -SubSpecies, names_to = "Muestra", values_to = "Prop")

AR_patient_long <- AR_patient_long %>%
  mutate(Group = case_when(
    Muestra %in% si ~ "SI",
    Muestra %in% no ~ "NO",
    TRUE ~ "Otro"  # Opcional, para manejar casos no incluidos
  ))

AR_patient_long <- AR_patient_long %>%
  mutate(Group = case_when(
    Muestra %in% acne ~ "Acne",
    Muestra %in% ninguna ~ "SinAfeccion",
    TRUE ~ "Otro"  # Opcional, para manejar casos no incluidos
  ))
AR_patient_long <- AR_patient_long %>%
  mutate(Group = case_when(
    Muestra %in% fototipo2 ~ "Fototipo 2",
    Muestra %in% fototipo3 ~ "Fototipo 3",
    Muestra %in% fototipo4 ~ "Fototipo 4",
    TRUE ~ "Otro"  # Opcional, para manejar casos no incluidos
  ))

AR_patient_long$Muestra <- factor(AR_patient_long$Muestra, levels = c(acne, ninguna))

AR_patient_long$Muestra <- factor(AR_patient_long$Muestra, levels = c(si, no))
AR_patient_long$Muestra <- factor(AR_patient_long$Muestra, levels = c(fototipo2, fototipo3, fototipo4))


library(ggplot2)
library(dplyr)
# Asegurarse de que 'Group' esté en el orden deseado
AR_patient_long <- AR_patient_long %>%
  mutate(Group = factor(Group, levels = c("SI", "NO")))

AR_patient_long <- AR_patient_long %>%
  mutate(Group = factor(Group, levels = c("Acne", "SinAfeccion")))

AR_patient_long <- AR_patient_long %>%
  mutate(Group = factor(Group, levels = c("Fototipo 2", "Fototipo 3", "Fototipo 4")))

AR_patient_long$Group
unique(AR_patient_long$SubSpecies)

# Ordenar muestras segun AR de una especie particular (mejor la que mas tiene);
orden_muestras <- AR_patient_long %>%
  filter(SubSpecies == "Bacillus cereus m1293") %>%
  arrange(Group, desc(Prop))

orden_muestras <- AR_patient_long %>%
  filter(SubSpecies == "Aggregatibacter aphrophilus NJ8700") %>%
  arrange(Group, desc(Prop))

orden_muestras <- AR_patient_long %>%
  filter(SubSpecies == "Propionibacterium phage PHL113M01") %>%
  arrange(Group, desc(Prop))

# Agregar un factor para reordenar las muestras según la proporción dentro de cada grupo
AR_patient_long <- AR_patient_long %>%
  mutate(Muestra = factor(Muestra, levels = orden_muestras$Muestra))

ggplot(AR_patient_long, aes(x = Muestra, y = Prop, fill = SubSpecies)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "Muestra", y = "Proporción", fill = "Subespecie", title = "AR por Muestra - 8 SubEspecies") +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.text = element_text(size = 8),  # Tamaño de texto más pequeño para la leyenda
    legend.title = element_text(size = 8),  # Tamaño del título de la leyenda
    legend.key.width = unit(0.5, "cm"),  # Ancho de la clave de la leyenda
    legend.key.height = unit(0.2, "cm"),  # Altura de la clave de la leyenda
    legend.spacing.y = unit(0.2, "cm"),  # Espacio vertical entre ítems de la leyenda
    axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1, size = 8)
  )  +
  guides(fill = guide_legend(ncol = 1)) +
  facet_wrap(~ Group, scales = "free_x", ncol = 3)  # Crear paneles separados para cada grupo


# PREDICTORES ---------------------------------------------------------------

