QCcontrol(patients_dir = "~/Daniela/Biota/Muestras/73m", de_host = "Bowtie")
RunKRAKEN(patients_dir = "~/Daniela/Biota/Muestras/73m", de_host = "Bowtie")
#RunKRAKEN(patients_dir = "~/Daniela/Biota/Muestras/73m/169", de_host = "Bowtie")
setwd("~/Daniela/Biota/Muestras/73m/169")
library(PipelineBiota)
otus <- generateOTUsTableGrupal(patients_dir = "~/Daniela/Biota/Muestras/73m", source = "KRAKEN", conEukaryota = TRUE, de_host = "Bowtie")
nivel <- "Species"
str(otus)
library(dplyr)
out <- group_TaxonomicLevels(patients_dir = "~/Daniela/Biota/Muestras/73m", tabla_otus = otus, source = "KRAKEN", de_host = "Bowtie", conEukaryota =  FALSE)
Conteos_totales <- out[[6]]
AR_totales <- out[[1]]
Conteos_Species_totales <- Conteos_totales[[8]]
AR_Species_totales <- AR_totales[[8]]
colnames(Conteos_Species_totales) <- gsub("_KRAKEN", "", colnames(Conteos_Species_totales))
colnames(AR_Species_totales) <- gsub("_KRAKEN", "", colnames(AR_Species_totales))

AR_Bo_Species_KRAKEN <- read_excel("~/Daniela/Biota/Muestras/73m/AR_Bo__Species_KRAKEN.xlsx")


noenblanco <- Conteos_Species_totales[which(Conteos_Species_totales$blanco == 0),]
sienblanco <- Conteos_Species_totales[which(Conteos_Species_totales$blanco != 0),]

counts_tax <- counts_Tax(patients_dir = "~/Daniela/Biota/Muestras/73m", source = "KRAKEN", de_host = "Bowtie", conEukaryota = TRUE)


library(tidyr)
library(ggplot2)
#Conteos_Species_totales <- Conteos_Species_totales[, -which(colnames(Conteos_Species_totales) %in% c("170", "169", "174", "184-1"))]
Conteos_Species_totales_scale <- as.data.frame(t(scale(t(Conteos_Species_totales[,-1]), center = TRUE, scale = TRUE)))

Conteos_Species_totales_scale <- as.data.frame(t(scale(t(AR_Species_totales[,-1]), center = TRUE, scale = TRUE)))
Conteos_Species_totales_scale <- cbind("Species" = Conteos_Species_totales$Species, Conteos_Species_totales_scale)

library(ggplot2)

# Prepara los datos para el PCA
data_pca <- t(Conteos_Species_totales_scale[,-1])  # Transpone los datos (excluyendo la columna de Species)
colnames(data_pca) <- Conteos_Species_totales_scale$Species  # Asigna nombres de filas como especies

# Realiza el PCA
pca_result <- prcomp(data_pca, center = TRUE, scale. = TRUE)

# Crea un data.frame con los resultados del PCA
pca_df <- data.frame(PC1 = pca_result$x[,1],
                     PC2 = pca_result$x[,2],
                     ID = rownames(data_pca))  # IDs de las muestras (nombres de las columnas originales)

# Agrega una columna para resaltar la muestra "blanco"
pca_df$highlight <- ifelse(pca_df$ID == "blanco", "Blanco", "Otros")

library(ggplot2)
# Graficar el PCA etiquetando todos los puntos y resaltando "blanco"
ggplot(pca_df, aes(x = PC1, y = PC2, color = highlight)) +
  geom_point(size = 2) +
  geom_text(aes(label = ID), nudge_y = 5, size = 3) +  # Etiqueta todos los puntos
  scale_color_manual(values = c("Blanco" = "red", "Otros" = "black")) +
  labs(title = "PCA: Resaltando la muestra 'blanco'",
       x = "PC1",
       y = "PC2") +
  theme_minimal()


# Transformar el dataframe al formato long
library(tidyr)
Conteos_Species_totales_long <- Conteos_Species_totales_scale %>%
  pivot_longer(cols = -Species, names_to = "Paciente", values_to = "Conteos")

# Transformar el dataframe al formato long
Conteos_Species_totales_long <- Conteos_Species_totales %>%
  pivot_longer(cols = -Species, names_to = "Paciente", values_to = "Conteos")

Conteos_Species_totales_long <- Conteos_Species_totales_long[which(Conteos_Species_totales_long$Species %in% c("Homo sapiens", "Cutibacterium acnes")),]

ggplot(Conteos_Species_totales_long, aes(x = Species, y = Conteos)) +
  geom_boxplot() +
  labs(
    #title = "Distribución conteo Cutibacterium acnes",
    title = "Distribución conteo escalado especies",
    x = "",
    y = "CPM escalado"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0)) +
  guides(fill = "none") +
  geom_point(data = subset(Conteos_Species_totales_long, Paciente == "blanco"),
             aes(x = Species, y = Conteos),
             color = "red", size = 2) +
  geom_point(data = subset(Conteos_Species_totales_long, Paciente %in% c("138", "41", "130", "90", "36", "144")),
             aes(x = Species, y = Conteos),
             color = "blue", size = 1)


# Crear el boxplot y destacar al paciente "blanco"
ggplot(Conteos_Species_totales_long, aes(x = Paciente, y = Conteos, fill = Paciente %in% c("blanco", "169", "170", "174"))) +
  geom_boxplot() +
  scale_fill_manual(values = c("FALSE" = "gray", "TRUE" = "red")) +
  labs(
    title = "Distribución de especies por paciente (Scaled)",
    x = "Paciente",
    #y = "Conteos Escalados"
    y = "AR escalado"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0)) +
  guides(fill = "none")

#Analizar 10 species más abundantes en piel:
Conteos_Species_totales <- Conteos_Species_totales[order(Conteos_Species_totales$'1', decreasing = TRUE),]
especies_top <- Conteos_Species_totales$Species[1:10]
Especies_imp_todos <- Conteos_Species_totales[which(Conteos_Species_totales$Species %in% especies_top),]

# Aplicar logaritmo y luego scale por columna
Conteos_Species_totales_log_scale <- Especies_imp_todos

# Transformar el dataframe al formato long
Conteos_Species_totales_long <- Conteos_Species_totales_log_scale %>%
  pivot_longer(cols = -Species, names_to = "Paciente", values_to = "Conteos")

# Crear el boxplot y destacar al paciente "blanco"
ggplot(Conteos_Species_totales_long, aes(x = Paciente, y = Conteos, fill = Paciente == "blanco")) +
  geom_boxplot() +
  scale_fill_manual(values = c("FALSE" = "gray", "TRUE" = "red")) +
  labs(
    title = "Distribución de especies por paciente (Log + Scale a media 0)",
    x = "Paciente",
    y = "Conteos Log-Normalizados"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  guides(fill = "none")

especies_imp <- c("Cutibacterium acnes", "Streptococcus oralis", "Corynebacterium striatum", "Staphylococcus epidermidis", "Micrococcus luteus", "Acinetobacter johnsonii", "Streptococcus mitis", "Staphylococcus aureus", "Cutibacterium granulosum")
Especies_imp_todos <- Conteos_Species_totales[which(Conteos_Species_totales$Species %in% especies_imp),]



# Generar tabla con tamaños de muestra:
QCcontrol(patients_dir = "~/Daniela/Biota/Muestras/73m", de_host = "Bowtie")

patients <- list.dirs(patients_dir, full.names = TRUE, recursive = FALSE)
ids <- list.dirs(patients_dir, full.names = FALSE, recursive = FALSE)
Sizes_patients <- data.frame("ID"= ids, "SizeR1" = NA, "SizeR2" = NA, "DHBoSizeR1" = NA, "DHBoSizeR2"= NA )

for (patient in patients) {
  #patient <- patients[1]
  id <- basename(patient)
  list_files <- list.files(patient, full.names = TRUE)

  fileR1 <- list_files[grepl("_L001_R1_001.fastq.gz", list_files)]
  fileR2 <- list_files[grepl("_L001_R2_001.fastq.gz", list_files)]
  Sizes_patients$SizeR1[which(Sizes_patients$ID == id)] <- file.info(fileR1)$size
  Sizes_patients$SizeR2[which(Sizes_patients$ID == id)] <- file.info(fileR2)$size

  patient <- paste(patient, "/trimmed", sep= "")
  print(patient)
  list_files <- list.files(patient, full.names = TRUE)


  fileR1 <- list_files[grepl(sprintf("%sDH%s_S04_L001_R1_001.fastq", id, de_host_file), list_files)]
  fileR2 <- list_files[grepl(sprintf("%sDH%s_S04_L001_R2_001.fastq", id, de_host_file), list_files)]
  Sizes_patients$DHBoSizeR1[which(Sizes_patients$ID == id)] <- file.info(fileR1)$size
  Sizes_patients$DHBoSizeR2[which(Sizes_patients$ID == id)] <- file.info(fileR2)$size

}

Sizes_patients$PromedioR1R2 <- (Sizes_patients$SizeR1 + Sizes_patients$SizeR2)/2
Sizes_patients$PromedioDH <- (Sizes_patients$DHBoSizeR1 + Sizes_patients$DHBoSizeR2)/2
Sizes_patients$DifRDH <- Sizes_patients$PromedioR1R2 - Sizes_patients$PromedioDH

#Comparar cantidad de lecturas:
Table_Basic_Stats <- plotFastQC_PBSQ(patients_dir = "~/Daniela/Biota/Muestras/73m")

Sizes_patients <- merge(Sizes_patients, Table_Basic_Stats, by = "ID")

str(Sizes_patients)

ggplot(Sizes_patients, aes(x = "", y = TotalSequencesMean)) +
  geom_boxplot() +
  labs(
    title = "Distribución de cantidad de lecturas",
    x = "",
    y = "Sizes"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0)) +
  guides(fill = "none") +
  geom_point(data = subset(Sizes_patients, ID == "blanco"),
             aes(x = "", y = TotalSequencesMean),
             color = "red", size = 3) +
  geom_point(data = subset(Sizes_patients, ID %in% c("36", "90", "41")),
             aes(x = "", y = TotalSequencesMean),
             color = "blue", size = 3)

library(tidyr)
library(ggplot2)

# Convertir el dataframe a formato largo
Sizes_patients_long <- Sizes_patients %>%
  pivot_longer(
    cols = c(PromedioR1R2, PromedioDH, DifRDH, TotalSequencesMean),
    names_to = "Variable",
    values_to = "Valor"
  )

ggplot(Sizes_patients_long, aes(x = Variable, y = Valor)) +
  geom_boxplot() +
  labs(
    title = "Distribución de tamaños de archivos + secuencias totales",
    x = "Variable",
    y = "Valores"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  guides(fill = "none") +
  geom_point(
    data = subset(Sizes_patients_long, ID == "blanco"),
    aes(x = Variable, y = Valor),
    color = "red", size = 3
  ) +
  geom_point(data = subset(Sizes_patients_long, ID %in% c("36", "90", "41")),
             aes(x = Variable, y = Valor),
             color = "blue", size = 2)


# Correr human al blanco

RunHuman(patient_dir = "~/Daniela/Biota/Muestras/73m/1", de_host = "Bowtie")

RunHuman(patient_dir = "~/Daniela/Biota/Muestras/73m/blanco", de_host = "Bowtie")
library(readr)
generatePathwayReport(id = "blanco", de_host = "Bowtie")



#Unificar todos los RsubreadBAM summary:
QCcontrol(patients_dir = "~/Daniela/Biota/Muestras/73m", de_host = "Bowtie")

list_dirs <- list.dirs("~/Daniela/Biota/Muestras/73m", full.names = TRUE, recursive = FALSE)
p <-1
Mapeos <- data.frame("ID" = c(), "PorcentajeMapeo" = c())

for (p in 1:length(list_dirs)) {
  patient_dir <- list_dirs[p]
  id <- basename(patient_dir)
  error_file <- sprintf("%s/trimmed/%s_Bo_Error.txt", patient_dir, id)

  if(file.exists(error_file)) {
    porcentaje_mapeo <- readLines(error_file)
    length(porcentaje_mapeo)
    porcentaje_mapeo <- porcentaje_mapeo[length(porcentaje_mapeo)]
    porcentaje_mapeo <-  sub("^(\\d+\\.\\d+)%.*$", "\\1", porcentaje_mapeo)
    if(length(nchar(porcentaje_mapeo)) != 0) {
      Mapeos_ind <- data.frame("ID" = id, "PorcentajeMapeo" = as.numeric(porcentaje_mapeo))
      Mapeos <- rbind(Mapeos, Mapeos_ind)
    } else {
      Mapeos_ind <- data.frame("ID" = id, "PorcentajeMapeo" = "ERROR")
      Mapeos <- rbind(Mapeos, Mapeos_ind)
    }
  } else {
    #QCcontrol(patients_dir = patient_dir, de_host = "Bowtie")
    message(sprintf("patient %s no tiene error file", id))
  }
}

str(Mapeos)
Mapeos <- Mapeos[-which(is.na(Mapeos$PorcentajeMapeo)),]
Mapeos$PorcentajeMapeo <- as.numeric(Mapeos$PorcentajeMapeo)

ggplot(Mapeos, aes(x = "", y = PorcentajeMapeo)) +
  geom_boxplot() +
  labs(
    title = "% de mapeos",
    x = "",
    y = "%s"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0)) +
  guides(fill = "none") +
  geom_point(data = subset(Mapeos, ID == "blanco"),
             aes(x = "", y = PorcentajeMapeo),
             color = "red", size = 3)



# VOLVER A HACER EL DEHOSTING ---------------------------------------
list_dirs <- list.dirs("~/Daniela/Biota/Muestras/73m", full.names = TRUE, recursive = FALSE)
p <-1

for (p in 1:length(list_dirs)) {
  patient_dir <- list_dirs[p]
  id <- basename(patient_dir)

  dir.create(sprintf("%s/trimmed/DeHostsVIEJOS", patient_dir))

  list_dehost_files <- c("Bo", "Rs", "bwa")

  for( dehost_file in list_dehost_files) {
    viejo_dhbo_r1 <- sprintf("%s/trimmed/%sDH%s_S04_L001_R1_001.fastq.gz", patient_dir, id, dehost_file)
    nuevo_dhbo_r1 <- sprintf("%s/trimmed/DeHostsVIEJOS/VIEJO%sDH%s_S04_L001_R1_001.fastq.gz", patient_dir, id, dehost_file)
    file.rename(from = viejo_dhbo_r1, to = nuevo_dhbo_r1)

    viejo_dhbo_r2 <- sprintf("%s/trimmed/%sDH%s_S04_L001_R2_001.fastq.gz", patient_dir, id, dehost_file)
    nuevo_dhbo_r2 <- sprintf("%s/trimmed/DeHostsVIEJOS/VIEJO%sDH%s_S04_L001_R2_001.fastq.gz", patient_dir, id, dehost_file)
    file.rename(from = viejo_dhbo_r2, to = nuevo_dhbo_r2)

  }

  QCcontrol(patient_dir, de_host = "Bowtie")

}


