
library(readr)
id ="115"
patients_dir <- "~/Daniela/Biota/Muestras/73m"
patient_dir <- "~/Daniela/Biota/Muestras/73m/99"


CPM_Vias_82p_SumaValoresVerdes <- read_excel("~/Daniela/Biota/Muestras/SubsetPathways/CPM_Vias_82p - SumaValoresVerdes.xlsx")
CPM_Vias_Clases_82p <- read_excel("~/Daniela/Biota/Muestras/SubsetPathways/Bo/CPM_Vias_Clases_82p.xlsx")
all(CPM_Vias_82p_SumaValoresVerdes$Pathway %in% CPM_Vias_Clases_82p$Via)
colnames(CPM_Vias_82p_SumaValoresVerdes)[1] <- "Via"
all(CPM_Vias_82p_SumaValoresVerdes$Via %in% CPM_Vias_Clases_82p$Via)

CPM_Vias_Clases_Annotacion <- merge(CPM_Vias_Clases_82p, CPM_Vias_82p_SumaValoresVerdes[, c(1:4, 87)], by ="Via", all = TRUE)
CPM_Vias_Clases_Annotacion <- CPM_Vias_Clases_Annotacion[, c(1,2,85,86,87, 88, 3:84)]

Vias_Clases_interesantes <- read_excel("~/Daniela/Biota/Muestras/SubsetPathways/Vias_Clases interesantes.xlsx")

colnames(Vias_Clases_interesantes)[1] <- "ClaseCPM"
colnames(Vias_Clases_interesantes)[4] <- "ClaseAR"

colnames(Vias_Clases_interesantes)[3] <- "DifSign_Clase_CPM"
colnames(Vias_Clases_interesantes)[6] <- "DifSign_Clase_AR"
Vias_Clases_interesantes$ClaseCPM <- gsub('"', '', Vias_Clases_interesantes$ClaseCPM)
Vias_Clases_interesantes$ClaseAR <- gsub('"', '', Vias_Clases_interesantes$ClaseAR)

all(Vias_Clases_interesantes$ClaseCPM %in% CPM_Vias_Clases_Annotacion$Clase)

clasesCPM <- unique(Vias_Clases_interesantes$ClaseCPM)
clasesCPM <- na.omit(clasesCPM)

Vias_Clases_Anotacion <- CPM_Vias_Clases_Annotacion
Vias_Clases_Anotacion$DifSign_Clase_CPM <- NA
Vias_Clases_Anotacion$DifSign_Clase_AR <- NA
c=1
for(c in 1:length(clasesCPM)) {
  clase <- clasesCPM[c]
  print(clase)
  if(clase %in% Vias_Clases_Anotacion$Clase) {
    dif_sign <- Vias_Clases_interesantes$DifSign_Clase_CPM[which(Vias_Clases_interesantes$ClaseCPM == clase)]
    if(length(dif_sign)>1) {
      dif_sign <- paste0(dif_sign, collapse=", ")
    }
    Vias_Clases_Anotacion$DifSign_Clase_CPM[which(Vias_Clases_Anotacion$Clase == clase)] <- dif_sign
  }
}

clasesAR <- unique(Vias_Clases_interesantes$ClaseAR)
clasesAR <- na.omit(clasesAR)
c=1
for(c in 1:length(clasesAR)) {
  clase <- clasesAR[c]
  print(clase)
  if(clase %in% Vias_Clases_Anotacion$Clase) {
    dif_sign <- Vias_Clases_interesantes$DifSign_Clase_AR[which(Vias_Clases_interesantes$ClaseAR == clase)]
    if(length(dif_sign)>1) {
      dif_sign <- paste0(dif_sign, collapse=", ")
    }
    Vias_Clases_Anotacion$DifSign_Clase_AR[which(Vias_Clases_Anotacion$Clase == clase)] <- dif_sign
  }
}

Vias_Clases_Anotacion <- Vias_Clases_Anotacion[, c(1:6, 89, 90, 7:88)]

write.xlsx(Vias_Clases_Anotacion, file = "~/Daniela/Biota/CPM_Vias_Clases_Anotacion_DifSign_82p.xlsx")

#acne:
QCcontrol(patients_dir = "~/Daniela/Biota/Muestras/acne", de_host = "Bowtie")

RunHuman(patient_dir = "~/Daniela/Biota/Muestras/acne/217", de_host= "Bowtie")
generatePathwayReport(id = id, de_host = "Bowtie")

RunHuman(patient_dir = "~/Daniela/Biota/Muestras/acne/25", de_host= "Bowtie")

list_dirs <- list.dirs(path = "~/Daniela/Biota/Muestras/73m", recursive = FALSE, full.names = TRUE)
p <- list_dirs[2]
for (p in list_dirs) {
  print(p)
  id <- basename(p)
  #RunHuman(patient_dir = p, de_host = "Bowtie")
  #generatePathwayReport(id = id, de_host = "Bowtie")

  RunHuman(patient_dir = p, de_host = "BWA")
  generatePathwayReport(id = id, de_host = "BWA")

  #RunHuman(patient_dir = p, de_host = "")
  #generatePathwayReport(id = id, de_host = "")

  #RunHuman(patient_dir = p, de_host = "RSubread")
  #generatePathwayReport(id = id, de_host = "RSubread")

  #de_host_file = "Bo"
  #vias <- as.data.frame(read_excel(sprintf("%s/%s/Vias/%s%s_vias.xlsx", patients_dir, id,id, de_host_file )))
  #vias <- vias[,-3]
  #colnames(vias)[3] <- paste0("Conteo_", id, sep="")
  #colnames(vias)[4] <- paste0("AR_", id, sep="")

}


## Mover todas las cpm tsv a una carpeta con el nombre del de_host_file:
list_tsv <- list.files("/home/daniela/Daniela/Biota/Muestras/SubsetPathways", full.names = TRUE, recursive = FALSE)
list_tsv <- list_tsv[which(grepl("_CPM_pathabundance.tsv", list_tsv))]
tsv <- list_tsv[2]
for (tsv in list_tsv) {
  print(tsv)
  de_host_file <- gsub("^\\d+([A-Za-z]+)_.*", "\\1", basename(tsv))
  id <- gsub("^([0-9]+).*", "\\1",basename(tsv))
  dir.create(sprintf("/home/daniela/Daniela/Biota/Muestras/SubsetPathways/%s", de_host_file))
  cpm_fileN <- sprintf("/home/daniela/Daniela/Biota/Muestras/SubsetPathways/%s/%s%s_CPM_pathabundance.tsv", de_host_file, id, de_host_file)
  file.rename(from = tsv, to = cpm_fileN)
}

#######################################################################################################################


######################################################################################
# HEATMAPS #######################################
#############################################

library(PipelineBiota)
length(list.dirs( "~/Daniela/Biota/Muestras/73m", recursive = FALSE))
#out <- generatePathwaysTable(patients_dir = "~/Daniela/Biota/Muestras/SubsetPathways/Bo", de_host = "Bowtie")
out <- generatePathwaysTable(patients_dir = "~/Daniela/Biota/Muestras/73m", de_host = "Bowtie")
CPM_conacne <- out[[1]]
AR_conacne <- out[[2]]
CPM_con_org <- out[[3]]

CPM_Vias_Clases_Anotacion_DifSign_82p <- read_excel("~/Daniela/Biota/Muestras/SubsetPathways/Bo/CPM_Vias_Clases_Anotacion_DifSign_82p.xlsx")
vias_colageno <- CPM_Vias_Clases_Anotacion_DifSign_82p$Via[which(grepl("Colágeno", CPM_Vias_Clases_Anotacion_DifSign_82p$CategoriaInteres))]
vias_acidosgrasos <- CPM_Vias_Clases_Anotacion_DifSign_82p$Via[which(grepl("Acidos grasos", CPM_Vias_Clases_Anotacion_DifSign_82p$CategoriaInteres))]
vias_cicatrizacion <- CPM_Vias_Clases_Anotacion_DifSign_82p$Via[which(grepl("Cicatrización", CPM_Vias_Clases_Anotacion_DifSign_82p$CategoriaInteres))]
vias_antioxidante <- CPM_Vias_Clases_Anotacion_DifSign_82p$Via[which(grepl("Antioxidante", CPM_Vias_Clases_Anotacion_DifSign_82p$CategoriaInteres))]
vias_acidohialuronico <- CPM_Vias_Clases_Anotacion_DifSign_82p$Via[which(grepl("Acido hialurónico", CPM_Vias_Clases_Anotacion_DifSign_82p$CategoriaInteres))]

# HACER HEATMAP DE VIAS Y BACTERIAS DEL PAPER: https://advanced.onlinelibrary.wiley.com/doi/epdf/10.1002/advs.202300050 ------------------------
colnames(CPM_Vias_Clases_Anotacion_DifSign_82p)[1] <- "Pathway"

CPM_con_org <- merge(CPM_con_org, CPM_Vias_Clases_Anotacion_DifSign_82p[, c("Pathway", "ClasifiaciónSegunPaper")], by = "Pathway")
CPM_con_org <- cbind(CPM_con_org[,1:5], CPM_con_org[,"ClasifiaciónSegunPaper"], CPM_con_org[,6:88])
colnames(CPM_con_org)[6] <- "ClasificacionSegunPaper"

CPM_con_org_sinna <- CPM_con_org

#Completar con info de filos:
tabla_otus <- generateOTUsTableGrupal(patients_dir = "~/Daniela/Biota/Muestras/73m", source = "KRAKEN", conEukaryota = FALSE, de_host = "Bowtie")
otus_unique <- tabla_otus[!duplicated(tabla_otus$Genus), c("Genus", "Phylum")]

colnames(CPM_con_org_sinna)[4] <- "Genus"
CPM_con_org_sinna_filo <- merge(CPM_con_org_sinna, otus_unique, by = "Genus", all.x = TRUE, all.y = TRUE)
colnames(CPM_con_org_sinna_filo)
CPM_con_org_sinna_filo <- cbind("Pathway" = CPM_con_org_sinna_filo$Pathway, CPM_con_org_sinna_filo[,c(3,4,90,5,1,6)], CPM_con_org_sinna_filo[, 7:89])

CPM_con_org_sinna_filo <- CPM_con_org_sinna_filo[-which(is.na(CPM_con_org_sinna_filo$Pathway)),]
any(is.na(CPM_con_org_sinna_filo$Organism))

CPM_sinorg <- CPM_con_org_sinna_filo[which(is.na(CPM_con_org_sinna_filo$Organism)),]

# ME quedo solo con los pathways que tienen clasificacion segun el paper:
CPM_sin_org_solopaper <- CPM_sinorg[-which(is.na(CPM_sinorg$ClasificacionSegunPaper)),]
CPM_con_org_solopaper <- CPM_con_org_sinna_filo[-which(is.na(CPM_con_org_sinna_filo$ClasificacionSegunPaper)),]

#Me quedo con todas las vías: ----
CPM_sin_org_solopaper <- CPM_sinorg
CPM_con_org_solopaper <- CPM_con_org_sinna_filo
CPM_con_org_solopaper$ClasificacionSegunPaper <- ifelse(is.na(CPM_con_org_solopaper$ClasificacionSegunPaper), "NoClasificado", CPM_con_org_solopaper$ClasificacionSegunPaper )
unique(CPM_con_org_solopaper$ClasificacionSegunPaper)
# -----------------------------------------

length(unique(CPM_con_org_solopaper$Pathway)) # 33 vias en total
length(unique(CPM_sin_org_solopaper$Pathway)) # 33 vias en total

CPM_con_org_solopaper <- CPM_con_org_solopaper[-which(is.na(CPM_con_org_solopaper$Organism)),]
CPM_con_org_solopaper <- CPM_con_org_solopaper[-which(CPM_con_org_solopaper$Organism == "unclassified"),]

#AR en vez de CPM:
CPM_con_org_solopaper[, 8:ncol(CPM_con_org_solopaper)] <- prop.table(as.matrix(CPM_con_org_solopaper[, 8:ncol(CPM_con_org_solopaper)]), margin = 2) * 100
colSums(CPM_con_org_solopaper[, 8:ncol(CPM_con_org_solopaper)] )
CPM_con_org_solopaper[, 8:ncol(CPM_con_org_solopaper)][is.na(CPM_con_org_solopaper[, 8:ncol(CPM_con_org_solopaper)])] <- 0

# HEATMAP COMO PAPER: ----------------------------------------------
MetadataB <- as.data.frame(read_excel("~/Daniela/Biota/PipelineBiota/data/Metadata-soloColumnasUsables.xlsx"))
rownames(MetadataB) <- MetadataB$ID
MetadataB <- MetadataB[colnames(CPM_con_org_solopaper), ]

ids_raros <- c("36", "41", "90", "138", "144")
jovenes <- MetadataB$ID[which(MetadataB$`Rango etario` == "18-35")]
medios <- MetadataB$ID[which(MetadataB$`Rango etario` == "35-55")]
adultos <- MetadataB$ID[which(MetadataB$`Rango etario` == ">55")]
ids <- colnames(CPM_con_org_solopaper)[8:ncol(CPM_con_org_solopaper)]
ids_sin_raros <- setdiff(ids, ids_raros)

CPM_con_org_solopaper$Promedio <- rowMeans(CPM_con_org_solopaper[,jovenes])
CPM_con_org_solopaper$Promedio <- rowMeans(CPM_con_org_solopaper[,adultos])
CPM_con_org_solopaper$Promedio <- rowMeans(CPM_con_org_solopaper[,8:ncol(CPM_con_org_solopaper)])
CPM_con_org_solopaper$Promedio <- rowMeans(CPM_con_org_solopaper[,ids_sin_raros])

library(dplyr)
library(ComplexHeatmap)
library(circlize)

CPM_con_org_solopaper_filtered <- CPM_con_org_solopaper %>%
  select(Pathway, ClasificacionSegunPaper, Especies, Phylum, Promedio)

if(any(is.na(CPM_con_org_solopaper_filtered$Especies))) {
  CPM_con_org_solopaper_filtered <- CPM_con_org_solopaper_filtered[-which(is.na(CPM_con_org_solopaper_filtered$Especies)),]
}

if(any(is.na(CPM_con_org_solopaper_filtered$Phylum))) {
  CPM_con_org_solopaper_filtered <- CPM_con_org_solopaper_filtered[-which(is.na(CPM_con_org_solopaper_filtered$Phylum)),]
}

# Species en las filas y los Pathways en las columnas
heatmap_data <- CPM_con_org_solopaper_filtered %>%
  select(Especies, Pathway, Promedio) %>%
  pivot_wider(names_from = Pathway, values_from = Promedio)

heatmap_matrix <- as.matrix(heatmap_data[, -1])
str(heatmap_matrix)
#heatmap_matrix <- apply(heatmap_matrix[,-1], 2, function(x) as.numeric(as.character(x)))
rownames(heatmap_matrix) <- heatmap_data$Especies

# Crear la anotación para las vías:
pathway_annotations <- CPM_con_org_solopaper %>%
  select(Pathway, ClasificacionSegunPaper) %>%
  distinct()  # Asegurarse de que no haya duplicados de Pathway y ClasificacionSegunPaper

column_split <- pathway_annotations$ClasificacionSegunPaper[match(colnames(heatmap_matrix), pathway_annotations$Pathway)]
pathway_colors <- structure(rainbow(length(unique(pathway_annotations$ClasificacionSegunPaper))),
                            names = unique(pathway_annotations$ClasificacionSegunPaper))
column_ha <- HeatmapAnnotation(
  df = data.frame(ClasificacionSegunPaper = factor(column_split, levels = unique(column_split))),
  col = list(ClasificacionSegunPaper = pathway_colors)
)

#Anotación especies
row_split <- CPM_con_org_solopaper$Phylum[match(rownames(heatmap_matrix), CPM_con_org_solopaper$Especies)]
phylum_colors <- structure(rainbow(length(unique(CPM_con_org_solopaper$Phylum))),
                           names = unique(CPM_con_org_solopaper$Phylum))
names(phylum_colors)[4] <- "NA"
row_ha <- rowAnnotation(
  Phylum = factor(row_split, levels = unique(row_split)),
  col = list(Phylum = phylum_colors)
)


#heatmap_matrix_normalized <- apply(heatmap_matrix, 2, function(x) (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
#heatmap_matrix_normalized <- scale(heatmap_matrix)
heatmap_matrix_normalized <- log10(heatmap_matrix + 1)

col_fun <- colorRamp2(c(0, 1), c("white", "purple"))
col_fun <- colorRamp2(c(min(heatmap_matrix_normalized, na.rm = TRUE), max(heatmap_matrix_normalized, na.rm = TRUE)),
                      c("#E6CCFF", "#4B0082")) # Lila clarito → Violeta oscuro

col_fun <- colorRamp2(c(min(heatmap_matrix, na.rm = TRUE), max(heatmap_matrix, na.rm = TRUE)),
                      c("#E6CCFF", "#4B0082"))
# Crear el heatmap
Heatmap(
  #heatmap_matrix,
  heatmap_matrix_normalized,
  name = "Promedio",
  cluster_rows = FALSE,  # No agrupar las filas (Species)
  cluster_columns = FALSE,  # No agrupar las columnas (Pathways)
  row_names_side = "left",  # Colocar los nombres de las filas a la izquierda
  column_names_side = "top",  # Colocar los nombres de las columnas arriba
  row_title = "Species",
  column_title = "Pathways",
  col = col_fun,  # Colores definidos
  na_col = "white",
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (!is.na(heatmap_matrix[i, j])) {  # Evita imprimir NA
      grid.text(
        round(heatmap_matrix[i, j], 2),  # Redondea a 2 decimales
        x, y, gp = gpar(fontsize = 5.5, col = "black") # Ajusta tamaño y color
      )
    }
  },
  row_names_gp = gpar(fontsize = 8),  # Reducir tamaño de las filas
  column_names_gp = gpar(fontsize = 7),  # Reducir tamaño de las columnas
  column_split = column_split,  # Agrupar columnas por Pathways
  row_split = row_split,
  top_annotation = column_ha,  # Mostrar las anotaciones de ClasificacionSegunPaper
  left_annotation = row_ha
)

# HEATMAP POR MUESTRA -----------------------------------------------------------
CPM_usar <- CPM_sin_org_solopaper

#Escalar los valores para que estén entre 0 y 1:
CPM_usar[, 8:ncol(CPM_usar)] <- apply(
  CPM_usar[, 8:ncol(CPM_usar)],
  2,
  function(x) {
    min_x <- min(x, na.rm = TRUE)
    max_x <- max(x, na.rm = TRUE)

    if (max_x == min_x) {
      return(rep(0, length(x)))  # Si todos los valores son iguales, deja todo en 0
    } else {
      return((x - min_x) / (max_x - min_x))  # Escala entre 0 y 1
    }
  }
)
colSums(CPM_usar[,8:ncol(CPM_usar)])


# Extraer solo los valores numéricos (todas las columnas de muestras)
heatmap_data <- as.matrix(CPM_usar[, 8:ncol(CPM_usar)])
rownames(heatmap_data) <- CPM_usar$Pathway
unique(rownames(heatmap_data))

# Asegurar que los nombres de las columnas coincidan con los de heatmap_data
MetadataB <- as.data.frame(read_excel("~/Daniela/Biota/PipelineBiota/data/Metadata-soloColumnasUsables.xlsx"))
rownames(MetadataB) <- MetadataB$ID
MetadataB <- MetadataB[colnames(heatmap_data), ]

# Definir colores para anotaciones
library(ComplexHeatmap)
library(circlize)
library(dplyr)
ann_colors <- list(
  Sexo = c("Femenino" = "red", "Masculino" = "blue"),
  `Rango etario` = c("18-35" = "lightgreen", "35-55" = "orange", ">55" = "violet")
)
# Crear anotación para las muestras
row_ha <- HeatmapAnnotation(df = MetadataB[, c("Sexo", "Rango etario")], col = ann_colors, which = "column")
MetadataB$`Rango etario`

# Crear dataframe con ClasificacionPaper
pathway_annotations <- data.frame(ClasificacionPaper = CPM_usar$ClasificacionSegunPaper)
rownames(pathway_annotations) <- rownames(heatmap_data)

# Definir colores para las categorías de ClasificacionPaper
pathway_colors <- structure(rainbow(length(unique(pathway_annotations$ClasificacionPaper))),
                            names = unique(pathway_annotations$ClasificacionPaper))

# Crear anotación para las vías metabólicas
col_ha <- HeatmapAnnotation(df = pathway_annotations, col = list(ClasificacionPaper = pathway_colors), which = "row")

#Todo sin outliers:
ids_raros <- c("36", "41", "90", "138", "144")
medios <- MetadataB$ID[which(MetadataB$`Rango etario` == "35-55")]
ids_raros <- c(ids_raros, medios)
heatmap_data_sin_outliers <- heatmap_data[, !(colnames(heatmap_data) %in% ids_raros)]
MetadataB_sin_outliers <- MetadataB[colnames(heatmap_data_sin_outliers), ]
row_ha_sin_outliers <- HeatmapAnnotation(df = MetadataB_sin_outliers[, c("Sexo", "Rango etario")],
                            col = ann_colors, which = "column")
pathway_annotations_sin_outliers <- data.frame(ClasificacionPaper = CPM_usar$ClasificacionSegunPaper)
rownames(pathway_annotations_sin_outliers) <- rownames(heatmap_data_sin_outliers)
pathway_colors_sin_outliers <- structure(rainbow(length(unique(pathway_annotations_sin_outliers$ClasificacionPaper))),
                                         names = unique(pathway_annotations_sin_outliers$ClasificacionPaper))

col_ha_sin_outliers <- HeatmapAnnotation(df = pathway_annotations_sin_outliers, col = list(ClasificacionPaper = pathway_colors_sin_outliers), which = "row")

Heatmap(
  #heatmap_data,
  heatmap_data_sin_outliers,
  name = "CPM escalado min-max",  # Nombre de la escala de colores
  cluster_rows = TRUE,  # Agrupar muestras
  cluster_columns = TRUE,  # Agrupar vías metabólicas
  row_names_side = "left",
  column_names_side = "top",
  row_dend_side = "left",
  column_dend_side = "top",
  row_title = "Muestras",
  column_title = "Vías metabólicas",

  #top_annotation = row_ha,  # Anotaciones de muestras
  #left_annotation = col_ha,  # Anotaciones de vías

  #sin outliers:
  top_annotation = row_ha_sin_outliers,
  left_annotation = col_ha_sin_outliers,

  col = colorRamp2(c(0, 0.5, 1), c("white", "yellow", "red")),  # Gradiente de colores
  row_split = pathway_annotations$ClasificacionPaper,
  #column_split = factor(MetadataB$`Rango etario`, levels = c("18-35", "35-55", ">55")),
  column_split = factor(MetadataB_sin_outliers$`Rango etario`, levels = c("18-35", "35-55", ">55")),  # Agrupar las columnas según Rango Etario

  row_names_gp = gpar(fontsize = 8),  # Reducir tamaño de la fuente para los nombres de las filas
  column_names_gp = gpar(fontsize = 8)  # Reducir tamaño de la fuente para los nombres de las columnas
)

#REPLICAR FIG 4 PAPER: -------------------------------------------------------------
#https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-024-01891-0/figures/4
library(PipelineBiota)
length(list.dirs( "~/Daniela/Biota/Muestras/73m", recursive = FALSE))
#out <- generatePathwaysTable(patients_dir = "~/Daniela/Biota/Muestras/SubsetPathways/Bo", de_host = "Bowtie")
out <- generatePathwaysTable(patients_dir = "~/Daniela/Biota/Muestras/73m", de_host = "Bowtie")
CPM_conacne <- out[[1]]
AR_conacne <- out[[2]]
colSums(AR_conacne[,-c(1:3)])

library(readxl)
CPM_Vias_Clases_Anotacion_DifSign_82p <- read_excel("~/Daniela/Biota/Muestras/SubsetPathways/Bo/CPM_Vias_Clases_Anotacion_DifSign_82p.xlsx")
unique(CPM_Vias_Clases_Anotacion_DifSign_82p$CategoriaInteres)
vias_colageno <- CPM_Vias_Clases_Anotacion_DifSign_82p$Via[which(grepl("Colágeno", CPM_Vias_Clases_Anotacion_DifSign_82p$CategoriaInteres))]
vias_acidosgrasos <- CPM_Vias_Clases_Anotacion_DifSign_82p$Via[which(grepl("Acidos grasos", CPM_Vias_Clases_Anotacion_DifSign_82p$CategoriaInteres))]
vias_cicatrizacion <- CPM_Vias_Clases_Anotacion_DifSign_82p$Via[which(grepl("Cicatrización", CPM_Vias_Clases_Anotacion_DifSign_82p$CategoriaInteres))]
vias_antioxidante <- CPM_Vias_Clases_Anotacion_DifSign_82p$Via[which(grepl("Antioxidante", CPM_Vias_Clases_Anotacion_DifSign_82p$CategoriaInteres))]
vias_acidohialuronico <- CPM_Vias_Clases_Anotacion_DifSign_82p$Via[which(grepl("Acido hialurónico", CPM_Vias_Clases_Anotacion_DifSign_82p$CategoriaInteres))]
vias_hidratantes <- CPM_Vias_Clases_Anotacion_DifSign_82p$Via[which(grepl("Factor hidratante natural", CPM_Vias_Clases_Anotacion_DifSign_82p$CategoriaInteres))]
vias_uv <- CPM_Vias_Clases_Anotacion_DifSign_82p$Via[which(grepl("UV protection", CPM_Vias_Clases_Anotacion_DifSign_82p$CategoriaInteres))]
vias_ROS <- CPM_Vias_Clases_Anotacion_DifSign_82p$Via[which(grepl("ROS", CPM_Vias_Clases_Anotacion_DifSign_82p$CategoriaInteres))]
vias_eccema <- CPM_Vias_Clases_Anotacion_DifSign_82p$Via[which(grepl("Eccema", CPM_Vias_Clases_Anotacion_DifSign_82p$CategoriaInteres))]
vias_dermatitis <- CPM_Vias_Clases_Anotacion_DifSign_82p$Via[which(grepl("Dermatitis Atopica", CPM_Vias_Clases_Anotacion_DifSign_82p$CategoriaInteres))]
vias_infeccion <- CPM_Vias_Clases_Anotacion_DifSign_82p$Via[which(grepl("Infeccion", CPM_Vias_Clases_Anotacion_DifSign_82p$CategoriaInteres))]


vias_anti_aging <- CPM_Vias_Clases_Anotacion_DifSign_82p$Via[which(grepl("Anti-Aging", CPM_Vias_Clases_Anotacion_DifSign_82p$ClasifiaciónSegunPaper))]
vias_anti_oxidation <- CPM_Vias_Clases_Anotacion_DifSign_82p$Via[which(grepl("Anti-oxidation", CPM_Vias_Clases_Anotacion_DifSign_82p$ClasifiaciónSegunPaper))]
vias_biosintesis <- CPM_Vias_Clases_Anotacion_DifSign_82p$Via[which(grepl("Biosintesis", CPM_Vias_Clases_Anotacion_DifSign_82p$ClasifiaciónSegunPaper))]
vias_inflamation <- CPM_Vias_Clases_Anotacion_DifSign_82p$Via[which(grepl("Inflamación", CPM_Vias_Clases_Anotacion_DifSign_82p$ClasifiaciónSegunPaper))]
vias_nad <- CPM_Vias_Clases_Anotacion_DifSign_82p$Via[which(grepl("NAD-consumo", CPM_Vias_Clases_Anotacion_DifSign_82p$ClasifiaciónSegunPaper))]
vias_glicosilacion <- CPM_Vias_Clases_Anotacion_DifSign_82p$Via[which(grepl("Oxidación/glicosilación", CPM_Vias_Clases_Anotacion_DifSign_82p$ClasifiaciónSegunPaper))]

#con vias de clasificación paper:
categoria = "Inflamación"
AR_conacne_inflammation <- AR_conacne[AR_conacne$Pathway %in% vias_inflamation,]

categoria = "Consumo NAD"
AR_conacne_inflammation <- AR_conacne[AR_conacne$Pathway %in% vias_nad,]

categoria = "Glicosilación"
AR_conacne_inflammation <- AR_conacne[AR_conacne$Pathway %in% vias_glicosilacion,]

categoria = "Anti-Aging"
AR_conacne_inflammation <- AR_conacne[AR_conacne$Pathway %in% vias_anti_aging,]

categoria = "Biosíntesis"
AR_conacne_inflammation <- AR_conacne[AR_conacne$Pathway %in% vias_biosintesis,]

categoria = "Oxidación"
AR_conacne_inflammation <- AR_conacne[AR_conacne$Pathway %in% vias_anti_oxidation,]

#con vias de clasificación LEO:
categoria = "Colágeno"
AR_conacne_inflammation <- AR_conacne[AR_conacne$Pathway %in% vias_colageno,]

categoria = "Acidos Grasos"
AR_conacne_inflammation <- AR_conacne[AR_conacne$Pathway %in% vias_acidosgrasos,]

categoria = "Cicatrización"
AR_conacne_inflammation <- AR_conacne[AR_conacne$Pathway %in% vias_cicatrizacion,]

categoria = "Anti-Oxidante"
AR_conacne_inflammation <- AR_conacne[AR_conacne$Pathway %in% vias_antioxidante,]

categoria = "Acido hialurónico"
AR_conacne_inflammation <- AR_conacne[AR_conacne$Pathway %in% vias_acidohialuronico,]

categoria = "Factor hidratante natural"
AR_conacne_inflammation <- AR_conacne[AR_conacne$Pathway %in% vias_hidratantes,]

categoria = "Protección UV"
AR_conacne_inflammation <- AR_conacne[AR_conacne$Pathway %in% vias_uv,]

categoria = "ROS"
AR_conacne_inflammation <- AR_conacne[AR_conacne$Pathway %in% vias_ROS,]

categoria = "Eccema"
AR_conacne_inflammation <- AR_conacne[AR_conacne$Pathway %in% vias_eccema,]

categoria = "Dermatitis Atopica"
AR_conacne_inflammation <- AR_conacne[AR_conacne$Pathway %in% vias_dermatitis,]

categoria = "Infeccion"
AR_conacne_inflammation <- AR_conacne[AR_conacne$Pathway %in% vias_infeccion,]


stats <- Boxplot_PCA_MDS(AR_conacne_inflammation = AR_conacne_inflammation, categoria = categoria)

Boxplot_PCA_MDS <- function(AR_conacne_inflammation, categoria, selec_vias= FALSE) {

  AR_conacne_inflammation <- as.data.frame(t(AR_conacne_inflammation[, -which(colnames(AR_conacne_inflammation) %in% c("Clase", "Description"))]))
  colnames(AR_conacne_inflammation) <- AR_conacne_inflammation[1,]
  AR_conacne_inflammation$ID <- rownames(AR_conacne_inflammation)
  AR_conacne_inflammation <- AR_conacne_inflammation[-1,]

  MetadataB <- as.data.frame(read_excel("~/Daniela/Biota/PipelineBiota/data/Metadata-soloColumnasUsables.xlsx"))
  colnames(MetadataB)
  AR_conacne_met <- merge(AR_conacne_inflammation, MetadataB[, c("ID", "Sexo", "Rango etario", "Edad", "FacilidadBroncearse")], by = "ID")
  str(AR_conacne_met)
  library(reshape2)
  AR_conacne_met_long <- melt(AR_conacne_met, id.vars = c("ID", "Sexo", "Rango etario", "Edad", "FacilidadBroncearse"),
                              variable.name = "Via_Metabolica", value.name = "Abundancia_Relativa")
  str(AR_conacne_met_long)
  AR_conacne_met_long$Abundancia_Relativa <- as.numeric( AR_conacne_met_long$Abundancia_Relativa)
  AR_conacne_met_long$`Rango etario` <- factor(AR_conacne_met_long$`Rango etario`, levels = c("18-35", "35-55", ">55"))
  AR_conacne_met_long$FacilidadBroncearse <- factor(AR_conacne_met_long$FacilidadBroncearse, levels = c("2", "3", "4"))

  library(ggplot2)

  #Grafico general todas las  vias:
  kruskal_results <- AR_conacne_met_long %>%
    summarise(p_value = kruskal.test(Abundancia_Relativa ~ `Rango etario`)$p.value) %>%
    mutate(p_label = paste0("p = ", signif(p_value, 3)))  # Formatear p-valores

  kruskal_results <- AR_conacne_met_long %>%
    summarise(p_value = kruskal.test(Abundancia_Relativa ~ FacilidadBroncearse)$p.value) %>%
    mutate(p_label = paste0("p = ", signif(p_value, 3)))  # Formatear p-valores

  boxp <- ggplot(AR_conacne_met_long, aes(#x = `Rango etario`,
                                          x = FacilidadBroncearse,
                                          y = Abundancia_Relativa)) +
    geom_boxplot() +
    scale_y_log10()+
    labs(x = "Rango etario",
         y = "Abundancia Relativa (log10)",
         y = "CPM (log10)",
         title = sprintf("%s", categoria)) +
    theme_minimal() +
    geom_text(data = kruskal_results, aes(x = 1.5, y = max(AR_conacne_met_long$Abundancia_Relativa, na.rm = TRUE),
                                          label = p_label), inherit.aes = FALSE, size = 4)


  print(boxp)

  # Boxplot por via:
  AR_conacne_met_long <- AR_conacne_met_long %>%
    mutate(Abundancia_Relativa_log = log10(Abundancia_Relativa + 1))  # Evita problemas con ceros

  kruskal_results <- AR_conacne_met_long %>%
    group_by(Via_Metabolica) %>%
    summarise(p_value = kruskal.test(Abundancia_Relativa ~ `Rango etario`)$p.value) %>%
    mutate(p_label = paste0("p = ", signif(p_value, 3)))  # Formatear p-valores

  kruskal_results <- AR_conacne_met_long %>%
    group_by(Via_Metabolica) %>%
    summarise(p_value = kruskal.test(Abundancia_Relativa ~ FacilidadBroncearse)$p.value) %>%
    mutate(p_label = paste0("p = ", signif(p_value, 3)))  # Formatear p-valores

  kruskal_results_log10 <- AR_conacne_met_long %>%
    group_by(Via_Metabolica) %>%
    summarise(p_value = kruskal.test(Abundancia_Relativa_log ~ `Rango etario`)$p.value) %>%
    mutate(p_label = paste0("p = ", signif(p_value, 3)))  # Formatear p-valores

  boxp <- ggplot(AR_conacne_met_long, aes(#x = `Rango etario`,
                                          x = FacilidadBroncearse,
                                          y = Abundancia_Relativa)) +
    geom_boxplot() +
    scale_y_log10()+
    labs(x = "Rango etario",
         y = "Abundancia Relativa (log10)",
         y = "CPM (log10)",
         title = sprintf("%s", categoria)) +
    theme_minimal() +
    #facet_grid(~Via_Metabolica, scales = "free_x") +
    facet_wrap(~Via_Metabolica, scales = "free_x", nrow = 3) + # Cambia "2" por la cantidad de filas que quieras
    geom_text(data = kruskal_results, aes(x = 1.5, y = max(AR_conacne_met_long$Abundancia_Relativa, na.rm = TRUE),
                                          label = p_label), inherit.aes = FALSE, size = 3)


  print(boxp)

  #Si quiero seleccionar algunas vias en particular:
  if( selec_vias == TRUE) {
    vias <- unique(AR_conacne_met_long$Via_Metabolica)
    vias_selec <- vias[6:8]
    vias_selec <- vias[c(2,3,4,5,6,7,10,12)]
    AR_vias_selec <- AR_conacne_met_long[which(AR_conacne_met_long$Via_Metabolica %in% vias_selec),]
    kruskal_results <- AR_vias_selec %>%
      group_by(Via_Metabolica) %>%
      summarise(p_value = kruskal.test(Abundancia_Relativa ~ `Rango etario`)$p.value) %>%
      mutate(p_label = paste0("p = ", signif(p_value, 3)))  # Formatear p-valores

    boxp <- ggplot(AR_vias_selec, aes(x = `Rango etario`, y = Abundancia_Relativa)) +
      geom_boxplot() +
      scale_y_log10()+
      labs(x = "Rango etario",
           y = "Abundancia Relativa (log10)",
           y = "CPM (log10)",
           title = sprintf("%s", categoria)) +
      theme_minimal() +
      #facet_grid(~Via_Metabolica, scales = "free_x") +
      facet_wrap(~Via_Metabolica, scales = "free_x", nrow = 2) +
      geom_text(data = kruskal_results, aes(x = 1.5, y = max(AR_vias_selec$Abundancia_Relativa, na.rm = TRUE),
                                            label = p_label), inherit.aes = FALSE, size = 4)

    print(boxp)

    wilcox_results <- AR_vias_selec[-which(AR_vias_selec$`Rango etario` == "35-55"),] %>%
      group_by(Via_Metabolica) %>%
      summarise(p_value = wilcox.test(Abundancia_Relativa ~ `Rango etario`)$p.value) %>%
      mutate(p_label = paste0("p = ", signif(p_value, 3)))  # Formatear p-valores

    boxp <- ggplot(AR_vias_selec[-which(AR_vias_selec$`Rango etario` == "35-55"),] , aes(x = `Rango etario`, y = Abundancia_Relativa)) +
      geom_boxplot() +
      #scale_y_log10()+
      labs(x = "Rango etario",
           y = "Abundancia Relativa (log10)",
           y = "CPM (log10)",
           title = sprintf("%s", categoria)) +
      theme_minimal() +
      #facet_grid(~Via_Metabolica, scales = "free_x") +
      facet_wrap(~Via_Metabolica, scales = "free_x", nrow = 2) +
      geom_text(data = wilcox_results, aes(x = 1.5, y = max(AR_vias_selec[-which(AR_vias_selec$`Rango etario` == "35-55"),]$Abundancia_Relativa, na.rm = TRUE),
                                           label = p_label), inherit.aes = FALSE, size = 4)
    print(boxp)

  }
  #----------------------------------

  #Wilcox entre grandes y chicos:
  wilcox_results <- AR_conacne_met_long[-which(AR_conacne_met_long$`Rango etario` == "35-55"),] %>%
    group_by(Via_Metabolica) %>%
    summarise(p_value = wilcox.test(Abundancia_Relativa ~ `Rango etario`)$p.value) %>%
    mutate(p_label = paste0("p = ", signif(p_value, 3)))  # Formatear p-valores

  wilcox_results <- AR_conacne_met_long[-which(AR_conacne_met_long$FacilidadBroncearse== "3"),] %>%
    group_by(Via_Metabolica) %>%
    summarise(p_value = wilcox.test(Abundancia_Relativa ~ FacilidadBroncearse)$p.value) %>%
    mutate(p_label = paste0("p = ", signif(p_value, 3)))  # Formatear p-valores

  boxp <- ggplot(AR_conacne_met_long[-which(AR_conacne_met_long$`Rango etario` == "35-55"),] , aes(x = `Rango etario`, y = Abundancia_Relativa)) +
    geom_boxplot() +
    scale_y_log10()+
    labs(x = "Rango etario",
         y = "Abundancia Relativa (log10)",
         y = "CPM (log10)",
         title = sprintf("%s", categoria)) +
    theme_minimal() +
    facet_grid(~Via_Metabolica, scales = "free_x") +
    geom_text(data = wilcox_results, aes(x = 1.5, y = max(AR_conacne_met_long[-which(AR_conacne_met_long$`Rango etario` == "35-55"),]$Abundancia_Relativa, na.rm = TRUE),
                                          label = p_label), inherit.aes = FALSE, size = 4)


  boxp <- ggplot(AR_conacne_met_long[-which(AR_conacne_met_long$FacilidadBroncearse == "3"),] , aes(x = FacilidadBroncearse, y = Abundancia_Relativa)) +
    geom_boxplot() +
    scale_y_log10()+
    labs(x = "Rango etario",
         y = "Abundancia Relativa (log10)",
         y = "CPM (log10)",
         title = sprintf("%s", categoria)) +
    theme_minimal() +
    facet_grid(~Via_Metabolica, scales = "free_x") +
    facet_wrap(~Via_Metabolica, scales = "free_x", nrow = 3) +
    geom_text(data = wilcox_results, aes(x = 1.5, y = max(AR_conacne_met_long[-which(AR_conacne_met_long$FacilidadBroncearse == "3"),]$Abundancia_Relativa, na.rm = TRUE),
                                         label = p_label), inherit.aes = FALSE, size = 3)


  print(boxp)

  library(dplyr)

  summary_stats_log <- AR_conacne_met_long %>%
    group_by(Via_Metabolica, `Rango etario`) %>%
    summarise(
      Media = round(mean(Abundancia_Relativa_log, na.rm = TRUE), 2),
      Mediana = round(median(Abundancia_Relativa_log, na.rm = TRUE), 2),
      Q1 = round(quantile(Abundancia_Relativa_log, 0.25, na.rm = TRUE), 2),
      Q3 = round(quantile(Abundancia_Relativa_log, 0.75, na.rm = TRUE), 2),
      .groups = "drop"
    )

  print(summary_stats_log)

  summary_stats <- AR_conacne_met_long %>%
    group_by(Via_Metabolica, `Rango etario`) %>%
    summarise(
      Media = round(mean(Abundancia_Relativa, na.rm = TRUE), 2),
      Mediana = round(median(Abundancia_Relativa, na.rm = TRUE), 2),
      Q1 = round(quantile(Abundancia_Relativa, 0.25, na.rm = TRUE), 2),
      Q3 = round(quantile(Abundancia_Relativa, 0.75, na.rm = TRUE), 2),
      .groups = "drop"
    )

  print(summary_stats)


  #PCA:
  AR <-  AR_conacne_met[-which(AR_conacne_met$`Rango etario` == "35-55"),]
  #AR <-  AR[, which(colnames(AR) %in% c(as.character(vias_selec), "ID", "Rango etario", "Sexo", "Edad")) ]
  #AR <- AR_conacne_met
  data_pca <- AR %>% select(-ID, -Sexo, -`Rango etario`, -Edad)
  str(data_pca)
  data_pca <- data.frame(lapply(data_pca, as.numeric))
  pca_result <- prcomp(data_pca, center = TRUE, scale. = TRUE)
  pca_df <- as.data.frame(pca_result$x)
  pca_df$`Rango etario` <- AR$`Rango etario`  # Agregar rango etario para colorear
  pca_df$`Rango etario` <- AR$`Rango etario`
  pca_df$ID <- AR$ID
  pca_df$Sexo <- AR$Sexo
  library(ggrepel)
  pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = `Rango etario`, shape = Sexo , label = ID)) +
    geom_point(size = 2, alpha = 0.8) +
    geom_text_repel(size = 3) +  # Evita superposiciones de etiquetas
    labs(title = sprintf("PCA AR - %s", categoria), x = "PC1", y = "PC2", color = "Rango etario") +
    theme_minimal()
  print(pca)

  #MDS:
  dist_matrix <- dist(data_pca)
  mds_result <- cmdscale(dist_matrix, k = 2)
  mds_df <- as.data.frame(mds_result)
  colnames(mds_df) <- c("Dim1", "Dim2")
  mds_df$`Rango etario` <- AR$`Rango etario`
  mds_df$Sexo <- AR$Sexo
  mds_df$ID <- AR$ID
  mds <- ggplot(mds_df, aes(x = Dim1, y = Dim2, color = `Rango etario`, shape = Sexo, label = ID)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_text_repel(size = 3) +  # Evita superposiciones de etiquetas
    labs(title = sprintf("MDS AR - %s", categoria), x = "Dim 1", y = "Dim 2", color = "Rango etario", shape = "Sexo") +
    theme_minimal()
  print(mds)

  return(summary_stats)
}


#-----------------------------------------------------------------------------------------------------------------------
#Me quedo solo con clase:
CPM_vias <- as.data.frame(read_excel("~/Daniela/Biota/Muestras/SubsetPathways/Bo/CPM_Vias_Clases_82p.xlsx"))
str(CPM_vias)
CPM_clases <- CPM_vias %>%
  group_by(Clase) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE), .groups = "drop")

CPM_vias <- CPM_clases
CPM_vias <- AR_conacne
colSums(CPM_vias[,-c(1,2)])
colSums(CPM_vias[,-c(1,2,3)])

CPM_vias <- CPM_vias[,-which(colnames(CPM_vias) == "36")]
#Me quedo con las vias de interes:
CPM_vias <- CPM_vias[which(CPM_vias$Pathway %in% vias_colageno),]

#AR:
AR_vias <- prop.table(as.matrix(CPM_vias[,-c(1,2)]), margin = 2) * 100
AR_vias <- prop.table(as.matrix(CPM_vias[,-c(1,2,3)]), margin = 2) * 100
colSums(AR_vias)
#CPM_vias <- as.data.frame(cbind("Clase" = CPM_vias[,1], AR_vias))
CPM_vias <- as.data.frame(cbind("Via" = CPM_vias[,1], "Clase" = CPM_vias[,2], AR_vias))
CPM_vias <- as.data.frame(cbind("Via" = CPM_vias[,1], "Clase" = CPM_vias[,2], "Descripcion" = CPM_vias[,3], AR_vias))

str(CPM_vias)
colSums(CPM_vias[,-c(1,2)])


#---------------------
library(vegan)

# Selecciona solo las columnas de abundancias relativas de las vías
CPM_vias[,-c(1,2)] <- lapply(CPM_vias[,-c(1,2)], as.numeric)
colSums(CPM_vias[, -c(1,2)])

CPM_vias[,-c(1,2,3)] <- lapply(CPM_vias[,-c(1,2,3)], as.numeric)
colSums(CPM_vias[, -c(1,2,3)])

CPM_viasT <- as.data.frame(t(CPM_vias))
colnames(CPM_viasT) <- CPM_viasT[1,]
CPM_viasT <- CPM_viasT[-c(1,2),]
CPM_viasT <- CPM_viasT[-c(1,2,3),]

#CPM_viasT <- CPM_viasT[-1,]

#Clacular distancia bray curtis:
str(CPM_viasT)
CPM_viasT[] <- sapply(CPM_viasT, as.numeric)
#Elimino outliers:
#CPM_viasT <- CPM_viasT[-which(rownames(CPM_viasT) == "36"),]
CPM_viasT <- CPM_viasT[-which(rownames(CPM_viasT) %in% c("138", "90", "44", "41", "144")),]

library(vegan)
bray_curtis_dist <- vegdist(CPM_viasT, method = "bray")
# Análisis PCoA
pcoa_result <- cmdscale(bray_curtis_dist, eig = TRUE, k = 2)  # k = 2 para dos dimensiones
# Extrae coordenadas principales
pcoa_coords <- as.data.frame(pcoa_result$points)
colnames(pcoa_coords) <- c("PCoA1", "PCoA2")

# Agrega metadata (si la tienes)
pcoa_coords$ID <- rownames(pcoa_coords)  # Identificador de muestra
MetadataB <- as.data.frame(read_excel("~/Daniela/Biota/PipelineBiota/data/Metadata-soloColumnasUsables.xlsx"))

pcoa_coords <- merge(pcoa_coords, MetadataB, by= "ID")
pcoa_coords$Group <- pcoa_coords$Sexo
pcoa_coords$Group <- pcoa_coords$`Rango etario`
pcoa_coords$Group <- pcoa_coords$AfeccionesPiel

ggplot(pcoa_coords, aes(x = PCoA1, y = PCoA2, color = Group)) +
  geom_point(size = 3) +
  #geom_line(aes(group = Group), alpha = 0.5) +  # Conecta puntos del mismo grupo
  geom_text(aes(label = ID), vjust = -1, size = 3) +
  theme_minimal() +
  labs(
    x = paste0("PCoA1: ", round(pcoa_result$eig[1] / sum(pcoa_result$eig) * 100, 1), "% variance"),
    y = paste0("PCoA2: ", round(pcoa_result$eig[2] / sum(pcoa_result$eig) * 100, 1), "% variance"),
    title = "PCoA Bray-Curtis"
  ) +
  theme(legend.position = "right")



#Agrego metadata:
CPM_viasT <- cbind("ID" = rownames(CPM_viasT), CPM_viasT)

MetadataB <- as.data.frame(read_excel("~/Daniela/Biota/PipelineBiota/data/Metadata-soloColumnasUsables.xlsx"))
MetadataB <- MetadataB[which(MetadataB$ID %in% CPM_viasT$ID),]
all(CPM_viasT$ID %in% MetadataB$ID)
all(CPM_viasT$ID == MetadataB$ID)

colnames(MetadataB)
CPM_vias_completo <- merge(CPM_viasT, MetadataB, by = "ID")
str(CPM_vias_completo)
n <- ncol(CPM_viasT)
colnames(CPM_vias_completo)[n]
colnames(CPM_vias_completo)[n+1]
CPM_vias_completo[,2:n] <- lapply(CPM_vias_completo[,2:n], as.numeric)
CPM_vias_completo[,(n+1):ncol(CPM_vias_completo)] <- lapply(CPM_vias_completo[,(n+1):ncol(CPM_vias_completo)], as.factor)



#Calculo pvalues:
p_por_categoria <- data.frame("Categoria" = c(), "Via" = c(), "Wilcoxon"= c(), "GrupoDominante" = c())

categorias <- colnames(CPM_vias_completo)[(n+1):ncol(CPM_vias_completo)]
vias <- colnames(CPM_vias_completo)[2:n]
c= n+1
j = 2
i=1

#c= n+4
for (c in (n+2):ncol(CPM_vias_completo)) {
  categoria <- colnames(CPM_vias_completo)[c]
  print(categoria)

  for (j in 2:(n)) {

    via <- colnames(CPM_vias_completo)[j]
    print(via)
    df_sin_na <- CPM_vias_completo[complete.cases(CPM_vias_completo[, categoria]), ]

    if (any(df_sin_na[,categoria] == "Desconocido")) {
      df_sin_na <- df_sin_na[-which(df_sin_na[,categoria] == "Desconocido"),]
    }

    if (any(is.na(df_sin_na[,categoria]))) {
      df_sin_na <- df_sin_na[-which(is.na(df_sin_na[,categoria])),]
    }

    if( categoria == "Maquillaje_Base") {
      df_sin_na <- df_sin_na[-which(df_sin_na$Sexo == "Masculino"),]
    }

    if( categoria == "AfeccionesPiel") {
      df_sin_na <- df_sin_na[-which(df_sin_na$`Rango etario` != "18-35"),]
    }

    # Calcular la media para cada grupo en la categoría para sacar el grupo dominante:
    mean_values <- aggregate(df_sin_na[, via], by = list(df_sin_na[, categoria]), FUN = mean)
    colnames(mean_values) <- c("Grupo", "Media")
    grupo_dominante <- mean_values$Grupo[which.max(mean_values$Media)]
    print(grupo_dominante)

    if( length(levels(df_sin_na[, categoria])) > 2 ) {
      k_test <- kruskal.test(df_sin_na[, via] ~ df_sin_na[, categoria], data = df_sin_na)
      p_v <- k_test$p.value
      p_valor <- p_v
      p_val <- p_v
      p_adj <- p_v
    } else if (length(levels(df_sin_na[, categoria])) == 2 ) {
      w_test <- wilcox.test(df_sin_na[,via] ~ df_sin_na[,categoria], df_sin_na)
      p_valor <- w_test$p.value
    } else {
      p_valor <- 1
    }

    p_por_categoria[i,"Categoria"] <- categoria
    p_por_categoria[i,"Via"] <- via
    p_por_categoria[i,"Wilcoxon"] <- p_valor
    p_por_categoria[i, "GrupoDominante"] <-  as.character(grupo_dominante)

    i <- i+1
  }
}

write.xlsx(p_por_categoria, file =  "~/Daniela/Biota/Muestras/SubsetPathways/Bo/Bo_p_por_categoria-274vias_CPM_82p.xlsx")


save(p_por_categoria, file = "~/Daniela/Biota/Muestras/SubsetPathways/Bo/p_por_categoria_65clases_82p.RData")
load("~/Daniela/Biota/Muestras/SubsetPathways/Bo/p_por_categoria_65clases_82p.RData")

save(p_por_categoria, file = "~/Daniela/Biota/Muestras/SubsetPathways/Bo/p_por_categoria_65clasesAR_82p.RData")

save(p_por_categoria, file = "~/Daniela/Biota/Muestras/SubsetPathways/Bo/p_por_categoria-274vias_82p.RData")
load("~/Daniela/Biota/Muestras/SubsetPathways/Bo/p_por_categoria-274vias_82p.RData")

save(p_por_categoria, file = "~/Daniela/Biota/Muestras/SubsetPathways/Bo/p_por_categoria-274vias_CPM_82p.RData")
load("~/Daniela/Biota/Muestras/SubsetPathways/Bo/p_por_categoria-274vias_CPM_82p.RData")

p_significativos <- p_por_categoria[which(p_por_categoria$Wilcoxon<0.05),]
p_significativos <- p_por_categoria
vias_sign <- unique(p_significativos$Via)

rowSums(CPM_vias_completo[,2:n])

write.xlsx(CPM_vias_completo, file = "~/Daniela/Biota/Muestras/SubsetPathways/Bo/CPM_viasAR_completo_sinoutliers.xlsx")
write.xlsx(CPM_vias_completo, file = "~/Daniela/Biota/Muestras/SubsetPathways/Bo/CPM_viasAR_completo.xlsx")
write.xlsx(CPM_vias_completo, file = "~/Daniela/Biota/Muestras/SubsetPathways/Bo/CPM_viasCPM_completo.xlsx")
write.xlsx(CPM_vias_completo, file = "~/Daniela/Biota/Muestras/SubsetPathways/Bo/CPM_clases_completo.xlsx")
write.xlsx(CPM_vias_completo, file = "~/Daniela/Biota/Muestras/SubsetPathways/Bo/CPM_clasesAR_completo.xlsx")

CPM_vias_completo <- as.data.frame(read_excel("~/Daniela/Biota/Muestras/SubsetPathways/Bo/CPM_viasAR_completo.xlsx"))

#Boxplots de las vias con dif significativas: -----------------------------------

unique(p_significativos$Categoria)
categoria = "Sexo"
categoria = "Rango etario"
categoria = "Alcohol"
categoria = "Tabaco"
categoria = "Tipodepiel"
categoria = "AplicacionProtectorSolar"
categoria = "CuándoLimpiaCara"
categoria = "CuándoMaquillaje"
categoria = "AfeccionesPiel"

print(categoria)
p <- p_significativos[which(p_significativos$Categoria == categoria),]
if(categoria == "Rango etario") {
  p$Categoria <- "RangoEtario"
  colnames(CPM_vias_completo)[which(colnames(CPM_vias_completo) == "Rango etario")] <- "RangoEtario"
}
p$Via
p$Categoria

i=1
lista_graficos <- list()
for (i in 1:nrow(p)) {
  categoria <- p[i, "Categoria"]
  via <- p[i, 2]

  df_sin_na <- CPM_vias_completo[complete.cases(CPM_vias_completo[,categoria]), ]
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

  if(categoria == "RangoEtario") {
    df_sin_na[,categoria] <- factor(df_sin_na[,categoria], levels = c("18-35", "35-55", ">55"))
  }

  if(categoria == "AfeccionesPiel") {
    df_sin_na[,categoria] <- factor(df_sin_na[,categoria], levels = c("Ninguna", "Acné"))
    df_sin_na <- df_sin_na[-which(df_sin_na$`Rango etario` != "18-35"),]
  }

  g <- ggplot(df_sin_na, aes_string(x = df_sin_na[,categoria], y = df_sin_na[,via])) +
    geom_boxplot() +
    theme(axis.text.x = element_text(size = 7),
          axis.text.y = element_text(size = 8)) +
    theme(axis.title.x = element_text(size = 7),
          axis.title.y = element_text(size = 8)) +
    labs(y = paste(via), x = paste(categoria))

  lista_graficos[[i]] <- g
  i <- i+1
}

length(lista_graficos)
grid.arrange(grobs = lista_graficos, ncol = 3)

grid.arrange(grobs = lista_graficos[1:20], ncol = 5)
grid.arrange(grobs = lista_graficos[21:41], ncol = 5)
grid.arrange(grobs = lista_graficos[42:69], ncol = 5)
grid.arrange(grobs = lista_graficos[63:89], ncol = 5)


grid.arrange(grobs = lista_graficos[c(1,2)], ncol = 2)
p$Via[c(1,2)]
grid.arrange(grobs = lista_graficos[c(23)], ncol = 4)
p$Via[c(4,7,8,13,15,17,19,20)]


#Sexo vias ---------------------------------------------------------------
  #aumento masculino
listado_masc <- c(1,3,8,12, 15,26, 29 , 31, 32)
p$Via[listado_masc]
grid.arrange(grobs = lista_graficos[listado_masc], ncol = 4)
  #aumento femenino
listado_fem <- c(2, 16, 20, 21)
grid.arrange(grobs = lista_graficos[listado_fem], ncol = 4)
p$Via[listado_fem]

  #Anotacion
CPM_vias$DifSign_Via_AR_Sexo <- NA
#CPM_vias <- CPM_vias[, c(1,2,ncol(CPM_vias), 3:84)]

CPM_vias$DifSign_Via_AR_Sexo[which(CPM_vias$Via %in% p$Via[listado_masc])] <- "Aumento en masculino"
CPM_vias$DifSign_Via_AR_Sexo[which(CPM_vias$Via %in% p$Via[listado_fem])] <- "Disminución en masculino"

  # PARA CPM:
listado_masc <- c(1,2,6,7,11,12,14,19, 25, 29, 32)
p$Via[listado_masc]
grid.arrange(grobs = lista_graficos[listado_masc], ncol = 6)

listado_fem <- c(15,20, 26, 28)
grid.arrange(grobs = lista_graficos[listado_fem], ncol = 4)
p$Via[listado_fem]

CPM_vias$DifSign_Via_CPM_Sexo <- NA
CPM_vias$DifSign_Via_CPM_Sexo[which(CPM_vias$Via %in% p$Via[listado_masc])] <- "Aumento en masculino"
CPM_vias$DifSign_Via_CPM_Sexo[which(CPM_vias$Via %in% p$Via[listado_fem])] <- "Disminución en masculino"


#Rango Etario vias ---------------------------------------------------------------
  #aumento >55
mas50 <- c(4,11,15,16,18, 23, 32, 51, 52, 54, 56,67, 80, 85, 87)
p$Via[mas50]
grid.arrange(grobs = lista_graficos[mas50], ncol = 5)
  #aumento <55
menos50 <- c(2,3,7,13,14,17,29, 41, 43, 49,62, 64, 66, 69, 86, 90)
grid.arrange(grobs = lista_graficos[menos50], ncol = 4)
p$Via[menos50]

CPM_vias$DifSign_Via_AR_Rango <- NA
CPM_vias$DifSign_Via_AR_Rango[which(CPM_vias$Via %in% p$Via[mas50])] <- "Aumento en >55 años"
CPM_vias$DifSign_Via_AR_Rango[which(CPM_vias$Via %in% p$Via[menos50])] <- "Disminución en >55 años"

  #PARA CPM:
mas50 <- c(7,11,12,14, 17, 18, 25, 27, 36, 37, 39, 40, 45, 54, 59, 61)
p$Via[mas50]
grid.arrange(grobs = lista_graficos[mas50], ncol = 4)

menos50 <- c(2,6,9,10,13,23,29, 30, 31, 32, 34, 35, 42, 44,47, 51,60)
grid.arrange(grobs = lista_graficos[menos50], ncol = 6)
p$Via[menos50]

CPM_vias$DifSign_Via_CPM_Rango <- NA
CPM_vias$DifSign_Via_CPM_Rango[which(CPM_vias$Via %in% p$Via[mas50])] <- "Aumento en >55 años"
CPM_vias$DifSign_Via_CPM_Rango[which(CPM_vias$Via %in% p$Via[menos50])] <- "Disminución en >55 años"


#Tipo de piel vias ---------------------------------------------------------------
  #aumento Piel grasa
masGrasa <- c(2, 17, 20, 27, 30, 37, 38, 40, 41, 42, 43,45, 46, 59)
p$Via[masGrasa]
grid.arrange(grobs = lista_graficos[masGrasa], ncol = 5)

#Agrego anotacion:
CPM_vias$DifSign_Via_AR_Piel <- NA
CPM_vias$DifSign_Via_AR_Piel[which(CPM_vias$Via %in% p$Via[masGrasa])] <- "Aumento en piel grasa"

  #PARA CPM:
masGrasa <- c(2, 16, 17, 20, 22, 24, 25, 26, 30, 32,33, 34, 35, 36, 37,44,47,57, 59, 60, 62,64,  66, 67, 68, 69 )
p$Via[masGrasa]
grid.arrange(grobs = lista_graficos[masGrasa], ncol = 7)

masSeca <- c(38)
p$Via[masSeca]
grid.arrange(grobs = lista_graficos[masSeca], ncol = 5)

#Agrego anotacion:
CPM_vias$DifSign_Via_CPM_Piel <- NA
CPM_vias$DifSign_Via_CPM_Piel[which(CPM_vias$Via %in% p$Via[masGrasa])] <- "Aumento en piel grasa"
CPM_vias$DifSign_Via_CPM_Piel[which(CPM_vias$Via %in% p$Via[masSeca])] <- "Aumento en piel seca"



# UNIENDO TODAS LAS ANOTACIONES:
CPM_vias <- CPM_vias %>%
  mutate(
    DifSign_Via_AR = case_when(
      !is.na(DifSign_Via_AR_Piel) & !is.na(DifSign_Via_AR_Sexo) & !is.na(DifSign_Via_AR_Rango) ~ paste(DifSign_Via_AR_Piel, DifSign_Via_AR_Sexo, DifSign_Via_AR_Rango, sep = ","),
      !is.na(DifSign_Via_AR_Piel) & !is.na(DifSign_Via_AR_Sexo) ~ paste(DifSign_Via_AR_Piel, DifSign_Via_AR_Sexo, sep = ","),
      !is.na(DifSign_Via_AR_Piel) & !is.na(DifSign_Via_AR_Rango) ~ paste(DifSign_Via_AR_Piel, DifSign_Via_AR_Rango, sep = ","),
      !is.na(DifSign_Via_AR_Sexo) & !is.na(DifSign_Via_AR_Rango) ~ paste(DifSign_Via_AR_Sexo, DifSign_Via_AR_Rango, sep = ","),
      !is.na(DifSign_Via_AR_Piel) ~ DifSign_Via_AR_Piel,
      !is.na(DifSign_Via_AR_Sexo) ~ DifSign_Via_AR_Sexo,
      !is.na(DifSign_Via_AR_Rango) ~ DifSign_Via_AR_Rango,
      TRUE ~ NA_character_
    )
  )

CPM_vias <- CPM_vias %>%
  mutate(
    DifSign_Via_CPM = case_when(
      !is.na(DifSign_Via_CPM_Piel) & !is.na(DifSign_Via_CPM_Sexo) & !is.na(DifSign_Via_CPM_Rango) ~ paste(DifSign_Via_CPM_Piel, DifSign_Via_CPM_Sexo, DifSign_Via_CPM_Rango, sep = ","),
      !is.na(DifSign_Via_CPM_Piel) & !is.na(DifSign_Via_CPM_Sexo) ~ paste(DifSign_Via_CPM_Piel, DifSign_Via_CPM_Sexo, sep = ","),
      !is.na(DifSign_Via_CPM_Piel) & !is.na(DifSign_Via_CPM_Rango) ~ paste(DifSign_Via_CPM_Piel, DifSign_Via_CPM_Rango, sep = ","),
      !is.na(DifSign_Via_CPM_Sexo) & !is.na(DifSign_Via_CPM_Rango) ~ paste(DifSign_Via_CPM_Sexo, DifSign_Via_CPM_Rango, sep = ","),
      !is.na(DifSign_Via_CPM_Piel) ~ DifSign_Via_CPM_Piel,
      !is.na(DifSign_Via_CPM_Sexo) ~ DifSign_Via_CPM_Sexo,
      !is.na(DifSign_Via_CPM_Rango) ~ DifSign_Via_CPM_Rango,
      TRUE ~ NA_character_
    )
  )


CPM_Vias_Clases_Anotacion_DifSign_82p <- read_excel("~/Daniela/Biota/CPM_Vias_Clases_Anotacion_DifSign_82p.xlsx")

all(CPM_vias$Via == CPM_Vias_Clases_Anotacion_DifSign_82p$Via)
CPM_vias <- CPM_vias[order(match(CPM_Vias_Clases_Anotacion_DifSign_82p$Via, CPM_vias$Via)),]

CPM_vias <- CPM_vias[,c(1,2, 88, 3:87)]
write.xlsx(CPM_vias, file = "~/Daniela/Biota/Muestras/SubsetPathways/AR_vias_clase_DifSignAR.xlsx")
write.xlsx(CPM_vias, file = "~/Daniela/Biota/Muestras/SubsetPathways/CPM_vias_clase_DifSignCPM.xlsx")

########################################################################
#Heatmap presencia/ausencia: --------------------------------------------

CPM_vias_completo <- as.data.frame(read_excel("~/Daniela/Biota/Muestras/SubsetPathways/Bo/CPM_clases_completo.xlsx"))

CPM_vias_completo <- as.data.frame(read_excel("~/Daniela/Biota/Muestras/SubsetPathways/Bo/CPM_vias_completo.xlsx"))
  #Eliminar vias que etsan presentes en todas las muestras: # Eliminar las columnas donde todos los valores son > 0
colnames(CPM_vias_completo)[275]
colnames(CPM_vias_completo)[67] #fechadenacimiento

CPM_onoff_completo <- CPM_vias_completo[, -which(colSums(CPM_vias_completo[, 2:275] > 0) == nrow(CPM_vias_completo))]
CPM_onoff_completo <- CPM_vias_completo[, -which(colSums(CPM_vias_completo[, 2:66] > 0) == nrow(CPM_vias_completo))]

# Eliminar las columnas donde todos los valores son > 0
cols_to_remove <- apply(CPM_vias_completo[, 2:275], 2, function(col) all(col > 0))
cols_to_remove <- apply(CPM_vias_completo[, 2:66], 2, function(col) all(col > 0))

cols_to_remove <- names(cols_to_remove[which(cols_to_remove == TRUE)])
cols_to_remove
CPM_onoff_completo <- CPM_vias_completo[, -which(colnames(CPM_vias_completo) %in% cols_to_remove)]

ncol(MetadataB)
colnames(CPM_onoff_completo)
colnames(CPM_onoff_completo)[ncol(CPM_onoff_completo)-31]


summary(CPM_onoff_completo[,2:274])
summary(CPM_onoff_completo[,2:63])

colnames(CPM_vias_completo)[275]

  #CPM:
CPM_onoff_completo_binary <- CPM_onoff_completo[, 2:63]
  #AR:
CPM_onoff_completo_binary <- prop.table(as.matrix(CPM_onoff_completo[,2:274]), margin = 1) * 100
CPM_onoff_completo_binary <- prop.table(as.matrix(CPM_onoff_completo[,2:63]), margin = 1) * 100

rowSums(CPM_onoff_completo_binary)
  #on/off:
CPM_onoff_completo_binary <- (CPM_onoff_completo[, 2:274])
CPM_onoff_completo_binary <- (CPM_onoff_completo[, 2:63])
CPM_onoff_completo_binary <-ifelse(CPM_onoff_completo_binary >0, 1,0)
  #me quedo con las vias significativas:
CPM_onoff_completo_binary <- CPM_onoff_completo[, which(colnames(CPM_onoff_completo) %in% vias_sign)]
#CPM_onoff_completo_binary <-ifelse(CPM_onoff_completo_binary >0, 1,0)


CPM_onoff_completo_binary <- t(CPM_onoff_completo_binary)
colnames(CPM_onoff_completo_binary) <- CPM_onoff_completo$ID

colnames(CPM_vias_completo)[276:ncol(CPM_vias_completo)]

jaccard_dist = function(x, y) 1 - sum(x & y)/sum(x | y)

library(ComplexHeatmap)
Ms <- t(scale(t(CPM_onoff_completo_binary)))
hm <- Heatmap(
  #Ms,
  CPM_onoff_completo_binary,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  clustering_distance_columns = jaccard_dist,
  clustering_distance_rows = jaccard_dist,
  column_split = CPM_vias_completo$`Rango etario`,
  top_annotation = HeatmapAnnotation(Rango = CPM_vias_completo$`Rango etario`,
                                     TipoPiel = CPM_vias_completo$Tipodepiel,
                                     ProtectorSolar = CPM_vias_completo$AplicacionProtectorSolar,
                                     CuandoMaquilla = CPM_vias_completo$CuándoMaquillaje,
                                     #ActividadFisica = CPM_vias_completo$ActividadFísica,
                                     Tabaco = CPM_vias_completo$Tabaco,
                                     Alcohol = CPM_vias_completo$Alcohol
                                     #TratamientoMedico = CPM_vias_completo$TratamientoMédico,
                                     #AntecedentesEnf = CPM_vias_completo$AntecedenteEnfermedad
                                     ),  # Asume que 'grupo' es tu columna de metadata
  col = colorRampPalette(c("red", "green"))(50),
  column_names_gp = gpar(fontsize = 7),  # Reducir tamaño de los nombres de las columnas
  row_names_gp = gpar(fontsize = 8),     # Reducir tamaño de los nombres de las filas
  #width = unit(10, "cm"),                # Ajustar el tamaño
  #height = unit(8, "cm")                 # Ajustar el tamaño
)

print(hm)

column_dist <- dist(t(Ms))  # Usar 't()' porque queremos trabajar con columnas
column_dist <- dist(Ms)
dend_H_columns <- hclust(column_dist)
plot(dend_H_columns, cex = 0.6, main = "Dendrograma de las columnas")
plot(dend_H_columns, cex = 0.6, main = "Dendrograma de las vias")
cut_dend <- cutree(dend_H_columns, k=2)
table(cutree(dend_H_columns, k=2))

# Agrupar vias en clases: -------------------------------------------------
#https://metacyc.org/group?id=biocyc14-14708-3818509268
clases_vias <- read.table("~/Daniela/Biota/Muestras/SubsetPathways/ontology_All_pathways_MetaCyc.txt", header = TRUE, sep = "\t")
CPM_vias_completo <- as.data.frame(read_excel("~/Daniela/Biota/Muestras/SubsetPathways/Bo/CPM_vias_completo.xlsx"))

colnames(CPM_vias_completo)[276]
all(colnames(CPM_vias_completo)[2:275] %in% clases_vias$Pathways)
solovias <- as.data.frame(t(CPM_vias_completo[,2:275]))
colnames(solovias) <- CPM_vias_completo$ID
solovias$Via <- rownames(solovias)

v=1
for (v in 1:nrow(solovias)) {
  via <- solovias$Via[v]
  print(via)
  clase <- clases_vias$Ontology...parents.of.class[which(clases_vias$Pathways == via)]
  solovias$Clase[v] <- clase
}

length(unique(solovias$Clase))
solovias <- cbind("Via" = solovias$Via, "Clase" = solovias$Clase, solovias[,-which(colnames(solovias) %in% c("Clase", "Via"))])
write.xlsx(solovias, file = "~/Daniela/Biota/Muestras/SubsetPathways/Bo/CPM_Vias_Clases_82p.xlsx")


#PLOT WITH HUMAN: ------------------------------------------------------------
de_host_file <- "Bo"
categoria <- "Rango"

pathway <- "PWY-6147"
pathway <- "HEME-BIOSYNTHESIS-II"
pathway <- "PWY-6125" #colagenos
pathway <- "1CMET2-PWY"
pathway <- "ARO-PWY"
pathway <- "PWY-5918"
pathway <- "PWY66-389"

cpm_table <- path.expand(sprintf("~/Daniela/Biota/Muestras/SubsetPathways/%s/subset_CPM_pathabundance.tsv", de_host_file))
path_abundance <- as.data.frame(read_tsv(cpm_table))
colnames(path_abundance)[1] <-"FEATURE"
colnames(path_abundance)[-1] <- gsub(sprintf("%s_concatR1R2_Abundance-CPM", de_host_file), "", colnames(path_abundance)[-1] )

MetadataB <- as.data.frame(read_excel("~/Daniela/Biota/PipelineBiota/data/Metadata-soloColumnasUsables.xlsx"))
all(colnames(path_abundance)[-1] %in% MetadataB$ID)
MetadataB <- MetadataB[which(MetadataB$ID %in% colnames(path_abundance)),]

# Reordenar las columnas (excepto la primera)
path_ab <- path_abundance[, c(1, match(MetadataB$ID, colnames(path_abundance)[-1]) + 1)]
all(colnames(path_ab)[-1] == MetadataB$ID)

categoria <- "Rango"
categoria <- "Alcohol"
colnames(MetadataB)
metadato_fila <- c(categoria, MetadataB[,categoria])
path_ab <- rbind(metadato_fila, path_ab)

rango_cpm_tsv <- sprintf("~/Daniela/Biota/Muestras/SubsetPathways/%s/%s_CPM_pathways.tsv", de_host_file, categoria)
write_tsv(path_ab, file = rango_cpm_tsv)

barplot_file <- sprintf("~/Daniela/Biota/Muestras/SubsetPathways/%s/CPM_%s-%s_barplot.png", de_host_file, categoria, pathway)
command <- paste(
  "bash -c 'source activate biobakery3 &&",  # Activar el entorno
  "humann_barplot",                                # Ejecutar humann
  sprintf("--input %s", rango_cpm_tsv),
  sprintf("--focal-metadata %s --last-metadata %s", categoria, categoria),
  sprintf("--output %s", barplot_file),
  sprintf("--focal-feature %s", pathway),
  "--sort sum metadata --scaling logstack'",
  sep = " "
)
system(command = command, intern = FALSE)



#Grafico individual:
rango_cpm_tsv <- "~/Daniela/Biota/Muestras/SubsetPathways/108Bo_CPM_pathabundance.tsv"
command <- paste(
  "bash -c 'source activate biobakery3 &&",  # Activar el entorno
  "humann_barplot",                                # Ejecutar humann
  sprintf("--input %s", rango_cpm_tsv),
  #sprintf("--focal-metadata %s --last-metadata %s", categoria, categoria),
  sprintf("--output %s", barplot_file),
  sprintf("--focal-feature %s'", pathway),
  #"--sort sum metadata --scaling logstack'",
  sep = " "
)
system(command = command, intern = FALSE)

#################################################################################





