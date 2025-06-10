#Agregar las muestras de acné:
otus <- generateOTUsTableGrupal(patients_dir = "~/Daniela/Biota/Muestras/73m", source = "KRAKEN", conEukaryota = FALSE, de_host = "Biota")

otus <- generateOTUsTableGrupal(patients_dir = "~/Daniela/Biota/Muestras/73m", source = "KRAKEN", conEukaryota = FALSE, de_host = "Biota")
AR_SubSpecies_KRAKEN <- group_TaxonomicLevels(patients_dir = "~/Daniela/Biota/Muestras/73m", tabla_otus = otus, source = "KRAKEN", conEukaryota = FALSE, nombre_extra = "conAcne")
#TENGO QUE VOLVER A GENERAR EL GROUP TAXONOMIC LEVELS DE LOS 83 PACIENTES ORIGINALES!!!!

AR_SubSpecies_KRAKEN <- AR_SubSpecies_KRAKEN[[1]]
AR_SubSpecies_KRAKEN <- AR_SubSpecies_KRAKEN[[9]]
colSums(AR_SubSpecies_KRAKEN[,-1])

#me quedo con las 138 especies:
load("~/Daniela/Biota/Muestras/73m/SubEspecies_AUSAR.RData")
AR_SubSpecies_KRAKEN <- AR_SubSpecies_KRAKEN[which(AR_SubSpecies_KRAKEN$SubSpecies %in% subespecies_ausar),]
AR_SubSpecies_K <-  as.data.frame(prop.table(as.matrix(AR_SubSpecies_KRAKEN[,-1]), margin = 2) * 100)
AR_SubSpecies_KRAKEN <- cbind(AR_SubSpecies_KRAKEN$SubSpecies, AR_SubSpecies_K)
colSums(AR_SubSpecies_K)
colSums(AR_SubSpecies_KRAKEN[,-1])


#AGREGARLOS A LA METADATAB
patient_dir <- "~/Daniela/Biota/Muestras/73m/25"
patient_dir <- "~/Daniela/Biota/Muestras/73m/217"

id <- basename(patient_dir)
patients_dir <- dirname(patient_dir)
MetadataB <- as.data.frame(read_excel("~/Daniela/Biota/PipelineBiota/data/Metadata-soloColumnasUsables.xlsx"))
BIOTALIFE_SKIN_Respuestas_ <- as.data.frame(read_excel("~/Daniela/Biota/BIOTALIFE SKIN  (Respuestas).xlsx", sheet = 4))
#BIOTALIFE_SKIN_Respuestas_ <- read_excel("~/Daniela/Biota/BIOTALIFE SKIN- ActualizadoOct.xlsx")

#Obtener info gral
BIOTALIFE_SKIN_Respuestas_$ID  <- gsub("\\.0$", "", as.character(BIOTALIFE_SKIN_Respuestas_$ID))
BIOTALIFE_SKIN_Respuestas_$ID[which(BIOTALIFE_SKIN_Respuestas_$ID == "118")] <- "118-1"
BIOTALIFE_SKIN_Respuestas_$ID[which(BIOTALIFE_SKIN_Respuestas_$ID == "184")] <- "184-1"


#Agrego info del paciente si es que NO esta en MetadabaB pero SI en BiotaLife respuestas: --------------------------
#MetadataB <- MetadataB[, -which(colnames(MetadataB) == "TipoPielViejo")]
#MetadataB <- MetadataB[, -which(colnames(MetadataB) == "Secuenciado")]
colnames(BIOTALIFE_SKIN_Respuestas_)[c(2,5,6,7,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,25,26,27,28,29,30,31,33,34,57, 60, 61, 32,35)] <- c("ID", "FechadeNacimiento", "Sexo", "Email",
                                                                                                                                               "ColorCabello", "ColorOjos", "FacilidadBroncearse", "Peso",
                                                                                                                                               "Altura", "AntecedenteEnfermedad", "CualEnfermedad", "AfeccionesPiel",
                                                                                                                                               "OtraAfeccionPiel", "TratamientoEstético3meses", "MétodoAnticonceptivo", "EmbarazadaoAmantando",
                                                                                                                                               "TratamientoMédico", "CualTratamientoMedico", "FrecuenciaTratamientoMedico", "Tabaco",
                                                                                                                                               "CigarrillosporSemana", "Alcohol", "AlcoholporSemana", "ActividadFísica",
                                                                                                                                               "CuándoLimpiaCara", "ConquéLimpiaCara", "AplicacionProtectorSolar", "CuándoMaquillaje",
                                                                                                                                               "Fecha Cita", "Edad", "Rango etario",
                                                                                                                                               "Tipodepiel", "Maquillaje_Base")

BIOTALIFE_SKIN_Respuestas_[which(BIOTALIFE_SKIN_Respuestas_$ID == id),]
if(!(id %in% MetadataB$ID) & id %in% BIOTALIFE_SKIN_Respuestas_$ID) {
  MetadataB <- rbind(MetadataB, BIOTALIFE_SKIN_Respuestas_[which(BIOTALIFE_SKIN_Respuestas_$ID == id), which(colnames(BIOTALIFE_SKIN_Respuestas_) %in% colnames(MetadataB))])
  #modifico lo de maquillaje base y lo de facilidad de bronceado
  MetadataB$FacilidadBroncearse <- sub("^(\\d+).*", "\\1", MetadataB$FacilidadBroncearse)
  MetadataB$FacilidadBroncearse[which(MetadataB$FacilidadBroncearse == "1")] <- "2"
  MetadataB$CuándoMaquillaje[which(MetadataB$CuándoMaquillaje == "1 o 2 veces por semana")] <- "Diariamente"
  MetadataB$AplicacionProtectorSolar[which(MetadataB$AplicacionProtectorSolar == "1 o 2 veces por semana")] <- "Solo ante exposición en verano"

  if(any(MetadataB$Maquillaje_Base != "0" | MetadataB$Maquillaje_Base != "1")) {
    maq_base <-  MetadataB$Maquillaje_Base[-which(MetadataB$Maquillaje_Base == "0" | MetadataB$Maquillaje_Base == "1")]
    MetadataB$Maquillaje_Base[which(MetadataB$Maquillaje_Base == maq_base)] <- ifelse(grepl("base", maq_base) | grepl("corrector", maq_base) | grepl("polvo", maq_base), "1", "0")
  }
  MetadataB$Maquillaje_Base[which(is.na(MetadataB$Maquillaje_Base))] <- "0"
}

write.xlsx(MetadataB, file = "~/Daniela/Biota/PipelineBiota/data/Metadata-soloColumnasUsables.xlsx")

#------------------------------------------------------------------------------------ -------------------

AR_SubSpecies_KRAKEN <- as.data.frame(read_excel("~/Daniela/Biota/Muestras/73m/AR__SubSpecies_KRAKEN.xlsx"))
AR_SubSpecies_KRAKEN <- as.data.frame(read_excel("~/Daniela/Biota/Muestras/73m/AR__Species_KRAKEN.xlsx"))
AR_SubSpecies_KRAKEN <- read_excel("~/Daniela/Biota/Muestras/73m/AR__Genus_KRAKEN.xlsx")


str(AR_SubSpecies_KRAKEN)
rowMeans(AR_SubSpecies_KRAKEN[which(AR_SubSpecies_KRAKEN$Species == "Betapapillomavirus 5"),-1])

#Quedarme con algunas subespecies
Biomarcadores_PielGrasa <- read_excel("~/Daniela/Biota/PipelineBiota/data/Biomarcadores_TipodePiel_sinoutliers.xlsx")
subespecies_piel <- unique(Biomarcadores_PielGrasa$Subespecie)

load("~/Daniela/Biota/p_por_categoria-138subspecies_82p.RData")
subespecies_piel <- unique(p_por_categoria$SubEspecie)

AR_SubSpecies_KRAKEN <- AR_SubSpecies_KRAKEN[which(AR_SubSpecies_KRAKEN$SubSpecies %in% subespecies_piel),]
AR_SubSpecies_K <-  as.data.frame(prop.table(as.matrix(AR_SubSpecies_KRAKEN[,-1]), margin = 2) * 100)
AR_SubSpecies_KRAKEN <- cbind(AR_SubSpecies_KRAKEN$SubSpecies, AR_SubSpecies_K)
colSums(AR_SubSpecies_K)
colSums(AR_SubSpecies_KRAKEN[,-1])

#-----------------------------

AR_subespeciesT <- as.data.frame(t(AR_SubSpecies_KRAKEN))
colnames(AR_subespeciesT) <- AR_subespeciesT[1,]
AR_subespeciesT <- AR_subespeciesT[-1,]
AR_subespeciesT <- cbind("ID" = rownames(AR_subespeciesT), AR_subespeciesT)
colnames(AR_subespeciesT)[1] <- "ID"
AR_subespeciesT$ID <- sub("_KRAKEN", "", AR_subespeciesT$ID)

MetadataB <- read_excel("~/Daniela/Biota/PipelineBiota/data/Metadata-soloColumnasUsables.xlsx")

AR_subespeciesT$ID[which(AR_subespeciesT$ID == "118-1")] <- "118"
#MetadataB$ID[which(MetadataB$ID == "118-1")] <- "118"
#write.xlsx(MetadataB, file ="~/Daniela/Biota/PipelineBiota/data/Metadata-soloColumnasUsables.xlsx")

setdiff(MetadataB$ID, AR_subespeciesT$ID)

df_completo <- merge(AR_subespeciesT, MetadataB, by = "ID")
df_completo$FacilidadBroncearse[which(df_completo$FacilidadBroncearse == "1")] <- "2"
df_completo$CuándoMaquillaje[which(df_completo$CuándoMaquillaje == "1 o 2 veces por semana")] <- "Diariamente"
df_completo$AplicacionProtectorSolar[which(df_completo$AplicacionProtectorSolar == "1 o 2 veces por semana")] <- "Solo ante exposición en verano"

n <- ncol(df_completo) - ncol(MetadataB)
colnames(df_completo)[n+2] #debe ser fehca de nacimiento
str(df_completo)

#Hacer factor las categorias y numericas las AR:
df_completo[,(n+2):ncol(df_completo)] <- lapply(df_completo[,(n+2):ncol(df_completo)], as.factor)
df_completo[,2: (n+1)] <- lapply(df_completo[,2: (n+1)], as.numeric)
rowSums(df_completo[,2:n+1])

#write.xlsx(df_completo, file = "~/Daniela/Biota/df_completo_85_138SubSpecies.xlsx")
write.xlsx(df_completo, file = "~/Daniela/Biota/df_completo_83_138SubSpecies.xlsx")
write.xlsx(df_completo, file = "~/Daniela/Biota/df_completo_83_11780Species.xlsx")
write.xlsx(df_completo, file = "~/Daniela/Biota/df_completo_82_142Genus.xlsx")


niveles <- c("SubSpecies", "Species", "Genus", "Phylum")
nivel <- niveles[1]


colnames(df_completo)[which(colnames(df_completo) == "Rango etario")] <- "Rango"
n <- ncol(df_completo) - ncol(MetadataB)


df_completo$AfeccionesPiel

resumen <- df_completo %>%
  #group_by(Rango) %>%  # Agrupar por rango etario
  #group_by(Tipodepiel) %>%  # Agrupar por rango etario
  group_by(AfeccionesPiel) %>%
  summarise(
    across(2: (n+1),  # Seleccionar las primeras 8 columnas (las subespecies)
           list(Media = ~mean(.),
                Mediana = ~median(.),
                Q1 = ~quantile(., 0.25),
                Q3 = ~quantile(., 0.75)),
           .names = "{col}_{fn}")) %>%  # Generar nombres de columna automáticamente
  #pivot_longer(cols = -Rango,  # Convertir todas las columnas excepto 'Rango' en filas
  #pivot_longer(cols = -Tipodepiel,
  pivot_longer(cols = -AfeccionesPiel,
               names_to = c(nivel, ".value"),  # Separar el nombre de la subespecie y la métrica
               names_sep = "_")
resumen <- resumen[,1:6]
resumen <- resumen[,c(2,1,3,4,5,6)]
resumen <- as.data.frame(resumen)
resumen <- resumen[-which(is.na(resumen$Tipodepiel)),]

write.xlsx(resumen, file = sprintf("~/Daniela/Biota/PipelineBiota/data/%sAbundantes_RangoEtario_sinoutliers.xlsx", nivel), overwrite = TRUE)
