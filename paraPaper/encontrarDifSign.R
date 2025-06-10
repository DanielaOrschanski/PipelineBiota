CPM_vias = AR_genK_Bo
de_host_file = "Bo"
CPM_vias <- CPM_vias[,-which(colnames(CPM_vias) == "Promedio")]
nivel = "Genero"
library(writexl)
library(openxlsx)
library(PipelineBiota)

generateOTUsTableGrupal(patients_dir = "~/Daniela/Biota/Muestras/73m", source = "KRAKEN", conEukaryota = FALSE, de_host = "Bowtie")

#POR ESPECIE:
lista_de_hosters <- c("Bowtie", "BWA","RSubread","","Bowtie", "BWA", "RSubread", "sinDH_PD","")
lista_source <- c("KRAKEN", "KRAKEN","KRAKEN","KRAKEN", "DRAGEN", "DRAGEN","DRAGEN","DRAGEN","DRAGEN" )
lista_esp <- mapply(extraerAR_Species,
                    patients_dir = "~/Daniela/Biota/Muestras/73m",
                    source = lista_source,
                    de_host = lista_de_hosters,
                    SIMPLIFY = FALSE)
names(lista_esp) <- c( "KBo","KBWA","KRs","K","DBo","DBWA", "DRs","DdhD", "D")
lista_p_esp <- mapply(generar_p_por_categoria,
                      CPM_vias = lista_esp,
                      de_host_file = names(lista_esp),
                      MoreArgs = list(nivel ="Species"),
                      SIMPLIFY = FALSE)

save(lista_p_esp, "~/Daniela/Biota/PipelineBiota/paraPaper/lista_p_esp.RData")

#GENERO ----
lista_de_hosters <- c("Bowtie", "BWA","RSubread","","Bowtie", "BWA", "RSubread", "sinDH_PD","")
lista_source <- c("KRAKEN", "KRAKEN","KRAKEN","KRAKEN", "DRAGEN", "DRAGEN","DRAGEN","DRAGEN","DRAGEN" )
lista_gen <- mapply(extraerAR_Generos,
                    patients_dir = "~/Daniela/Biota/Muestras/73m",
                    source = lista_source,
                    de_host = lista_de_hosters,
                    SIMPLIFY = FALSE)
save(lista_gen, file = "~/Daniela/Biota/PipelineBiota/paraPaper/lista_gen_nueva.RData")
load("~/Daniela/Biota/PipelineBiota/paraPaper/lista_gen_nueva.RData")

names(lista_gen) <- c( "KBo","KBWA","KRs","K","DBo","DBWA", "DRs","DdhD", "D")
dfK <- lista_gen[[1]]
dfD <- lista_gen[[5]]
setdiff(dfK$Genus, dfD$Genus)
lista_p_gen <- mapply(generar_p_por_categoria,
                      CPM_vias = lista_gen,
                      de_host_file = names(lista_gen),
                      MoreArgs = list(nivel ="Genus"),
                      SIMPLIFY = FALSE)

save(lista_p_gen, file = "~/Daniela/Biota/PipelineBiota/paraPaper/lista_p_gen_nueva.RData")
load("~/Daniela/Biota/PipelineBiota/paraPaper/lista_p_gen_nueva.RData")


library(dplyr)
library(purrr)

#VIAS ---- ----------------------------
source("~/Daniela/Biota/PipelineBiota/R/RunHuman.R", echo=TRUE)
lista_de_hosters <- c("Bowtie", "BWA","RSubread","")
lista_vias <- mapply(function(host) {
  generatePathwaysTable(patients_dir = "~/Daniela/Biota/Muestras/73m", de_host = host)[[2]]
}, lista_de_hosters, SIMPLIFY = FALSE)

df1 <- lista_vias[[1]]

names(lista_vias) <- c( "BoH","bwaH","RsH","-H")
lista_p_vias <- mapply(generar_p_por_categoria,
                      CPM_vias = lista_vias,
                      de_host_file = names(lista_vias),
                      MoreArgs = list(nivel ="Pathway"),
                      SIMPLIFY = FALSE)

# Lista con los dataframes que salieron de genera_p_por_categoria -------

dfs <- lista_p_gen
dfs <- lista_p_vias
names(dfs)
nombres <- names(dfs)
# Función para renombrar columnas
library(dplyr)
library(purrr)

colnames(dfs[[1]])
renombrar_columnas <- function(df, nombre) {
  df %>%
    rename_with(~ paste0(.x, "_", nombre), c("Wilcoxon", "GrupoDominante"))
    #rename_with(~ paste0(.x, "_", nombre), c("Sexo","Edad", "Sexo*Edad", "GrupoDominante"))
}

# Renombrar y unir por "Categoria" y "Via"
str(dfs)
dfs_renombrados <- map2(dfs, nombres, renombrar_columnas)

df_final <- dfs_renombrados[[1]]

for (i in 2:length(dfs_renombrados)) {
  df_final <- merge(df_final, dfs_renombrados[[i]], by = c("Categoria", "Via"), all = TRUE)
}

df_final <- na.omit(df_final)

#df_final_sexo <- df_final[df_final$Categoria == "Sexo",]
#length(unique(df_final_sexo$Via))

#df_final_rango <- df_final[df_final$Categoria == "Rango etario",]
#length(unique(df_final_rango$Via))
#unique(df_final$Categoria)

# Filtrar casos donde una metodología encuentra p<0.05 y otra p>0.05
df_significativo <- df_final %>%
  filter(
    rowSums(select(., starts_with("Wilcoxon_")) < 0.05, na.rm = TRUE) > 0 &
      rowSums(select(., starts_with("Wilcoxon_")) > 0.05, na.rm = TRUE) > 0
  )

df_significativo_KvsD <- df_final %>%
  filter(
    rowSums(select(., starts_with("Wilcoxon_K")) < 0.05, na.rm = TRUE) > 0 &
      rowSums(select(., starts_with("Wilcoxon_D")) > 0.05, na.rm = TRUE) > 0
  )

#Todas las de K son signif y las de D no:
df_significativo_todasKvstodasD <- df_final %>%
  filter(
    if_all(starts_with("Wilcoxon_K"), ~ . < 0.05) &  # Todas las Wilcoxon_K < 0.05
      if_all(starts_with("Wilcoxon_D"), ~ . > 0.05)    # Todas las Wilcoxon_D > 0.05
  )
#Al reves:
df_significativo_todasDvstodasK <- df_final %>%
  filter(
    if_all(starts_with("Wilcoxon_D"), ~ . < 0.05) &  # Todas las Wilcoxon_K < 0.05
      if_all(starts_with("Wilcoxon_K"), ~ . > 0.05)    # Todas las Wilcoxon_D > 0.05
  )

#Todas signif excepto alguna de BWA:
df_significativo_todasvsBWA <- df_final %>%
  filter(
    if_any(contains("BWA") & where(is.numeric), ~ . > 0.05) &  # Al menos UNA de las columnas con "BWA" > 0.05
      if_all(starts_with("Wilcoxon_") & !contains("BWA") & where(is.numeric), ~ . < 0.05) # TODAS las "Wilcoxon_" excepto "BWA" < 0.05
  )


df_sign_sexo_todasKvstodasD <- df_significativo_todasKvstodasD[df_significativo_todasKvstodasD$Categoria == "Sexo",]
length(unique(df_sign_sexo_todasKvstodasD$Via))

df_sign_sexo_todasDvstodasK <- df_significativo_todasDvstodasK[df_significativo_todasDvstodasK$Categoria == "Sexo",]
df_sign_sexo_todasvsBWA <- df_significativo_todasvsBWA[df_significativo_todasvsBWA$Categoria == "Sexo",]
length(unique(df_sign_sexo_todasvsBWA$Via))

df_sign_sexo_KvsD <- df_significativo_KvsD[df_significativo_KvsD$Categoria == "Sexo",]
write.xlsx(df_sign_sexo_KvsD, file = "~/Daniela/Biota/PipelineBiota/paraPaper/DifSign_Genero_Sexo_KsiDno.xlsx")
length(unique(df_sign_sexo_KvsD$Via))

df_sign_sexo <- df_significativo[df_significativo$Categoria == "Sexo",]
length(unique(df_sign_sexo$Via))
write.xlsx(df_sign_sexo, file = "~/Daniela/Biota/PipelineBiota/paraPaper/DifSign_Genero_Sexo.xlsx")
write.xlsx(df_sign_sexo, file = "~/Daniela/Biota/PipelineBiota/paraPaper/DifSign_Familia_Sexo.xlsx")
write.xlsx(df_sign_sexo, file = "~/Daniela/Biota/PipelineBiota/paraPaper/DifSign_Filo_Sexo.xlsx")
write.xlsx(df_sign_sexo, file = "~/Daniela/Biota/PipelineBiota/paraPaper/DifSign_Vias_Sexo.xlsx")

#rango:
df_sign_rango <- df_significativo[df_significativo$Categoria == "Rango etario",]
length(unique(df_sign_rango$Via))
write.xlsx(df_sign_rango, file = "~/Daniela/Biota/PipelineBiota/paraPaper/DifSign_Genero_RangoEtario.xlsx")
write.xlsx(df_sign_rango, file = "~/Daniela/Biota/PipelineBiota/paraPaper/DifSign_Familia_RangoEtario.xlsx")
write.xlsx(df_sign_rango, file = "~/Daniela/Biota/PipelineBiota/paraPaper/DifSign_Filo_RangoEtario.xlsx")
write.xlsx(df_sign_rango, file = "~/Daniela/Biota/PipelineBiota/paraPaper/DifSign_Vias_RangoEtario.xlsx")

df_sign_rango_todasDvstodasK <- df_significativo_todasDvstodasK[df_significativo_todasDvstodasK$Categoria == "Rango etario",]
length(unique(df_sign_rango_todasDvstodasK$Via))

df_sign_rango_todasKvstodasD <- df_significativo_todasKvstodasD[df_significativo_todasKvstodasD$Categoria == "Rango etario",]
length(unique(df_sign_rango_todasKvstodasD$Via))

df_sign_rango_todasvsBWA <- df_significativo_todasvsBWA[df_significativo_todasvsBWA$Categoria == "Rango etario",]
length(unique(df_sign_rango_todasvsBWA$Via))

df_sign_rango_KvsD <- df_significativo_KvsD[df_significativo_KvsD$Categoria == "Rango etario",]
length(unique(df_sign_rango_KvsD$Via))

"" %in% df_sign_rango_todasvsBWA$Via

#Hacer como un upset pero con las especies que dieron dif sign entre sexo:

data_list <- list(
  BoK = unique(df_sign_sexo$Via[df_sign_sexo$Wilcoxon_KBo < 0.05]),
  RsK = unique(df_sign_sexo$Via[df_sign_sexo$Wilcoxon_KRs < 0.05]),
  bwaK = unique(df_sign_sexo$Via[df_sign_sexo$Wilcoxon_KBWA < 0.05]),
  '-K' = unique(df_sign_sexo$Via[df_sign_sexo$Wilcoxon_K < 0.05]),

  BoD = unique(df_sign_sexo$Via[df_sign_sexo$Wilcoxon_DBo < 0.05]),
  RsD = unique(df_sign_sexo$Via[df_sign_sexo$Wilcoxon_DRs < 0.05]),
  bwaD = unique(df_sign_sexo$Via[df_sign_sexo$Wilcoxon_DBWA < 0.05]),
  '-D' = unique(df_sign_sexo$Via[df_sign_sexo$Wilcoxon_D < 0.05]),
  DdhD = unique(df_sign_sexo$Via[df_sign_sexo$Wilcoxon_DdhD< 0.05])
)


colnames(df_sign_rango)

data_list <- list(
  BoK = unique(df_sign_rango$Via[df_sign_rango$Wilcoxon_KBo < 0.05]),
  RsK = unique(df_sign_rango$Via[df_sign_rango$Wilcoxon_KRs < 0.05]),
  bwaK = unique(df_sign_rango$Via[df_sign_rango$Wilcoxon_KBWA < 0.05]),
  '-K' = unique(df_sign_rango$Via[df_sign_rango$Wilcoxon_K < 0.05]),

  BoD = unique(df_sign_rango$Via[df_sign_rango$Wilcoxon_DBo < 0.05]),
  RsD = unique(df_sign_rango$Via[df_sign_rango$Wilcoxon_DRs < 0.05]),
  bwaD = unique(df_sign_rango$Via[df_sign_rango$Wilcoxon_DBWA < 0.05]),
  '-D' = unique(df_sign_rango$Via[df_sign_rango$Wilcoxon_D < 0.05]),
  DdhD = unique(df_sign_rango$Via[df_sign_rango$Wilcoxon_DdhD< 0.05])
)

#vias: sexo
data_list <- list(
  BoH = unique(df_sign_sexo$Via[df_sign_sexo$Wilcoxon_BoH < 0.05]),
  RsH = unique(df_sign_sexo$Via[df_sign_sexo$Wilcoxon_RsH < 0.05]),
  bwaH = unique(df_sign_sexo$Via[df_sign_sexo$Wilcoxon_bwaH < 0.05]),
  '-H' = unique(df_sign_sexo$Via[df_sign_sexo$`Wilcoxon_-H` < 0.05])
)
#vias:rango
data_list <- list(
  BoH = unique(df_sign_rango$Via[df_sign_rango$Wilcoxon_BoH < 0.05]),
  RsH = unique(df_sign_rango$Via[df_sign_rango$Wilcoxon_RsH < 0.05]),
  bwaH = unique(df_sign_rango$Via[df_sign_rango$Wilcoxon_bwaH < 0.05]),
  '-H' = unique(df_sign_rango$Via[df_sign_rango$`Wilcoxon_-H` < 0.05])
)

# El upset plot no es para
library(UpSetR)
library(ComplexHeatmap)
data_matrix <- fromList(data_list)

#sexo:
bar_colors <- rep("grey40", 45)
bar_colors[c(9, 10, 16)] <- "black"  # Colores para destacar
#bar_colors[c(9, 11, 18)] <- "darkred"  # Colores para destacar
bar_colors[c(1,3)] <- "grey15"  # Colores para destacar

#rango:
bar_colors <- rep("grey40", 47)
bar_colors[c(5,11,14)] <- "black"  # Colores para destacar
bar_colors[c(1,2)] <- "grey15"  # Colores para destacar


upset(
  data_matrix,
  sets = names(data_list),
  order.by = "freq",
  #mainbar.y.label = "Nº Genera signif dif between sexes",
  mainbar.y.label = "Nº Genera signif dif between age ranges",
  sets.x.label = "Nº Signif Genera Detected",
  text.scale = c(1.5, 1.3, 1.2, 1, 1.3, 1.3),  # Aumenta la legibilidad del texto
  point.size = 2,  # Puntos más grandes en la matriz
  line.size = 0,  # Líneas más gruesas para mayor visibilidad
  keep.order = TRUE,  # Mantiene el orden original de los conjuntos
  sets.bar.color = "grey50",  # Todas las barras de conjuntos en negro
  matrix.color = "black",  # Matriz en negro
  main.bar.color = bar_colors,
  shade.color = "gray90",  # Fondo de sombreado más claro para contraste
  number.angles = 0,  # Mejor legibilidad de números
  mb.ratio = c(0.7, 0.3)
)


##############################
#Grafico con puntos de tamaño segun el p-value y color segun sexo dominante:

#sexo:
gen_sign_sexo <- c("Propionibacterium", "Elizabethkingia", "Rhodopseudomonas", "Dermabacter", "Enterococcus",
                                  "Aeromonas", "Flavobacterium", "Gemella")

unique(df_sign_sexo$Via)
unique(df_sign_rango$Via)
gen_sign_sexo <- c("ARGININE-SYN4-PWY","RIBOSYN2-PWY",
                   "PWY-6922", "SER-GLYSYN-PWY", "P4-PWY",
                   "THRESYN-PWY", "PWY-7761")

#gen_sign_sexo <- unique(df_sign_sexo$Via)
df_gen_sign_sexo <- df_sign_sexo[df_sign_sexo$Via %in% gen_sign_sexo,]
df_gen_sign_sexo <- df_gen_sign_sexo[,-1]

#rango:
gen_sign_rango <- c("HEMESYN2-PWY", "PWY-5345", "PWY4FS-7", "PWY4FS-8",
                    "PWY-6385", "PWY-6386", "PWY-8187",
                    "THRESYN-PWY", "PWY-7220", "PWY-7222",
                    "P41-PWY", "TCA-GLYOX-BYPASS"
                    )
#gen_sign_rango <- unique(df_sign_rango$Via)
df_gen_sign_rango <- df_sign_rango[df_sign_rango$Via %in% gen_sign_rango,]
df_gen_sign_rango <- df_gen_sign_rango[,-1]

df_gen_sign_ambos <- rbind(df_gen_sign_rango, df_gen_sign_sexo)

#PAra paper: no necesario para analsiis ---
df_gen_sign_sexo <- df_gen_sign_rango

colnames(df_gen_sign_sexo)
df <- df_gen_sign_sexo[, c(1,2,3, 4,6,8)]
df[,-1] <- round(df[,-1], 3)
CPM_Vias_Clases_Anotacion_DifSign_82p <- read_excel("~/Daniela/Biota/CPM_Vias_Clases_Anotacion_DifSign_82p.xlsx")
colnames(CPM_Vias_Clases_Anotacion_DifSign_82p)
df <- merge(df, CPM_Vias_Clases_Anotacion_DifSign_82p[, c("Via", "Clase", "Description")], by = "Via")
colnames(df)
df <- df[, c(1,2,5,4,3,6,7)]
write.xlsx(df, file = "~/Daniela/Biota/PipelineBiota/paraPaper/Vias_Sign_Rango.xlsx")
library(ggplot2)
library(dplyr)
library(tidyr)

#-------
# Convertir a formato largo

df_long <- df_gen_sign_ambos %>%
  pivot_longer(cols = -c(Via, Categoria),
               names_to = c(".value", "Method"),
               names_sep = "_")

df_long <- df_gen_sign_sexo %>%
  pivot_longer(cols = -Via,
               names_to = c(".value", "Method"),
               names_sep = "_")

df_long <- df_gen_sign_rango %>%
  pivot_longer(cols = -Via,
               names_to = c(".value", "Method"),
               names_sep = "_")

# Definir colores para los grupos dominantes
color_palette <- c("Femenino" = "grey", "Masculino" = "black")
color_palette <- c("18-35" = "black", "35-55" = "grey30", ">55" = "grey70")


df_long$pValor <- ifelse(df_long$Wilcoxon<0.05, -log10(df_long$Wilcoxon), NA)
unique(df_long$Method)

df_long$Method[df_long$Method == "KBo"] <- "BoK"
df_long$Method[df_long$Method == "KBWA"] <- "bwaK"
df_long$Method[df_long$Method == "KRs"] <- "RsK"
df_long$Method[df_long$Method == "K"] <- "-K"
df_long$Method[df_long$Method == "DBo"] <- "BoD"
df_long$Method[df_long$Method == "DBWA"] <- "bwaD"
df_long$Method[df_long$Method == "DRs"] <- "RsD"
df_long$Method[df_long$Method == "D"] <- "-D"


ggplot(df_long, aes(x = Method,
                    y = Via,
                    size = pValor,
                    fill = GrupoDominante
                    )) +
  geom_point(color = "black", shape = 21, alpha = 0.8, stroke = 0.5) +
  scale_size_continuous(
    name = "p-value",
    breaks = c(-log10(0.049), -log10(0.01), -log10(0.0015)),  # Puntos clave de significancia
    labels = c("0.05", "0.01", "0.001"),  # Mostrar valores originales de Wilcoxon
    range = c(2, 12)  # Controla el tamaño de los puntos
  ) +
  scale_y_discrete(limits = rev(c("Propionibacterium",
                               "Dermabacter","Rhodopseudomonas",
                              "Enterococcus", "Elizabethkingia",
                              "Aeromonas", "Flavobacterium", "Gemella"))) +
  scale_x_discrete(
    labels = c("-D", "DdhD", "BoD", "RsD", "bwaD", " ", "-K", "BoK", "RsK", "bwaK"),  # “ ” es un espacio unicode más grande
    limits = c("-D", "DdhD", "BoD", "RsD", "bwaD", " ", "-K", "BoK", "RsK", "bwaK")
  ) +

  scale_fill_manual(name = "Dominant Group",
                    values = c("Femenino" = "grey50", "Masculino" = "black"),
                    labels = c("Femenino" = "Female", "Masculino" = "Male")) +
  theme_minimal() +
  labs(x = "Methodology", y = "Genus", title = "") +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 11),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),      # Tamaño del texto de las categorías
    legend.title = element_text(size = 13)
  ) + guides(
    fill = guide_legend(
      override.aes = list(size = 5)
    )
  )

#Para vias - sexo:
summary(df_long$Wilcoxon)
ggplot(df_long, aes(x = Method,
                    y = Via,
                    size = pValor,  # Controla el tamaño con -log10(Wilcoxon)
                    fill = GrupoDominante,
                    shape = GrupoDominante)) +  # Forma según GrupoDominante
  geom_point(color = "black", alpha = 0.8) +
  scale_size_continuous(
    name = "p-value",
    breaks = c(-log10(0.043), -log10(0.01), -log10(0.0035)),  # Puntos clave de significancia
    labels = c("0.05", "0.01", "0.001"),  # Mostrar valores originales de Wilcoxon
    range = c(2, 12)  # Controla el tamaño de los puntos
  ) +
  scale_y_discrete(limits = rev(c("PWY-6922", "ARGININE-SYN4-PWY","SER-GLYSYN-PWY",
                                  "RIBOSYN2-PWY", "P4-PWY",
                                  "THRESYN-PWY","PWY-7761" )
  )) +
  scale_x_discrete(limits = (c("-H", "BoH", "RsH" ,"bwaH"))) +
  scale_fill_manual(values = color_palette,
                    #guide = "none",
                    name = "Sex",
                    labels = c("Femenino" = "Female", "Masculino" = "Male")) +
  scale_shape_manual(name = "Sex",
                     values = c("Femenino" = 21, "Masculino" = 21),
                     labels = c("Femenino" = "Female", "Masculino" = "Male"),
                     guide = "none") +
  theme_minimal() +
  labs(x = "Methodology", y = "Pathways",
       fill = "Dominant Group",
       title = "") +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size=10),
    legend.position = "bottom",
    legend.box = "vertical"
  )+ guides(
    fill = guide_legend(
      override.aes = list(size = 5)
    )
  )


#PAra vias - ambos:
summary(df_long$Wilcoxon)
df_long$GrupoDominante[df_long$GrupoDominante == "35-55"] <- "18-35"
unique(df_long$GrupoDominante)
df_long$GrupoDominante <- factor(df_long$GrupoDominante, levels = c("18-35", ">55", "Masculino"))

ggplot(df_long, aes(x = Method,
                    y = Via,
                    size = pValor,  # Controla el tamaño con -log10(Wilcoxon)
                    fill = GrupoDominante)) +  # Forma según GrupoDominante
  geom_point(shape = 21, color = "black", alpha = 0.8) +
  scale_size_continuous(
    name = "p-value",
    breaks = c(-log10(0.043), -log10(0.01), -log10(0.0035)),  # Puntos clave de significancia
    labels = c("0.05", "0.01", "0.001"),  # Mostrar valores originales de Wilcoxon
    range = c(2, 12)  # Controla el tamaño de los puntos
  ) +
  scale_y_discrete(limits = rev(c("PWY-6922", "ARGININE-SYN4-PWY","SER-GLYSYN-PWY",
                                  "RIBOSYN2-PWY", "P4-PWY","PWY-7761", "THRESYN-PWY",

                                  "HEMESYN2-PWY", "PWY-5345", "PWY4FS-7", "PWY4FS-8",
                                  "PWY-6385", "PWY-6386", "PWY-8187",
                                  "PWY-7220", "PWY-7222",
                                  "P41-PWY", "TCA-GLYOX-BYPASS")
                                )) +
  scale_x_discrete(limits = (c("-H", "BoH", "RsH" ,"bwaH"))) +
  scale_fill_manual(
    name = "Most Abundant Group",
    values = c(">55" = "grey50", "18-35" = "black", "35-55" = "black", "Masculino" = "white"),
    labels = c("Femenino" = "Female", "Masculino" = "Male",
               ">55" = "Older", "18-35" = "Younger", "35-55" = "Younger")) +
  theme_minimal() +
  labs(x = "Methodology", y = "Pathways", title = "") +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size=10)
  )+ guides(
    fill = guide_legend(
      override.aes = list(size = 5)
    )
  )

#para vias - amnbos con pattern: no puedohacer que ande!

#Para vias - RANGO:
df_long$GrupoDominante[df_long$GrupoDominante == "35-55"] <- "18-35"
df_long$GrupoDominante <- factor(df_long$GrupoDominante, levels = c("18-35", ">55"))

ggplot(df_long, aes(x = Method,
                    y = Via,
                    size = pValor,  # Controla el tamaño con -log10(Wilcoxon)
                    fill = GrupoDominante)) +  # Forma según GrupoDominante
  geom_point( shape = 21, color = "black", alpha = 0.8) +
  scale_size_continuous(
    name = "p-value",
    breaks = c(-log10(0.043), -log10(0.01), -log10(0.0035)),  # Puntos clave de significancia
    labels = c("0.05", "0.01", "0.001"),  # Mostrar valores originales de Wilcoxon
    range = c(2, 10)  # Controla el tamaño de los puntos
  ) +
  scale_fill_manual(
    name = "Most Abundant Group",
    values = c("18-35" = "black", ">55" = "grey50"),
    labels = c("18-35" = "Younger",  ">55" = "Older")
  ) +
  scale_y_discrete(limits = rev(c( "HEMESYN2-PWY", "PWY-5345", "PWY4FS-7", "PWY4FS-8",
                                   "PWY-6385", "PWY-6386", "PWY-8187",
                                   "THRESYN-PWY", "PWY-7220", "PWY-7222",
                                   "P41-PWY", "TCA-GLYOX-BYPASS" ))) +
  scale_x_discrete(limits = (c("-H", "BoH", "RsH" ,"bwaH"))) +
  theme_minimal() +
  labs(x = "Methodology", y = "Pathways", title = "") +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size=10),
    legend.position = "bottom",
    legend.box = "vertical"
  )+ guides(
    fill = guide_legend(
      override.aes = list(size = 5)
    )
  )


#Vias: -----------------------

df1 <- lista_vias[[1]]
df_combined <- bind_rows(lapply(names(lista_vias), function(nombre) {
  lista_vias[[nombre]] %>%
    mutate(Metodologia = nombre)
}), .id = "ID")
df_combined <- df_combined[,-1]

  #sexo:
vias_select <- df_sign_sexo$Via
vias_select <- c("Propionibacterium", "Elizabethkingia", "Rhodopseudomonas", "Dermabacter", "Enterococcus",
                 "Aeromonas", "Flavobacterium", "Gemella")
vias_select <- vias_select[c(1,9, 11, 12, 14, 16)]
  #rango:
vias_select <- df_sign_rango$Via
vias_select <- vias_select[c(2,4,5,6,9,11,13,14,15,16,17,19,23, 24,25, 28)]


df_filtered <- df_combined[df_combined$Pathway %in% vias_select,]
#df_filtered <- df_combined[df_combined$Genus %in% vias_select,]
colnames(df_filtered)[1] <- "Pathway"

df_long <- df_filtered %>%
    pivot_longer(cols = -c(Pathway, Clase, Description, Metodologia),
                 names_to = "ID",
                 values_to = "Abundancia")
library(tidyverse)
#df_long <- df_filtered %>%
#  pivot_longer(cols = -c(Pathway, Metodologia),
#               names_to = "ID",
#               values_to = "Abundancia")

  # Paso 2: Unir df_long con MetadataB para obtener la información de sexo y metodología
  MetadataB <- as.data.frame(read_excel("~/Daniela/Biota/PipelineBiota/data/Metadata-soloColumnasUsables.xlsx"))
  df_merged <- df_long %>%
    left_join(MetadataB[, c("ID", "Sexo", "Rango etario")], by = "ID")  # Asegúrate de que MetadataB tenga una columna "Muestra"

  str(df_merged)
  df_merged$Sexo <- as.factor(df_merged$Sexo)
  df_merged$`Rango etario` <- factor( df_merged$`Rango etario` , levels = c("18-35", "35-55", ">55"))
  df_merged$Metodologia <- as.factor(df_merged$Metodologia)

  # Paso 3: Crear el boxplot - SEXO
  ggplot(df_merged, aes(x =  Metodologia, y = Abundancia, fill = Sexo)) +
    geom_boxplot() +
    labs(title = sprintf(""),
         x = "Sex and Methodology",
         y = "Relative Abundance") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_wrap(~ Pathway, scales = "free_y")

  p_values_df <- df_merged %>%
    group_by(Pathway, Metodologia) %>%
    summarise(p_value = wilcox.test(Abundancia ~ Sexo)$p.value) %>%
    mutate(p_label = sprintf("p = %.3f", p_value))  # Formatear el p-valor
  p_values_df$simbolo <- ifelse(p_values_df$p_value < 0.05, ifelse(p_values_df$p_value < 0.01,"**", "*"), "-")


  # Paso 3: Crear el boxplot - RANGO
  ggplot(df_merged, aes(x =  Metodologia, y = Abundancia, fill = `Rango etario`)) +
    geom_boxplot() +
    labs(title = sprintf(""),
         x = "Ageing and Methodology",
         y = "Relative Abundance") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_wrap(~ Pathway, scales = "free_y")

  p_values_df <- df_merged %>%
    group_by(Pathway, Metodologia) %>%
    summarise(p_value = kruskal.test(Abundancia ~ `Rango etario`)$p.value) %>%
    mutate(p_label = sprintf("p = %.3f", p_value))  # Formatear el p-valor
  p_values_df$simbolo <- ifelse(p_values_df$p_value < 0.05, ifelse(p_values_df$p_value < 0.01,"**", "*"), "-")


  # Calcular el valor máximo de Abundancia por Pathway
  max_abundance <- df_merged %>%
    group_by(Pathway) %>%
    summarise(max_value = max(Abundancia, na.rm = TRUE))

  p_values_df <- p_values_df %>%
    left_join(max_abundance, by = "Pathway") %>%
    mutate(y_position = max_value * 1.05)  # Ajustar la posición del p-valor

  # Graficar y añadir los p-valores al gráfico - SEXO:
  ggplot(df_merged, aes(x = Metodologia, y = Abundancia, fill = Sexo)) +
    geom_boxplot() +
    labs(title = "",
         x = "Sex and Methodology",
         y = "Relative Abundance",
         fill = "Sex") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Título eje x con ángulo
          axis.text = element_text(size = 8),  # Reducir tamaño de texto de los ejes
          axis.title = element_text(size = 10),  # Reducir tamaño de títulos de ejes
          plot.title = element_text(size = 10),  # Reducir tamaño del título principal
          legend.title = element_text(size = 8),  # Reducir tamaño del título de la leyenda
          legend.text = element_text(size = 8),  # Reducir tamaño de las etiquetas de la leyenda
          strip.text = element_text(size = 8),  # Reducir tamaño de los textos de las facetas
          legend.position = "right") +
    facet_wrap(~ Pathway, scales = "free_y", ncol = 4) +
    scale_fill_discrete(labels = c("Femenino" = "Female","Masculino" = "Male")) +
    geom_text(data = p_values_df,
              aes(x = Metodologia,
                  y = y_position,
                  #label = p_label),
                  label = simbolo),
              size = 4, color = "black", inherit.aes = FALSE)

  # Graficar y añadir los p-valores al gráfico - Rango:
  ggplot(df_merged, aes(x = Metodologia, y = Abundancia, fill = `Rango etario`)) +
    geom_boxplot() +
    labs(title = "",
         x = "Ageing and Methodology",
         y = "Relative Abundance",
         fill = "Age Ranges") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Título eje x con ángulo
          axis.text = element_text(size = 7),  # Reducir tamaño de texto de los ejes
          axis.title = element_text(size = 9),  # Reducir tamaño de títulos de ejes
          plot.title = element_text(size = 9),  # Reducir tamaño del título principal
          legend.title = element_text(size = 7),  # Reducir tamaño del título de la leyenda
          legend.text = element_text(size = 7),  # Reducir tamaño de las etiquetas de la leyenda
          strip.text = element_text(size = 7),  # Reducir tamaño de los textos de las facetas
          legend.position = "right") +
    facet_wrap(~ Pathway, scales = "free_y", ncol = 4) +
    geom_text(data = p_values_df,
              aes(x = Metodologia,
                  y = y_position,
                  #label = p_label),
                  label = simbolo),
              size = 4, color = "black", inherit.aes = FALSE)

#Para tabla paper:
p_values_wide <- p_values_df[, -which(colnames(p_values_df) %in% c("p_label", "simbolo", "max_value", "y_position"))] %>%
    pivot_wider(names_from = Metodologia,
                values_from = p_value,
                names_sort = TRUE)
CPM_Vias_Clases_Anotacion_DifSign_82p <- read_excel("~/Daniela/Biota/CPM_Vias_Clases_Anotacion_DifSign_82p.xlsx")
colnames(CPM_Vias_Clases_Anotacion_DifSign_82p)[1] <- "Pathway"
p_values_wide <- merge(p_values_wide, CPM_Vias_Clases_Anotacion_DifSign_82p[, c("Pathway", "Clase", "Description")], by = "Pathway")
p_values_wide[, c(2:5)] <- round(p_values_wide[,c(2:5)], 3)
p_values_wide$Pathway

write.xlsx(p_values_wide, file = "~/Daniela/Biota/PipelineBiota/paraPaper/p_values_wide_Rango_Vias.xlsx")

#---------
generar_p_por_categoria <- function(CPM_vias, de_host_file, nivel) {

  CPM_viasT <- as.data.frame(t(CPM_vias))
  colnames(CPM_viasT) <- CPM_viasT[1,]
  CPM_viasT <- CPM_viasT[-1,]

  #Agrego metadata:
  CPM_viasT <- cbind("ID" = rownames(CPM_viasT), CPM_viasT)

  MetadataB <- as.data.frame(read_excel("~/Daniela/Biota/PipelineBiota/data/Metadata-soloColumnasUsables.xlsx"))
  MetadataB <- MetadataB[which(MetadataB$ID %in% CPM_viasT$ID),]
  all(CPM_viasT$ID %in% MetadataB$ID)
  all(CPM_viasT$ID == MetadataB$ID)

  colnames(MetadataB)
  CPM_vias_completo <- merge(CPM_viasT, MetadataB[, c("ID", "Sexo", "Rango etario")], by = "ID")
  str(CPM_vias_completo)
  n <- ncol(CPM_viasT)
  colnames(CPM_vias_completo)[n]
  colnames(CPM_vias_completo)[n+1]
  CPM_vias_completo[,2:n] <- lapply(CPM_vias_completo[,2:n], as.numeric)
  CPM_vias_completo[,(n+1):ncol(CPM_vias_completo)] <- lapply(CPM_vias_completo[,(n+1):ncol(CPM_vias_completo)], as.factor)

  #write.xlsx(CPM_vias_completo, file =  sprintf("~/Daniela/Biota/Muestras/SubsetPathways/%s/%s_%svias_completo_AR_83p.xlsx", de_host_file, de_host_file, length(2:n)))


  #Calculo pvalues:
  p_por_categoria <- data.frame("Categoria" = c(), "Via" = c(), "Wilcoxon"= c(), "GrupoDominante" = c())

  categorias <- colnames(CPM_vias_completo)[(n+1):ncol(CPM_vias_completo)]
  vias <- colnames(CPM_vias_completo)[2:n]
  c= n+1
  j = 2
  i=1

  #c= n+1
  for (c in (n+1):ncol(CPM_vias_completo)) {
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

      if( categoria == "Maquillaje_Base" | categoria == "CuándoMaquillaje") {
        df_sin_na <- df_sin_na[-which(df_sin_na$Sexo == "Masculino"),]
      }

      #if( categoria == "AfeccionesPiel") {
      #  df_sin_na <- df_sin_na[-which(df_sin_na$`Rango etario` != "18-35"),]
      #}

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

      p_por_categoria[i, "GrupoDominante"] <-  ifelse(length(grupo_dominante) == 0, NA, as.character(grupo_dominante))

      i <- i+1
    }
  }

  unique(p_por_categoria$Categoria)
  p_por_categoria <- p_por_categoria[p_por_categoria$Categoria %in% c("Sexo", "Rango etario",
                                                                      "Tipodepiel", "Maquillaje_Base",
                                                                      "FacilidadBroncearse", "ActividadFísica",
                                                                      "CuandoMaquillaje"), ]
  write.xlsx(p_por_categoria, file =  sprintf("~/Daniela/Biota/PipelineBiota/paraPaper/%s_%s_p_por_categoria_83p.xlsx", nivel, de_host_file))
  return(p_por_categoria)

}

###########################################################################################

# Calcular p-values con ANOVA de dos vías
generar_p_por_categoria_ANOVA <- function(CPM_vias, de_host_file, nivel) {

  CPM_viasT <- as.data.frame(t(CPM_vias))
  colnames(CPM_viasT) <- CPM_viasT[1,]
  CPM_viasT <- CPM_viasT[-1,]

  #Agrego metadata:
  CPM_viasT <- cbind("ID" = rownames(CPM_viasT), CPM_viasT)

  MetadataB <- as.data.frame(read_excel("~/Daniela/Biota/PipelineBiota/data/Metadata-soloColumnasUsables.xlsx"))
  MetadataB <- MetadataB[which(MetadataB$ID %in% CPM_viasT$ID),]
  all(CPM_viasT$ID %in% MetadataB$ID)
  all(CPM_viasT$ID == MetadataB$ID)

  colnames(MetadataB)
  CPM_vias_completo <- merge(CPM_viasT, MetadataB[, c("ID", "Sexo", "Rango etario")], by = "ID")
  str(CPM_vias_completo)
  n <- ncol(CPM_viasT)
  colnames(CPM_vias_completo)[n]
  colnames(CPM_vias_completo)[n+1]
  CPM_vias_completo[,2:n] <- lapply(CPM_vias_completo[,2:n], as.numeric)
  CPM_vias_completo[,(n+1):ncol(CPM_vias_completo)] <- lapply(CPM_vias_completo[,(n+1):ncol(CPM_vias_completo)], as.factor)

  p_por_categoria_ANOVA <- data.frame("Categoria" = c(), "Via" = c(), "Sexo"= c(), "Edad" = c(), "Sexo*Edad" = c(), "GrupoDominante" = c())

  # Iterar sobre las categorías y vías
  categorias <- colnames(CPM_vias_completo)[(n+1):ncol(CPM_vias_completo)]
  vias <- colnames(CPM_vias_completo)[2:n]
  i <- 1

  for (c in (n+2):ncol(CPM_vias_completo)) {
    categoria <- colnames(CPM_vias_completo)[c]
    print(categoria)

    for (j in 2:(n)) {

      via <- colnames(CPM_vias_completo)[j]
      print(via)
      df_sin_na <- CPM_vias_completo[complete.cases(CPM_vias_completo[, categoria]), ]

      # Filtrar categorías si contienen "Desconocido" o NA
      if (any(df_sin_na[,categoria] == "Desconocido")) {
        df_sin_na <- df_sin_na[-which(df_sin_na[,categoria] == "Desconocido"),]
      }

      if (any(is.na(df_sin_na[,categoria]))) {
        df_sin_na <- df_sin_na[-which(is.na(df_sin_na[,categoria])),]
      }

      # Filtrar si la categoría es "Maquillaje_Base" o "CuándoMaquillaje" para "Sexo" == "Masculino"
      if( categoria == "Maquillaje_Base" | categoria == "CuándoMaquillaje") {
        df_sin_na <- df_sin_na[-which(df_sin_na$Sexo == "Masculino"),]
      }

      # Calcular la media para cada grupo en la categoría para sacar el grupo dominante:
      mean_values <- aggregate(df_sin_na[, via], by = list(df_sin_na[, categoria]), FUN = mean)
      colnames(mean_values) <- c("Grupo", "Media")
      grupo_dominante <- mean_values$Grupo[which.max(mean_values$Media)]
      print(grupo_dominante)

      # Realizar el ANOVA de dos vías (Sexo, Edad, y la interacción Sexo*Edad)
      aov_result <- aov(df_sin_na[, via] ~ df_sin_na$Sexo * df_sin_na$`Rango etario`)

      # Extraer los valores p de los tres factores
      p_sex <- summary(aov_result)[[1]][1, "Pr(>F)"]
      p_age <- summary(aov_result)[[1]][2, "Pr(>F)"]
      p_interaction <- summary(aov_result)[[1]][3, "Pr(>F)"]

      # Almacenar los resultados
      p_por_categoria_ANOVA[i, "Categoria"] <- categoria
      p_por_categoria_ANOVA[i, "Via"] <- via
      p_por_categoria_ANOVA[i, "Sexo"] <- p_sex
      p_por_categoria_ANOVA[i, "Edad"] <- p_age
      p_por_categoria_ANOVA[i, "Sexo*Edad"] <- p_interaction
      p_por_categoria_ANOVA[i, "GrupoDominante"] <- ifelse(length(grupo_dominante) == 0, NA, as.character(grupo_dominante))

      i <- i + 1
    }
  }

  unique(p_por_categoria_ANOVA$Categoria)
  write.xlsx(p_por_categoria_ANOVA, file =  sprintf("~/Daniela/Biota/PipelineBiota/paraPaper/%s_%s_p_por_categoria_ANOVA_83p.xlsx", nivel, de_host_file))
  return(p_por_categoria_ANOVA)

}


###############################################################################


generar_graficos <- function(p_significativos, de_host_file) {

  CPM_vias_completo <- as.data.frame(read_excel(sprintf("~/Daniela/Biota/Muestras/SubsetPathways/%s/%s_239vias_completo_AR_83p.xlsx", de_host_file, de_host_file)))
  CPM_vias_completo <- CPM_vias_completo[1:83,]

  n <- which(colnames(CPM_vias_completo) == "FechadeNacimiento")

  unique(p_significativos$Categoria)
  categoria = "Sexo"
  categoria = "Rango etario"

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

    df_sin_na <- CPM_vias_completo
    #df_sin_na <- CPM_vias_completo[complete.cases(CPM_vias_completo[,categoria]), ]
    #df_sin_na <- df_sin_na[, order(colnames(df_sin_na))]
    #df_sin_na <- df_sin_na[, unique(colnames(df_sin_na))]

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
    #if( categoria == "FacilidadBroncearse") {
    #df_sin_na <- df_sin_na[-which(df_sin_na$ID %in% ids_out),]
    #}

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

    if(categoria == "Sexo") {
      df_sin_na[,categoria] <- as.factor(df_sin_na[,categoria])
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
  grid.arrange(grobs = lista_graficos, ncol = 5)

}


# BOXPLOTS PARA COMPARAR LAS DISTRIBUCIONES DE LA AR DE LOS GENEROS QUE DAN DIF SIGN CON UNA METODOLOGIA PERO NO CON OTRA:

# Identificar columnas que contienen "_K" y "_D"
cols_K <- grep("Wilcoxon_K", colnames(df_sign_sexo), value = TRUE)
cols_D <- grep("Wilcoxon_D", colnames(df_sign_sexo), value = TRUE)

# Filtrar las filas donde:
# - Todos los valores en columnas "_K" sean < 0.05
# - Todos los valores en columnas "_D" sean > 0.05
Ksi_Dno <- df_sign_sexo %>%
  filter(if_all(all_of(cols_K), ~ . < 0.05) & if_all(all_of(cols_D), ~ . > 0.05))
Dsi_Kno <- df_sign_sexo %>%
  filter(if_all(all_of(cols_K), ~ . > 0.05) & if_all(all_of(cols_D), ~ . < 0.05))


"Helicobacter"
genero_select = "Corynebacterium"
genero_select = "Staphylococcus"
genero_select = "Neisseria"
genero_select = "Lactobacillus"

genero_select = "Chlamydia"
genero_select = "Rhodopseudomonas"
genero_select = "Actinobacillus"
genero_select = "Aeromonas"
genero_select = "Burkholderia"
genero_select = "Klebsiella"
genero_select = "Serratia"
genero_select = "Yersinia"

list_generos <- c( "Chlamydia", "Rhodopseudomonas", "Actinobacillus",
                   "Burkholderia","Serratia", "Yersinia")

list_generos <- Ksi_Dno$Via

lista_df <- list("KBWA" = AR_genK_BWA, "KBo" = AR_genK_Bo, "KRs" = AR_genK_Rs, "K" = AR_genK_sin,
                 "DBWA" = AR_genD_BWA, "DBo" = AR_genD_Bo, "DRs" = AR_genD_Rs, "DdhD" = AR_genD_sin,
                 "D" = ARgenD_sinDH_PD)
df_combined <- bind_rows(lapply(names(lista_df), function(nombre) {
  lista_df[[nombre]] %>%
    mutate(Metodologia = nombre)
}), .id = "ID")
df_combined <- df_combined[,-1]
df_filtered <- df_combined
df_filtered$Source <- ifelse(grepl("K", df_filtered$Metodologia), "KRAKEN", "DRAGEN")
df_filtered$Metodologia <- factor(df_filtered$Metodologia, levels = names(lista_df))

genero_select <- list_generos[1]
list_plots <- list()
for (g in list_generos ){
  plot <- boxplotPorMetodologia_porSexo_porGenero(df_filtered = df_filtered, genero_select =  g)
  list_plots <- append(list_plots, list(plot))
}

library(gridExtra)
grid.arrange(grobs = list_plots, ncol= 4)

#Achicar los tamaños de las letras:
list_plots_red <- lapply(list_plots, function(plot) {
  plot +
    theme(
      text = element_text(size = 8),
      axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 6),
      plot.title = element_text(size = 8, face = "bold"),
      legend.position = "none"  # Oculta la leyenda
    )
})

grid.arrange(grobs = list_plots_red[1:9], ncol = 3)

boxplotPorMetodologia_porSexo_porGenero <- function(df_filtered, genero_select) {
  df_long <- df_filtered %>%
    pivot_longer(cols = -c(Genus, Promedio, Metodologia, Source),
                 names_to = "ID",
                 values_to = "Abundancia") %>%
    filter(Genus == genero_select)  # Filtrar por el género Corynebacterium

  # Paso 2: Unir df_long con MetadataB para obtener la información de sexo y metodología
  MetadataB <- as.data.frame(read_excel("~/Daniela/Biota/PipelineBiota/data/Metadata-soloColumnasUsables.xlsx"))
  df_merged <- df_long %>%
    left_join(MetadataB[, c("ID", "Sexo", "Rango etario")], by = "ID")  # Asegúrate de que MetadataB tenga una columna "Muestra"

  # Paso 3: Crear el boxplot
  ggplot(df_merged, aes(x =  Metodologia, y = Abundancia, fill = Sexo)) +
    geom_boxplot() +
    labs(title = sprintf("AR %s por Sexo y Metodología", genero_select),
         x = "Sexo y Metodología",
         y = "Abundancia Relativa") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  str(df_merged)
  df_merged$Sexo <- as.factor(df_merged$Sexo)
  p_values <- df_merged %>%
    group_by(Metodologia) %>%
    summarise(p_value = wilcox.test(Abundancia ~ Sexo)$p.value)

  # Unir los p-valores con el dataframe original
  df_merged_with_p <- df_merged %>%
    left_join(p_values, by = "Metodologia")

  # Graficar y añadir los p-valores al gráfico
  str(df_merged_with_p)
  plot <- ggplot(df_merged_with_p, aes(x = Metodologia, y = Abundancia, fill = Sexo)) +
    geom_boxplot() +
    labs(title = sprintf("%s", genero_select),
         x = "",
         y = "Relative Abundance") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_discrete(name = "Sex") +  # Cambia el título de la leyenda
    # Agregar los p-valores como texto
    geom_text(data = p_values,
              aes(x = 1:length(p_values$Metodologia), y = max(df_merged$Abundancia) + 0.005,
                  label = sprintf("p = %.3f", p_value)),
              inherit.aes = FALSE, color = "black", size = 2, vjust = 0)

  return(plot)
}


