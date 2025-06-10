
vens <- list.files("~/Daniela/Biota/PipelineBiota")
vens <- vens[which(grepl("VennDiagram", vens))]

file.remove(vens)

list.dirs <- list.dirs("~/Daniela/Biota/Muestras/73m", recursive = FALSE)
length(list.dirs)




T_vias <- generatePathwaysTable(patients_dir = "~/Daniela/Biota/Muestras/73m", de_host = "")
T_vias_CPM <- T_vias[[1]]
T_vias_AR <- T_vias[[2]]

Bo_vias <- generatePathwaysTable(patients_dir = "~/Daniela/Biota/Muestras/73m", de_host = "Bowtie")
Bo_vias_CPM <- Bo_vias[[1]]
Bo_vias_AR <- Bo_vias[[2]]

Rs_vias <- generatePathwaysTable(patients_dir = "~/Daniela/Biota/Muestras/73m", de_host = "RSubread")
Rs_vias_CPM <- Rs_vias[[1]]
Rs_vias_AR <- Rs_vias[[2]]

Bwa_vias <- generatePathwaysTable(patients_dir = "~/Daniela/Biota/Muestras/73m", de_host = "BWA")
Bwa_vias_CPM <- Bwa_vias[[1]]
Bwa_vias_AR <- Bwa_vias[[2]]

length(unique(Bwa_vias_AR$Pathway))
#T_vias <- as.data.frame(read_excel("~/Daniela/Biota/Muestras/73m/1/Vias/1T_vias.xlsx"))

#T_vias$Method <- "T"
#Bo_vias$Method <- "Bo"

common_pathways <- intersect(unique(T_vias_CPM$Pathway), unique(Bo_vias_CPM$Pathway))
common_pathways <- intersect(unique(Rs_vias_CPM$Pathway), common_pathways)
common_pathways <- intersect(unique(Bwa_vias_CPM$Pathway), common_pathways)

#common_pathways <- intersect(unique(clasified_vias$Pathway), common_pathways)

#Diferencias en pathways identificados ----------------------------------------

library(ggvenn)

# Convertir los datos a una lista
data_list <- list(
  T = unique(T_vias_CPM$Pathway),
  Bo = unique(Bo_vias_CPM$Pathway),
  Rs = unique(Rs_vias_CPM$Pathway),
  BWA = unique(Bwa_vias_CPM$Pathway)
  #,Clasiffied = unique(clasified_vias$Pathway)
)

# Crear el diagrama de Venn
ggvenn(data_list,
       #fill_color = c("red", "blue", "pink", "green", "purple", "black"),
       fill_color = c("black", "red", "blue", "green"),
       stroke_size = 0.5,
       set_name_size = 5,
       text_size = 5)

setdiff(unique(T_vias$Pathway), unique(Bo_vias$Pathway))
setdiff(unique(T_vias$Pathway), unique(clasified_vias$Pathway))


# Convertir listas a un tibble con presencia/ausencia
vias_data <- tibble(
  Vias = unique(unlist(data_list)),
  Bo = as.integer(Vias %in% data_list$Bo),
  BWA = as.integer(Vias %in% data_list$BWA),
  Rs = as.integer(Vias %in% data_list$Rs),
  SinDH = as.integer(Vias %in% data_list$T)
)

# Agrupar por combinaciones y contar especies en cada intersección
vias_summary <- vias_data %>%
  group_by(Bo, BWA, Rs, SinDH) %>%
  summarise(Vias_List = list(Vias), Count = n(), .groups = "drop")

#Las 34 que se pierden con BWA:
bwa_pierde <- vias_summary$Vias_List[which(vias_summary$BWA == 0)][[3]]
bwa_pierde_df <- T_vias_AR[which(T_vias_AR$Pathway %in% bwa_pierde), c(1,2,3) ]

todos_menossinDH_pierden <- vias_summary$Vias_List[which(vias_summary$SinDH == 1)][[1]]
todos_menossinDH_pierden_df <- T_vias_AR[which(T_vias_AR$Pathway %in% todos_menossinDH_pierden), c(1,2,3) ]

bosi_bwano <- Bo_vias_CPM[which(Bo_vias_CPM$Pathway %in% bwa_pierde),]

#Sumar info de papaers:
library(readxl)
Vias_Clases_PAPERS_DifSign_82p <- as.data.frame(read_excel("~/Daniela/Biota/PipelineBiota/data/Vias_Clases_PAPERS_DifSign_82p.xlsx"))
colnames(Vias_Clases_PAPERS_DifSign_82p)
Vias_Clases_PAPERS_pierdebwa <-  Vias_Clases_PAPERS_DifSign_82p[which(Vias_Clases_PAPERS_DifSign_82p$Via %in% bwa_pierde), c("Via", "Clase", "Description", "ClasifiaciónSegunPaper",
                                                                                                                             "NoCoincide","CategoriaInteres", "DifSign_Via_AR", "DifSign_Via_CPM")]

#Cuantas muestras tienen esas especies que bwa pierde:
bosi_bwano$ON <- apply(bosi_bwano[,-c(1,2,3)], 1, function(row) sum(row > 0))

ggplot(bosi_bwano, aes(x = "", y = ON)) +
  #geom_boxplot() +
  geom_boxplot(fill = "turquoise") +
  labs(title = "Cant Muestras Rs contienen 34 Esp  BWA no") +
  theme_minimal() +
  theme(
    text = element_text(size = 12),        # Ajuste del tamaño del texto
    axis.text.x = element_text(size = 12), # Tamaño de texto en eje X
    axis.text.y = element_text(size = 12), # Tamaño de texto en eje Y
    plot.title = element_text(size = 10),   # Tamaño del título del gráfico
    legend.position = "none"
  )

#Como son los conteos de esas vias que bwa pierde:
bosi_bwano_long <- bosi_bwano %>%
  pivot_longer(cols = -c(1,2,3),  # Excluir la primera columna (nombre de la especie)
               names_to = "Muestra",
               values_to = "CPM")
ggplot(bosi_bwano_long, aes(x = Pathway, y = CPM, fill = Pathway)) +
  geom_boxplot() +
  labs(title = "CPM 34 Vias perdidas con BWA",
       x = "Via",
       y = "CPM") +
  theme_minimal() +
  theme(legend.position = "none") + # Eliminar la leyendakraken_Ksi_clean <- kraken_Ksi[, -ncol(kraken_Ksi)]  # Eliminar última columna 'ON'
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(bosi_bwano_long, aes(x = "", y = CPM)) +
  geom_boxplot() +
  labs(title = "CPM 34 Vias perdidas con BWA",
       x = "Via",
       y = "CPM") +
  theme_minimal() +
  theme(legend.position = "none") + # Eliminar la leyendakraken_Ksi_clean <- kraken_Ksi[, -ncol(kraken_Ksi)]  # Eliminar última columna 'ON'
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Diferencias en cantidad de vias identificadas ----------------------------------------------------------
library(dplyr)

#Comparacion de cantidad de vias unicas identificadas por cada muestra segun cada dehosting:
  # Función para calcular el conteo de vías identificadas (conteos > 0):
count_identified <- function(df) {
  counts <- colSums(df[, -(1:3)] > 0) # Ignoramos las primeras 3 columnas (Pathway, Clase, Description)
  data.frame(Sample = names(counts), Count = counts)
}

T_counts <- count_identified(T_vias_CPM)
T_counts$Method <- "T"
#length(T_vias[which(T_vias$'105' >0), '105'])
Bo_counts <- count_identified(Bo_vias_CPM)
Bo_counts$Method <- "Bo"
BWA_counts <- count_identified(Bwa_vias_CPM)
BWA_counts$Method <- "BWA"
Rs_counts <- count_identified(Rs_vias_CPM)
Rs_counts$Method <- "Rs"

# Combinar todos los resultados en un solo dataframe
datos_combinados <- rbind(T_counts, Bo_counts, BWA_counts, Rs_counts)

library(ggplot2)
ggplot(datos_combinados, aes(x = Method, y = Count, fill = Method)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) + # Añadimos el boxplot
  geom_point(aes(color = Method), position = position_dodge(width = 0)) + # Puntos de datos
  geom_line(
    aes(group = Sample), # Conecta puntos del mismo Pathway
    #position = position_dodge(width = 0.3),
    color = "black",
    linetype = "dotted"
  ) +
  theme_minimal() +
  labs(
    x = "Method",
    y = "Count",
    title = "Number of identified pathways"
  ) +
  theme(
    text = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16)
  )

###################################################################################
# Diferencias en CPM de las vias en comun -----------------------------
library(dplyr)

  # Función para transformar los datos de cada método:

transform_data <- function(df, method_name) {

  df <- df[df$Pathway %in% common_pathways, ]
  df_long <- df %>%
    pivot_longer(
      cols = -c(1,2,3), # solo las muestras
      names_to = "Sample",
      values_to = "Value"
    )
  df_long$Method <- method_name
  return(df_long)
}

# Transformar cada dataframe por separado
T_long <- transform_data(T_vias_CPM, "T")
Bo_long <- transform_data(Bo_vias_CPM, "Bo")
BWA_long <- transform_data(Bwa_vias_CPM, "BWA")
Rs_long <- transform_data(Rs_vias_CPM, "Rs")

# Combinar todos los datos en un solo dataframe
datos_combinados <- bind_rows(T_long, Bo_long, BWA_long, Rs_long)
datos_combinados <- bind_rows(Bo_long, BWA_long)

library(tidyr)
library(dplyr)

# Filtrar solo vías comunes
common_pathways <- intersect(Bo_vias_CPM$Pathway, Bwa_vias_CPM$Pathway)

Bo_vias_common <- Bo_vias_CPM %>% filter(Pathway %in% common_pathways)
BWA_vias_common <- Bwa_vias_CPM %>% filter(Pathway %in% common_pathways)

# Asegurar el mismo orden en ambas tablas (por Pathway y Sample si aplica)
Bo_vias_common <- Bo_vias_common %>% arrange(Pathway)
BWA_vias_common <- BWA_vias_common %>% arrange(Pathway)

# Calcular el promedio de los valores de CPM por cada vía en Bo_vias_common
Bo_vias_common[,-c(1,2,3)] <- lapply(Bo_vias_common[,-c(1,2,3)], as.numeric)
Bo_vias_common_avg <- Bo_vias_common %>%
  rowwise() %>%
  mutate(Average_CPM = mean(c_across(4:ncol(Bo_vias_common)), na.rm = TRUE)) %>%
  ungroup()

# Seleccionar las 50 vías más abundantes en promedio
top_50_vias <- Bo_vias_common_avg %>%
  arrange(desc(Average_CPM)) %>%
  slice(1:50)

# Filtrar ambas tablas para quedarse con las vías top 50
Bo_vias_common_top50 <- Bo_vias_common %>%
  filter(Pathway %in% top_50_vias$Pathway)

BWA_vias_common_top50 <- BWA_vias_common %>%
  filter(Pathway %in% top_50_vias$Pathway)

#AR en vez de CPM:
Bo_vias_common_top50 <- Bo_vias_common
BWA_vias_common_top50 <- BWA_vias_common


AR_vias <- prop.table(as.matrix(Bo_vias_common_top50[,-c(1,2,3)]), margin = 2) * 100
colSums(AR_vias)
Bo_vias_common <-as.data.frame(cbind("Pathway" = Bo_vias_common_top50[,1],
                                     "Clase" = Bo_vias_common_top50[,2],
                                     "Description" = Bo_vias_common_top50[,3],
                                     AR_vias))

AR_vias <- prop.table(as.matrix(BWA_vias_common_top50[,-c(1,2,3)]), margin = 2) * 100
colSums(AR_vias)
BWA_vias_common <-as.data.frame(cbind("Pathway" = BWA_vias_common_top50[,1],
                                      "Clase" = BWA_vias_common_top50[,2],
                                      "Description" = BWA_vias_common_top50[,3],
                                      AR_vias))

# Transformar las tablas al formato largo
Bo_long <- Bo_vias_common %>%
  pivot_longer(-c(Pathway, Clase, Description), names_to = "Sample", values_to = "Bo_Value")

BWA_long <- BWA_vias_common %>%
  pivot_longer(-c(Pathway, Clase, Description), names_to = "Sample", values_to = "BWA_Value")

# Unir ambas tablas en un solo dataframe
bland_altman_data <- Bo_long %>%
  inner_join(BWA_long, by = c("Pathway", "Clase","Description", "Sample"))
str(bland_altman_data)
bland_altman_data[, c(5,6)] <- lapply(bland_altman_data[, c(5,6)], as.numeric)

bland_altman_data <- bland_altman_data %>%
  mutate(
    Mean = (Bo_Value + BWA_Value) / 2,
    Difference = Bo_Value - BWA_Value
  )

# Calcular media y desviación estándar de las diferencias
Diff_Mean <- mean(bland_altman_data$Difference, na.rm = TRUE)
Diff_SD <- sd(bland_altman_data$Difference, na.rm = TRUE)
Lower_Limit <- Diff_Mean - 1.96 * Diff_SD
Upper_Limit <- Diff_Mean + 1.96 * Diff_SD

comparacion_long <- bland_altman_data

# Graficar Bland-Altman sin discriminar por via
  #Nombrando a los puntos que estan por fuera del rango: ------------------------
comparacion_long$Color <- "Dentro de los límites"
comparacion_long$Color <- ifelse(comparacion_long$Difference > Upper_Limit, "UP", comparacion_long$Color)
comparacion_long$Color <- ifelse(comparacion_long$Difference < Lower_Limit, "DOWN", comparacion_long$Color)

out <-  comparacion_long[which(comparacion_long$Color != "Dentro de los límites"),]
out_freq <- as.data.frame(table(out$Pathway))

muestras_up <- comparacion_long[which(comparacion_long$Color == "UP"),]

# voy a pintar aquellos puntos que resultaron out para >10% de las muestras
out_freq_species <- out_freq$Var1[which(out_freq$Freq>=10)]
out_freq_species <- out_freq$Var1

# Paso 3: Graficar Bland-Altman con etiquetas para puntos fuera de los límites
library(ggplot2)
library(ggrepel)
library(dplyr)

length(unique(comparacion_long$Pathway))
#out_freq_species <- out_freq_species[-4]

ggplot(comparacion_long, aes(x = Mean, y = Difference)) +
  geom_point(aes(color = ifelse(Pathway %in% out_freq_species & Color != "Dentro de los límites", Pathway, NA)), size = 2.5) +
  #geom_point(aes(color = ifelse( Color != "Dentro de los límites", Species, NA)), size = 2) +
  geom_hline(yintercept = Diff_Mean, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = Lower_Limit, linetype = "dotted", color = "red") +
  geom_hline(yintercept = Upper_Limit, linetype = "dotted", color = "red") +

  labs(title = "Bland-Altman Plot - CPM 239 Pathways",
       x = "Mean (Bowtie & BWA)",
       y = "Difference (Bowtie - BWA)",
       color = "Pathway") +
  #scale_y_log10() +  # Escala logarítmica para el eje y
  #scale_y_continuous(trans = 'sqrt') +  # Aplica una transformación cuadrática
  theme_minimal() +
  # Ajustar tamaño de los textos de los ejes y la leyenda
  theme(
    axis.title = element_text(size = 16),   # Tamaño de títulos de los ejes
    axis.text = element_text(size = 14),    # Tamaño de los valores de los ejes
    legend.title = element_text(size = 10), # Tamaño del título de la leyenda
    legend.text = element_text(size = 8),   # Tamaño de los valores de la leyenda
    plot.title = element_text(size = 18)
  )





############################################################
# BUSCO DIF SIGN P VALORES CON LAS DIF TABLAS DE AR DE VIAS DE DEHOST

generatePathwaysTable(patients_dir = "~/Daniela/Biota/Muestras/73m", de_host = "Bowtie")
generatePathwaysTable(patients_dir = "~/Daniela/Biota/Muestras/73m", de_host = "RSubread")
generatePathwaysTable(patients_dir = "~/Daniela/Biota/Muestras/73m", de_host = "BWA")
generatePathwaysTable(patients_dir = "~/Daniela/Biota/Muestras/73m", de_host = "")

AR_vias_Bo <- as.data.frame(read_excel("~/Daniela/Biota/Muestras/SubsetPathways/Bo/Bo_AR_Vias_83p.xlsx"))
AR_vias_T <- as.data.frame(read_excel("~/Daniela/Biota/Muestras/SubsetPathways/T/T_AR_Vias_83p.xlsx"))
AR_vias_BWA <- as.data.frame(read_excel("~/Daniela/Biota/Muestras/SubsetPathways/bwa/bwa_AR_Vias_83p.xlsx"))
AR_vias_Rs <- as.data.frame(read_excel("~/Daniela/Biota/Muestras/SubsetPathways/Rs/Rs_AR_Vias_83p.xlsx"))

vias_comun <- intersect(AR_vias_Bo$Pathway, AR_vias_BWA$Pathway)
AR_vias_Bo <- AR_vias_Bo[which(AR_vias_Bo$Pathway %in% vias_comun),]
AR_vias_BWA <- AR_vias_BWA[which(AR_vias_BWA$Pathway %in% vias_comun),]

AR_vias_Bo[,-c(1,2,3)] <- prop.table(as.matrix(AR_vias_Bo[,-c(1,2,3)]), margin = 2) * 100
colSums(AR_vias_Bo[,-c(1,2,3)])
AR_vias_Bo[, "36"] <- 0

AR_vias_BWA[,-c(1,2,3)] <- prop.table(as.matrix(AR_vias_BWA[,-c(1,2,3)]), margin = 2) * 100
colSums(AR_vias_BWA[,-c(1,2,3)])
AR_vias_BWA[, "36"] <- 0


p_por_categoriaBo <- generar_p_por_categoria(CPM_vias = AR_vias_Bo, de_host_file = "Bo")
p_por_categoriaBWA <- generar_p_por_categoria(CPM_vias = AR_vias_BWA, de_host_file = "bwa")
p_por_categoriaRs <- generar_p_por_categoria(CPM_vias = AR_vias_Rs, de_host_file = "Rs")
p_por_categoriaT <- generar_p_por_categoria(CPM_vias = AR_vias_T, de_host_file = "T")

p_por_categoriaBo <- as.data.frame(read_excel("~/Daniela/Biota/Muestras/SubsetPathways/Bo/Bo_p_por_categoria-274vias_AR_83p.xlsx"))
p_signBo <- p_por_categoriaBo[which(p_por_categoriaBo$Wilcoxon<0.05),]

p_por_categoriaBWA <- as.data.frame(read_excel("~/Daniela/Biota/Muestras/SubsetPathways/bwa/bwa_p_por_categoria-274vias_AR_83p.xlsx"))
p_signBWA <- p_por_categoriaBWA[which(p_por_categoriaBWA$Wilcoxon<0.05),]

p_por_categoriaRs <- as.data.frame(read_excel("~/Daniela/Biota/Muestras/SubsetPathways/Rs/Rs_p_por_categoria-274vias_AR_83p.xlsx"))
p_signRs <- p_por_categoriaRs[which(p_por_categoriaRs$Wilcoxon<0.05),]

p_por_categoriaT <- as.data.frame(read_excel("~/Daniela/Biota/Muestras/SubsetPathways/T/T_p_por_categoria-274vias_AR_83p.xlsx"))
p_signT <- p_por_categoriaT[which(p_por_categoriaT$Wilcoxon<0.05),]


# Filas en Bowtie que no están en BWA
colnames(p_signBo)
loquepierdeBWA <- anti_join(p_signBo, p_signBWA, by = c("Categoria" = "Categoria", "Via" = "Via", "GrupoDominante" = "GrupoDominante"))
colnames(loquepierdeBWA)[3] <- "P_Bo"
loquepierdeBWA <- merge(loquepierdeBWA, p_por_categoriaBWA, by = c("Categoria" = "Categoria", "Via" = "Via", "GrupoDominante" = "GrupoDominante"), all = F)
colnames(loquepierdeBWA)[5] <- "P_BWA"

unique(loquepierdeBWA$Categoria)
loquepierdeBWA_imp <- loquepierdeBWA[which(loquepierdeBWA$Categoria %in% c("Sexo", "Rango etario")),]
all(loquepierdeBWA_imp$Via %in% AR_vias_BWA$Pathway)

loquepierdeBWA_imp$DifP <- abs(loquepierdeBWA_imp$P_Bo - loquepierdeBWA_imp$P_BWA)


##############################################################################################

de_host_file = "bwa"
de_host_file = "Bo"
ncol(MetadataB)

CPM_vias_completo_Bo <- as.data.frame(read_excel("~/Daniela/Biota/Muestras/SubsetPathways/Bo/Bo_vias_completo_AR_83p.xlsx"))
CPM_vias_completo_bwa <- as.data.frame(read_excel("~/Daniela/Biota/Muestras/SubsetPathways/bwa/bwa_vias_completo_AR_83p.xlsx"))

MetadataB <- as.data.frame(read_excel("~/Daniela/Biota/PipelineBiota/data/Metadata-soloColumnasUsables.xlsx"))

CPM_vias_completo_Bo <- as.data.frame(read_excel("~/Daniela/Biota/Muestras/SubsetPathways/Bo/Bo_239vias_completo_AR_83p.xlsx"))
CPM_vias_completo_bwa <- as.data.frame(read_excel("~/Daniela/Biota/Muestras/SubsetPathways/bwa/bwa_239vias_completo_AR_83p.xlsx"))


df_comparacion <- merge(CPM_vias_completo_Bo, CPM_vias_completo_bwa, by = colnames(MetadataB), all = TRUE)
colnames(df_comparacion) <- gsub("\\.x$", "_Bo", colnames(df_comparacion))
colnames(df_comparacion) <- gsub("\\.y$", "_bwa", colnames(df_comparacion))
library(tidyr)
library(ggplot2)

# Seleccionar las columnas necesarias: Sexo y las vías de interés
vias_interes <- colnames(df_comparacion)[which(grepl("_Bo|_bwa", colnames(df_comparacion)))]
df_vias <- df_comparacion[, c("ID", "Sexo", "Rango etario", vias_interes)]

# Reorganizar a formato largo
df_vias_long <- pivot_longer(
  df_vias,
  cols = vias_interes, # Columnas de las vías
  names_to = c("Via", "Metodo"),
  names_sep = "_", # Separar por "_"
  values_to = "Valor"
)

# Definir las vías de interés
p$Via[c(4,5,9,10,11,12,13,20)]

vias_interes <- c("P41-PWY", "PANTO-PWY","PWY-6385", "PWY-6386", "PWY-702","PWY-7220",  "PWY-7222",  "THRESYN-PWY")
vias_interes <- c("RIBOSYN2-PWY", "ARGININE-SYN4-PWY", "PWY-7198", "SER-GLYSYN-PWY") # Añadir las vías que desees graficar

# Filtrar el dataframe para incluir únicamente las vías seleccionadas
df_vias_filtradas <- df_vias_long %>%
  filter(Via %in% vias_interes)

# Crear el boxplot
df_vias_filtradas$`Rango etario`
df_vias_filtradas$Sexo

ggplot(df_vias_filtradas, aes(#x = Sexo,
                              x = `Rango etario`,
                              y = Valor, fill = Metodo)) +
  geom_boxplot() +
  facet_wrap(~ Via, scales = "free_y") + # Una faceta por vía
  labs(
    title = "Comparación de varias vías por dehost y Rango Etario",
    x = "Rango Etario",
    y = "Valor",
    fill = "Método"
  ) +
  theme_minimal()

p_vias_interes <- p_significativos[p_significativos$Via %in% vias_interes,]
p_vias_interes[,c(4:6)] <- round(p_vias_interes[,c(4:6)], 4)


p_significativos <- loquepierdeBWA_imp
de_host_file <- "Bo"

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


##############################################################################################



generar_p_por_categoria <- function(CPM_vias, de_host_file) {

  CPM_viasT <- as.data.frame(t(CPM_vias))
  colnames(CPM_viasT) <- CPM_viasT[1,]
  CPM_viasT <- CPM_viasT[-c(1,2,3),]

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

  write.xlsx(CPM_vias_completo, file =  sprintf("~/Daniela/Biota/Muestras/SubsetPathways/%s/%s_%svias_completo_AR_83p.xlsx", de_host_file, de_host_file, length(2:n)))


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

      p_por_categoria[i, "GrupoDominante"] <-  ifelse(length(grupo_dominante) == 0, NA, as.character(grupo_dominante))

      i <- i+1
    }
  }

  write.xlsx(p_por_categoria, file =  sprintf("~/Daniela/Biota/Muestras/SubsetPathways/%s/%s_p_por_categoria-274vias_AR_83p.xlsx", de_host_file, de_host_file))
  return(p_por_categoria)

}
