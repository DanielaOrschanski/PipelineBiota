# Cómo son las AR de los generos mas importantes en cada metodología? (segun paper) ---------------------------------------------
patients_dir = "~/Daniela/Biota/Muestras/73m"
source = "KRAKEN"
source = "DRAGEN"
de_host = "Bowtie"

extraerAR_Generos <- function(patients_dir, source, de_host ) {
  tabla_otus <- generateOTUsTableGrupal(patients_dir = patients_dir, source = source, conEukaryota = FALSE, de_host = de_host)

  if(any(tabla_otus$Species == "Cutibacterium modestum")) {
    tabla_otus$Species[which(tabla_otus$Species == "Cutibacterium modestum")] <- "Propionibacterium sp. oral taxon 193"
    tabla_otus$Genus[which(tabla_otus$Species ==  "Propionibacterium sp. oral taxon 193")] <- "Propionibacterium"
  }

  out <- group_TaxonomicLevels(patients_dir = patients_dir, tabla_otus = tabla_otus, source = source, de_host = de_host, conE = FALSE)
  list_AR <- out[[1]]
  #AR_filos <- list_AR[[3]]
  AR_generos <- list_AR[[7]]
  #AR_species <- list_AR[[8]]
  colnames(AR_generos)[-1] <- gsub(sprintf("_%s", source), "", colnames(AR_generos)[-1])

  #MetadataB <- as.data.frame(read_excel("~/Daniela/Biota/PipelineBiota/data/Metadata-soloColumnasUsables.xlsx"))
  #femeninas <- MetadataB$ID[which(MetadataB$Sexo == "Masculino")]
  #AR_species_fem <- AR_species[, which(colnames(AR_species) %in% c("Species", femeninas))]

  #AR_generos$Promedio <- rowMeans(AR_generos[,-1])
  #otus_unique <- tabla_otus[!duplicated(tabla_otus$Species), c("Species", "Phylum")]
  #AR_species_fem <- merge(AR_species_fem, otus_unique, by = "Species", all.x = TRUE, all.y = FALSE)
  #AR_generos <- AR_generos[order(AR_generos$Promedio, decreasing = TRUE),]
  #AR_generos <- AR_species_fem[-which(AR_species_fem$Phylum == "-"),]

  return(AR_generos)
}

extraerAR_Filos <- function(patients_dir, source, de_host ) {
  tabla_otus <- generateOTUsTableGrupal(patients_dir = patients_dir, source = source, conEukaryota = FALSE, de_host = de_host)
  if(any(tabla_otus$Species == "Cutibacterium modestum")) {
    tabla_otus$Species[which(tabla_otus$Species == "Cutibacterium modestum")] <- "Propionibacterium sp. oral taxon 193"
    tabla_otus$Genus[which(tabla_otus$Species ==  "Propionibacterium sp. oral taxon 193")] <- "Propionibacterium"
  }
  out <- group_TaxonomicLevels(patients_dir = patients_dir, tabla_otus = tabla_otus, source = source, de_host = de_host, conE = FALSE)
  list_AR <- out[[1]]
  AR_filos <- list_AR[[3]]
  #AR_generos <- list_AR[[7]]
  #AR_species <- list_AR[[8]]
  colnames(AR_filos)[-1] <- gsub(sprintf("_%s", source), "", colnames(AR_filos)[-1])

  return(AR_filos)
}

extraerAR_Familias <- function(patients_dir, source, de_host ) {
  tabla_otus <- generateOTUsTableGrupal(patients_dir = patients_dir, source = source, conEukaryota = FALSE, de_host = de_host)
  tabla_otus$Species[which(tabla_otus$Species == "Cutibacterium modestum")] <- "Propionibacterium sp. oral taxon 193"
  tabla_otus$Genus[which(tabla_otus$Species ==  "Propionibacterium sp. oral taxon 193")] <- "Propionibacterium"

  out <- group_TaxonomicLevels(patients_dir = patients_dir, tabla_otus = tabla_otus, source = source, de_host = de_host, conE = FALSE)
  list_AR <- out[[1]]
  #AR_filos <- list_AR[[3]]
  AR_familias <- list_AR[[6]]
  colnames(AR_familias)[-1] <- gsub(sprintf("_%s", source), "", colnames(AR_familias)[-1])
  return(AR_familias)
}

library(PipelineBiota)

path.expand("~/Daniela/Biota")
#K:
AR_genK_Bo <- extraerAR_Generos(patients_dir = "~/Daniela/Biota/Muestras/73m", source = "KRAKEN", de_host = "Bowtie" )
AR_genK_BWA <- extraerAR_Generos(patients_dir = "~/Daniela/Biota/Muestras/73m", source = "KRAKEN", de_host = "BWA" )
AR_genK_Rs <- extraerAR_Generos(patients_dir = "~/Daniela/Biota/Muestras/73m", source = "KRAKEN", de_host = "RSubread" )
AR_genK_sin <- extraerAR_Generos(patients_dir = "~/Daniela/Biota/Muestras/73m", source = "KRAKEN", de_host = "" )

#D:
AR_genD_Bo <- extraerAR_Generos(patients_dir = "~/Daniela/Biota/Muestras/73m", source = "DRAGEN", de_host = "Bowtie" )
AR_genD_BWA <- extraerAR_Generos(patients_dir = "~/Daniela/Biota/Muestras/73m", source = "DRAGEN", de_host = "BWA" )
AR_genD_Rs <- extraerAR_Generos(patients_dir = "~/Daniela/Biota/Muestras/73m", source = "DRAGEN", de_host = "RSubread" )
AR_genD_sin <- extraerAR_Generos(patients_dir = "~/Daniela/Biota/Muestras/73m", source = "DRAGEN", de_host = "" )
ARgenD_sinDH_PD <-  extraerAR_Generos(patients_dir = "~/Daniela/Biota/Muestras/73m", source = "DRAGEN", de_host = "sinDH_PD")

lista_df <- list("bwaK" = AR_genK_BWA, "BoK" = AR_genK_Bo, "RsK" = AR_genK_Rs, "-K" = AR_genK_sin,
                 "bwaD" = AR_genD_BWA, "BoD" = AR_genD_Bo, "RsD" = AR_genD_Rs, "DdhD" = AR_genD_sin,
                 "-D" = ARgenD_sinDH_PD)

colSums(AR_genK_BWA[,-1])
colSums(ARgenD_sinDH_PD[,-1])

# Unir dataframes en uno solo, agregando una columna con el nombre de la metodología
df_combined <- bind_rows(lapply(names(lista_df), function(nombre) {
  lista_df[[nombre]] %>%
    mutate(Metodologia = nombre)
}), .id = "ID")
df_combined <- df_combined[,-1]


# Filtrar solo los géneros de interés
generos_imp <- c("Staphylococcus", "Streptococcus", "Acinetobacter", "Cutibacterium", "Corynebacterium", "Micrococcus")
#generos_imp <- c("Staphylococcus", "Streptococcus", "Acinetobacter", "Cutibacterium", "Corynebacterium", "Bacteroides", "Bacillus", "Lactobacillus", "Pseudomonas", "Citrobacter")

df_filtered <- df_combined %>%
  filter(Genus %in% generos_imp)
length(unique(df_filtered$Genus))

#df_filtered <- df_combined
df_filtered$Source <- ifelse(grepl("K", df_filtered$Metodologia), "Kraken", "DRAGEN")
#df_filtered$Genus <- factor(df_filtered$Genus, levels = rev(generos_imp))
df_filtered$Metodologia <- factor(df_filtered$Metodologia, levels = names(lista_df))

#Agregar promedio fem y promedio masc:
MetadataB <- as.data.frame(read_excel("~/Daniela/Biota/PipelineBiota/data/Metadata-soloColumnasUsables.xlsx"))
#femeninas <- MetadataB$ID[MetadataB$Sexo == "Femenino"]
#masculinas <- MetadataB$ID[MetadataB$Sexo == "Masculino"]

#df_fem <- df_filtered[, colnames(df_filtered) %in% femeninas]
#df_mas <- df_filtered[, colnames(df_filtered) %in% masculinas]

# Calcular el promedio por fila para cada grupo
#df_filtered$PromedioFem <- rowMeans(df_fem, na.rm = TRUE)
#df_filtered$PromedioMas <- rowMeans(df_mas, na.rm = TRUE)


library(tidyr)

# Reshape los datos para tener una columna de Sexo (PromedioFem o PromedioMas)
#df_long <- df_filtered %>%
#  pivot_longer(cols = c("PromedioFem", "PromedioMas"),
#               names_to = "Sexo",
#               values_to = "Prom")

# Ahora crea el gráfico de barras con posición dodge para separar las barras de cada sexo por metodología
ggplot(df_long, aes(x = Metodologia, y = Prom, fill = Genus, group = interaction(Metodologia, Sexo))) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  labs(title = "Distribución de géneros por metodología y sexo",
       x = "Metodología",
       y = "Abundancia relativa promedio",
       fill = "Género") +
  scale_fill_brewer(palette = "Set3") +  # Paleta de colores para los géneros
  scale_y_continuous(breaks = seq(0, 110, 10)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotar etiquetas del eje X
  facet_grid(~ Source, scales = "free_x", space = "free_x")  # Facetas por Source



# Gráfico de barras apiladas como paper:
ggplot(df_filtered, aes(x = Promedio, y = Metodologia, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  labs(title = "Comparación géneros",
       x = "Abundancia relativa promedio",
       y = "Metodología",
       fill = "Género") +
  #scale_fill_brewer(palette = "Set2") +
  scale_fill_manual(values = c(
    "Staphylococcus" = "#9370DB",  # Violeta pastel
    "Streptococcus"  = "#6959CD",  # Lila pastel
    "Acinetobacter"  = "#A6BDDB",  # Azul pastel
    "Cutibacterium"  = "#C6DBEF",  # Celeste pastel
    "Corynebacterium"= "#FEE391",  # Amarillo pastel
    "Micrococcus"    = "#F4AEB5"   # Rosa pastel
  )) +
  scale_x_continuous(breaks = seq(0, 110,  10)) +
  theme_minimal() +
  theme(strip.text.y = element_blank()) +  # Oculta los nombres del facet
  facet_grid(Source ~ ., scales = "free_y", space = "free_y")
#+  facet_grid(~ Source, scales = "free_x", space = "free_x")

df_comparacion <- df_filtered[, c("Genus", "Promedio", "Metodologia", "Source")]
df_comparacion[df_comparacion$Genus == "Staphylococcus",]
df_comparacion[df_comparacion$Genus == "Cutibacterium",]
df_comparacion[df_comparacion$Genus == "Streptococcus",]
df_comparacion[df_comparacion$Genus == "Acinetobacter",]

# Calcular el coeficiente de Pearson:
df_wide <- df_comparacion %>%
  select(Genus, Promedio, Metodologia) %>%
  pivot_wider(names_from = Metodologia, values_from = Promedio)

cor_test_result <- cor.test(df_wide$Ksin, df_wide$KBWA, method = "pearson")
cor_test_result$p.value #la correlacion es significativa
cor_test_result$estimate # cercano a 1 es alta correlacion

cor_test_result <- cor.test(df_wide$Ksin, df_wide$KBo, method = "pearson")
cor_test_result <- cor.test(df_wide$KBo, df_wide$DBo, method = "pearson")
cor_test_result <- cor.test(df_wide$KBWA, df_wide$DBWA, method = "pearson")

###################################################################
# Bland altman KBo vs KBWA: --------------------------------------

# Unir por genero: se queda con los generos en comun
df_wide <- inner_join(AR_genK_Bo, AR_genK_BWA, by = "Genus", suffix = c("_KBo", "_KBWA"))
generos_comun <- intersect(AR_genK_Bo$Genus, AR_genK_BWA$Genus)
all(generos_comun == df_wide$Genus)

# Transformar de formato ancho a largo para facilitar el cálculo
df_long <- df_wide %>%
  pivot_longer(cols = -Genus,
               names_to = c("Muestra", "Metodologia"),
               names_pattern = "(.*)_(KBo|KBWA)",
               values_to = "Valor") %>%
  pivot_wider(names_from = Metodologia, values_from = Valor)

# Calcular el promedio y la diferencia
df_long <- df_long %>%
  mutate(Promedio = (KBo + KBWA) / 2,
         Diferencia = KBo - KBWA)

# Graficar Bland-Altman
ggplot(df_long, aes(x = Promedio, y = Diferencia)) +
  geom_point(alpha = 0.5) +  # Puntos con transparencia
  geom_hline(yintercept = mean(df_long$Diferencia), linetype = "dashed", color = "red") +  # Media de la diferencia
  geom_hline(yintercept = mean(df_long$Diferencia) + 1.96 * sd(df_long$Diferencia), linetype = "dotted", color = "blue") +  # Límites de acuerdo
  geom_hline(yintercept = mean(df_long$Diferencia) - 1.96 * sd(df_long$Diferencia), linetype = "dotted", color = "blue") +
  labs(title = "Bland-Altman Plot",
       x = "Promedio (KBo + KBWA) / 2",
       y = "Diferencia (KBo - KBWA)") +
  theme_minimal()

# FIGURE 4 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# Por genero: ---------------------
generos_imp <- c("Staphylococcus", "Streptococcus", "Acinetobacter", "Cutibacterium", "Corynebacterium", "Micrococcus")

library(tidyverse)
df_KBo <- AR_genK_Bo %>% filter(Genus %in% generos_imp)
df_KBo$Method <- "BoK"
df_KBWA <- AR_genK_BWA %>% filter(Genus %in% generos_imp)
df_KBWA$Method <- "bwaK"
df_KRs <- AR_genK_Rs %>% filter(Genus %in% generos_imp)
df_KRs$Method <- "RsK"
df_K <- AR_genK_sin %>% filter(Genus %in% generos_imp)
df_K$Method <- "-K"

df_DBo <- AR_genD_Bo %>% filter(Genus %in% generos_imp)
df_DBo$Method <- "BoD"
df_DBWA <- AR_genD_BWA %>% filter(Genus %in% generos_imp)
df_DBWA$Method <- "bwaD"
df_DRs <- AR_genD_Rs %>% filter(Genus %in% generos_imp)
df_DRs$Method <- "RsD"
df_DdhD <- AR_genD_sin %>% filter(Genus %in% generos_imp)
df_DdhD$Method <- "DdhD"
df_D <- ARgenD_sinDH_PD %>% filter(Genus %in% generos_imp)
df_D$Method <- "-D"

#Recalcular AR solo con los generos:

pasarAR <- function(df) {
  str(df)
  AR <- prop.table(as.matrix(df[,-c(1,85,86)]), margin = 2)*100
  colSums(AR)
  df <- cbind(df[,c(1,86)], AR)
  return(df)
}

#df_wide <- inner_join(df_KBo, df_KBWA, by = "Genus", suffix = c("_KBo", "_KBWA"))
dfs <- list(df_KBo, df_KBWA, df_KRs, df_K, df_DBo, df_DBWA, df_DRs, df_DdhD, df_D)
#dfs <- lapply(dfs, pasarAR)

df <- dfs[[1]]
summary(colSums(df[,-c(1,85)]))
str(df)

suffixes <- map_chr(dfs, ~ unique(.x$Method))

# Función para limpiar y renombrar columnas antes de unir
limpiar_df <- function(df, suffix) {
  df %>%
    #select(-Promedio, -Method) %>%  # Eliminamos 'promedio' y 'Method'
    select( -Method) %>%
    rename_with(~ paste0(., "_", suffix), -Genus) # Agregamos sufijo a todas las columnas excepto 'Genus'
}

dfs_limpios <- map2(dfs, suffixes, limpiar_df)

df_wide <- reduce(dfs_limpios, inner_join, by = "Genus")

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(purrr)

# Definir pares de métodos a comparar
#Fig4a:
comparisons <- list(
  c("BoK", "BoD"),
  c("bwaK", "bwaD"),
  c("RsK", "RsD"),
  c("-K", "-D")
)

#Fig4b:
comparisons <- list(
  c("BoK", "-K"),
  c("BoK", "RsK"),
  c("BoK", "bwaK")
)

#Fig4c:
comparisons <- list(
  c("DdhD", "-D"),
  c("DdhD", "BoD"),
  c("DdhD", "RsD"),
  c("DdhD", "bwaD")
)


df_long <- df_wide %>%
  pivot_longer(cols = -Genus,
               names_to = c("Sample", "Methodology"),
               names_sep = "_",
               values_to = "Value")

unique(df_long$Methodology)
# Generar los datos para Bland-Altman con todas las comparaciones
df_bland <- bind_rows(lapply(comparisons, function(par) {
  df_long %>%
    filter(Methodology %in% par) %>%
    pivot_wider(names_from = Methodology, values_from = Value) %>%
    mutate(
      Comparison = paste(par[1], "vs", par[2]),
      Mean = (as.numeric(.data[[par[1]]]) + as.numeric(.data[[par[2]]])) / 2,
      Difference = as.numeric(.data[[par[1]]]) - as.numeric(.data[[par[2]]]),
      DifferencePerc = (as.numeric(.data[[par[1]]]) - as.numeric(.data[[par[2]]])) / as.numeric(.data[[par[1]]]) * 100
    )
}))

# Calcular límites de Bland-Altman para cada comparación
df_lims <- df_bland %>%
  group_by(Genus, Comparison) %>%
  summarise(
    Diff_mean = mean(Difference, na.rm = TRUE),
    Diff_SD = sd(Difference, na.rm = TRUE),
    Upper_limit = Diff_mean + 1.96 * Diff_SD,
    Lower_limit = Diff_mean - 1.96 * Diff_SD
  )

# Unir los límites con los datos originales
df_final <- df_bland %>%
  left_join(df_lims, by =c("Genus", "Comparison"))

df_final <- df_final %>%
  mutate(Comparison = as.factor(Comparison))


# Calcular límites globales del eje Y
#y_min <- min(df_final$Diferencia, na.rm = TRUE)
#y_max <- max(df_final$Diferencia, na.rm = TRUE)
#x_min <- min(df_final$Promedio, na.rm = TRUE)
#x_max <- max(df_final$Promedio, na.rm = TRUE)

library(ggplot2)
library(gridExtra)
library(dplyr)
library(purrr)

# Crear los gráficos con mejor estética
plots <- df_final %>%
  split(.$Genus) %>%
  map(~ ggplot(.x, aes(x = Mean, y = Difference, color = Comparison, shape = Comparison)) +
        geom_point(size = 1.5, alpha = 0.8) +
        geom_hline(aes(yintercept = Diff_mean, color = Comparison), linetype = "dashed", linewidth = 0.4) +
        geom_hline(aes(yintercept = Upper_limit, color = Comparison), linetype = "dotted", linewidth = 0.4) +
        geom_hline(aes(yintercept = Lower_limit, color = Comparison), linetype = "dotted", linewidth = 0.4) +
        labs(title = unique(.x$Genus)) +
        theme_classic() +
        theme(
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none",
          plot.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
        )
      +
        scale_color_manual(values = c("black", "blue", "red", "green")) +
        scale_shape_manual(values = c(16, 17, 15, 18))
      #+
      #  scale_color_manual(values = c( "blue", "red", "black")) +
      #  scale_shape_manual(values = c(17, 15, 16))
  )

# Extraer la leyenda de un gráfico
legend_plot <- ggplot(df_final, aes(x = Mean, y = Difference, color = Comparison, shape = Comparison)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("black", "blue", "red", "green")) +
  scale_shape_manual(values = c(16, 17, 15, 18)) +
  #scale_color_manual(values = c("blue", "red", "black")) +
  #scale_shape_manual(values = c(17, 15, 16)) +
  theme_minimal() +
  theme(legend.position = "right")

legend <- cowplot::get_legend(legend_plot)


plots <- lapply(plots, function(p) {
  p + scale_y_continuous(labels = scales::number_format(accuracy = 0.01))  # Ajusta los decimales en el eje y
})


# Ahora los gráficos deberían alinearse correctamente
grid.arrange(
  grobs = plots,
  ncol = 2,
  bottom = textGrob("Average Relative Abundance", gp = gpar(fontsize = 12)),
  #left = textGrob("Difference", gp = gpar(fontsize = 12, fontface = "bold"), rot = 90),
  left = textGrob("Difference: Kraken - DRAGEN", gp = gpar(fontsize = 12), rot = 90),
  #left = textGrob("Difference: BoK vs Other Kraken Methods", gp = gpar(fontsize = 12, fontface = "bold"), rot = 90),
  #left = textGrob("Difference: DdhD vs Other DRAGEN Methods", gp = gpar(fontsize = 12, fontface = "bold"), rot = 90),

  #top = textGrob("Kraken vs DRAGEN Comparison", gp = gpar(fontsize = 16, fontface = "bold")),
  #top = textGrob("Bowtie+Kraken(BoK) vs Other Kraken Methods", gp = gpar(fontsize = 16, fontface = "bold")),
  #top = textGrob("de-hostDRAGEN+DRAGEN(DdhD) vs Other DRAGEN Methods", gp = gpar(fontsize = 16, fontface = "bold")),
  right = legend
)

unique(df_final$Comparison)
unique(df_final$Genus)
library(dplyr)


#Overestimated by kraken:
df_final %>%
  filter(Genus == "Cutibacterium", DifferencePerc > 0) %>%
  group_by(Comparison) %>%
  summarise(n_samples_over = n_distinct(Sample))
df_final %>%
  filter(Genus == "Cutibacterium") %>%
  group_by(Comparison) %>%
  summarise(
    n_samples = n_distinct(Sample),
    avg_diff_perc = mean(DifferencePerc, na.rm = TRUE),
    .groups = "drop"
  )

df_final %>%
  filter(Genus == "Staphylococcus", DifferencePerc > 0) %>%
  group_by(Comparison) %>%
  summarise(n_samples_over = n_distinct(Sample))
df_final %>%
  filter(Genus == "Staphylococcus") %>%
  group_by(Comparison) %>%
  summarise(
    n_samples = n_distinct(Sample),
    avg_diff_perc = mean(DifferencePerc, na.rm = TRUE),
    .groups = "drop"
  )

#Underestimated by kraken:
df_final %>%
  filter(Genus == "Acinetobacter", DifferencePerc < 0) %>%
  group_by(Comparison) %>%
  summarise(n_samples_over = n_distinct(Sample))
df_final %>%
  filter(Genus == "Acinetobacter") %>%
  group_by(Comparison) %>%
  summarise(
    n_samples = n_distinct(Sample),
    avg_diff_perc = mean(DifferencePerc, na.rm = TRUE),
    .groups = "drop"
  )

# Función para calcular p-valor de Wilcoxon
calcular_pvalor <- function(df, metodo1, metodo2) {
  df %>%
    group_by(Genus) %>%
    summarise(
      p_value = ifelse(
        sum(!is.na(get(metodo1))) > 1 & sum(!is.na(get(metodo2))) > 1, # Asegura que haya más de un valor
        wilcox.test(get(metodo1), get(metodo2), paired = FALSE)$p.value,
        NA
      )
    ) %>%
    mutate(Comparacion = paste(metodo1, "vs", metodo2))
}

# Aplicar la función a todas las comparaciones
resultados <- lapply(comparaciones, function(par) calcular_pvalor(df_final, par[1], par[2]))

# Combinar los resultados en un solo dataframe
df_pvalores <- bind_rows(resultados)
df_pvalores

# Calcular MSE por cada comparación dentro de cada Genus
mse_results <- df_final %>%
  group_by(Genus, Comparacion) %>%
  summarise(
    MSE = mean(Diferencia^2, na.rm = TRUE),  # MSE = promedio de los errores al cuadrado
    RMSE = sqrt(MSE),  # Raíz del MSE para interpretarlo en las mismas unidades
    .groups = "drop"
  )

# Mostrar resultados
print(mse_results)

# Mostrar resultados
###################################

plots <- df_final %>%
  split(.$Genus) %>%
  map(~ ggplot(.x, aes(x = Promedio, y = Diferencia)) +
        geom_point(alpha = 0.5) +
        geom_hline(aes(yintercept = Diff_mean), linetype = "dashed", color = "red") +
        geom_hline(aes(yintercept = Upper_limit), linetype = "dotted", color = "blue") +
        geom_hline(aes(yintercept = Lower_limit), linetype = "dotted", color = "blue") +
        # Agregar etiquetas SOLO a los puntos fuera de los límites
        geom_text_repel(aes(label = ifelse(Diferencia > Upper_limit | Diferencia < Lower_limit, Muestra, "")),
                        color = "black", size = 3, max.overlaps = 15) +

        labs(title = paste("", unique(.x$Genus)),
             #y = "Diferencia (KBo - KBWA)") +
             x = "Promedio",
             #y = "Diferencia (KBo - DBo)") +
             y = "Diferencia (KBo - Ksin)") +
        theme_minimal())

library(gridExtra)
grid.arrange(grobs = plots, ncol = 3, nrow = 2)


library(ggplot2)
library(dplyr)

# Ejemplo de cómo estructurar las diferencias y medias
df1 <- df_KBo
df2 <- df_DBo
df3 <- df_KBWA
df4 <- df_DBWA
df5 <- df_KRs
df6 <- df_DRs
df

bland_altman <- function(df1, df2) {

  #suf1 <- strsplit(deparse(substitute(df_1)), split ="_")[[1]][2]
  #suf2 <- strsplit(deparse(substitute(df_2)), split ="_")[[1]][2]
  suf1 <- unique(df1$Method)
  suf2 <- unique(df2$Method)
  suf3 <- unique(df3$Method)
  suf4 <- unique(df4$Method)

  #Para hacer solo entre 2:
  #df_wide <- inner_join(df1[,!colnames(df1) %in% "Method"],
  #                      df2[,!colnames(df2) %in% "Method"], by = "Genus", suffix = c(suf1, suf2))

  # Lista de dataframes sin la columna "Method"
  df_list <- list(df1, df2, df3, df4) %>%
    map(~ select(.x, -Method))

  # Unir iterativamente usando reduce() y inner_join()
  df_wide <- reduce(df_list, inner_join, by = "Genus", suffix = c(suf1, suf2, suf3, suf4))

  df_long <- df_wide %>%
    pivot_longer(cols = -Genus,
                 names_to = c("Muestra", "Metodologia"),
                 names_pattern = sprintf("(.*)(%s|%s)", suf1, suf2),
                 values_to = "Valor") %>%
    pivot_wider(names_from = Metodologia, values_from = Valor)

  # Calcular el promedio y la diferencia
  str(df_long)
  df_long <- as.data.frame(df_long)
  df_long <- df_long %>%
    mutate(Promedio = (as.numeric(df_long[,3]) + as.numeric(df_long[,4])) / 2,
           Diferencia = as.numeric(df_long[,3]) - as.numeric(df_long[,4]))
  #df_long <- as.data.frame(df_long)
  #colnames(df_long)[c(5,6)] <- c("Promedio", "Diferencia")
  #str(df_long)

  # Calcular los límites para cada género
  df_lims <- df_long %>%
    group_by(Genus) %>%
    summarise(Diff_mean = mean(Diferencia, na.rm = TRUE),
              Diff_SD = sd(Diferencia, na.rm = TRUE)) %>%
    mutate(Upper_limit = Diff_mean + 1.96 * Diff_SD,
           Lower_limit = Diff_mean - 1.96 * Diff_SD)

  # Unir límites con los datos originales
  df_final <- df_long %>%
    left_join(df_lims, by = "Genus")

  # Crear gráficos independientes
  library(ggrepel)

  plots <- df_final %>%
    split(.$Genus) %>%
    map(~ ggplot(.x, aes(x = Promedio, y = Diferencia)) +
          geom_point(alpha = 0.5) +
          geom_hline(aes(yintercept = Diff_mean), linetype = "dashed", color = "red") +
          geom_hline(aes(yintercept = Upper_limit), linetype = "dotted", color = "blue") +
          geom_hline(aes(yintercept = Lower_limit), linetype = "dotted", color = "blue") +
          # Agregar etiquetas SOLO a los puntos fuera de los límites
          geom_text_repel(aes(label = ifelse(Diferencia > Upper_limit | Diferencia < Lower_limit, Muestra, "")),
                          color = "black", size = 3, max.overlaps = 15) +
          labs(title = paste("", unique(.x$Genus)),
               x = "Promedio",
               y = sprintf("%s - %s", suf1, suf2)) +
          theme_minimal())

  library(gridExtra)
  grid.arrange(grobs = plots, ncol = 3, nrow = 2)

  return(plots)
}

plots[1]
# Generar Bland-Altman para cada comparación
bland_altman_DBo_KBo <- bland_altman(df1 = df_DBo, df2 = df_KBo)
bland_altman_DBWA_KBWA <- bland_altman(df1 = df_DBWA, df2 = df_KBWA)
bland_altman_DRs_KRs <- bland_altman(df1 = df_DRs, df2 = df_KRs)
bland_altman_D_K <- bland_altman(df1 = df_D, df2 = df_K)

bland_altman_DdhD_D <- bland_altman(df1 = df_DdhD, df2 = df_D)
bland_altman_DBo_D <- bland_altman(df1 = df_DBo, df2 = df_D)
bland_altman_DBo_DdhD <- bland_altman(df1 = df_DBo, df2 = df_DdhD)
bland_altman_DBo_DdhD[1]
bland_altman_DBo_DRs <- bland_altman(df1 = df_DBo, df2 = df_DRs)
bland_altman_DBo_DRs[1]
bland_altman_DBo_DBWA <- bland_altman(df1 = df_DBo, df2 = df_DBWA)
bland_altman_DBo_DBWA[1]

bland_altman_KBo_K <- bland_altman(df1 = df_KBo, df2 = df_K)
bland_altman_KBo_K[1]
bland_altman_KBo_KBWA <- bland_altman(df1 = df_KBo, df2 = df_KBWA)
bland_altman_KBo_KBWA[1]
bland_altman_KBo_KRs <- bland_altman(df1 = df_KBo, df2 = df_KRs)
bland_altman_KBo_KRs[1]

i=1
for(i in 1:6) {

  plots <- list(
    bland_altman_DBo_KBo[[i]],
    bland_altman_DBWA_KBWA[[i]],
    bland_altman_DBWA_KBWA[[i]],
    bland_altman_DRs_KRs[[i]],
    bland_altman_D_K[[i]],

    #bland_altman_DdhD_D[[i]],
    bland_altman_DBo_D[[i]],
    bland_altman_DBo_DdhD[[i]],
    bland_altman_DBo_DRs[[i]],
    bland_altman_DBo_DBWA[[i]],

    bland_altman_KBo_K[[i]],
    bland_altman_KBo_KBWA[[i]],
    bland_altman_KBo_KRs[[i]]
  )

  # Extraer datos de cada gráfico
  data_list <- map(plots, ~ .x$data)
  all_data <- bind_rows(data_list)

  # Determinar los límites de los ejes
  x_min <- min(all_data$Promedio, na.rm = TRUE)
  x_max <- max(all_data$Promedio, na.rm = TRUE)
  y_min <- min(all_data$Diferencia, na.rm = TRUE)
  y_max <- max(all_data$Diferencia, na.rm = TRUE)


  plots_fixed <- map(plots, ~ .x +
                       coord_cartesian(xlim = c(x_min, x_max), ylim = c(y_min, y_max))
  )

  # Mostrar en un grid
  library(gridExtra)
  grid.arrange(grobs = plots_fixed, ncol = 4)

}




# VER DIFERENCIAS SIGNIFICATIVAS -------------
# DataFrame para almacenar los resultados del test de Wilcoxon y el grupo dominante

comparacion <- AR_genK_Bo

p_por_categoria <- data.frame("Categoria" = c(), "Genero" = c(), "Wilcoxon"= c(),"Wilcoxon-Paired" = c(),  "Dominante" = c())
i = 1

categoria <- colnames(comparacionT)[ncol(comparacionT)]
print(categoria)
for (j in 1:(ncol(comparacionT)-1)) {
  genero <- colnames(comparacionT)[j]
  df_sin_na <- comparacionT[complete.cases(comparacionT[, categoria]), ]

  # Excluir "Desconocido"
  if (any(df_sin_na[,categoria] == "Desconocido")) {
    df_sin_na <- df_sin_na[-which(df_sin_na[,categoria] == "Desconocido"),]
  }

  # Wilcoxon test
  w_test <- wilcox.test(df_sin_na[,genero] ~ df_sin_na[,categoria], df_sin_na)
  p_valor <- w_test$p.value
  w_testP <- wilcox.test(df_sin_na[,genero] ~ df_sin_na[,categoria], df_sin_na, paired = TRUE)
  p_valorP <- w_testP$p.value

  # Calcular el promedio de los valores para cada grupo
  promedio_kraken <- mean(df_sin_na[df_sin_na$Categoria == "KRAKEN", genero], na.rm = TRUE)
  promedio_dragen <- mean(df_sin_na[df_sin_na$Categoria == "DRAGEN", genero], na.rm = TRUE)

  # Determinar el grupo dominante
  dominante <- ifelse(promedio_kraken > promedio_dragen, "KRAKEN", "DRAGEN")

  # Almacenar resultados
  p_por_categoria[i,"Categoria"] <- categoria
  p_por_categoria[i,"Especie"] <- genero
  p_por_categoria[i,"Wilcoxon"] <- p_valor
  p_por_categoria[i,"Wilcoxon-Paired"] <- p_valorP
  p_por_categoria[i,"Dominante"] <- dominante
  i <- i+1
}

table(p_por_categoria$Dominante)

p_signP <- p_por_categoria[which(p_por_categoria$`Wilcoxon-Paired`<0.05),]
table(p_signP$Dominante)
species_dif_sign <- p_signP$Especie
