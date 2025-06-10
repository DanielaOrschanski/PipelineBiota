#K:
AR_speciesK_Bo <- extraerAR_Species(patients_dir = "~/Daniela/Biota/Muestras/73m", source = "KRAKEN", de_host = "Bowtie")
AR_speciesK_BWA <- extraerAR_Species(patients_dir = "~/Daniela/Biota/Muestras/73m", source = "KRAKEN", de_host = "BWA")
AR_speciesK_Rs <- extraerAR_Species(patients_dir = "~/Daniela/Biota/Muestras/73m", source = "KRAKEN", de_host = "RSubread" )
AR_speciesK_sin <- extraerAR_Species(patients_dir = "~/Daniela/Biota/Muestras/73m", source = "KRAKEN", de_host = "" )

#D:
AR_speciesD_BWA <- extraerAR_Species(patients_dir = "~/Daniela/Biota/Muestras/73m", source = "DRAGEN", de_host = "BWA" )
AR_speciesD_Bo <- extraerAR_Species(patients_dir = "~/Daniela/Biota/Muestras/73m", source = "DRAGEN", de_host = "Bowtie" )
AR_speciesD_Rs <- extraerAR_Species(patients_dir = "~/Daniela/Biota/Muestras/73m", source = "DRAGEN", de_host = "RSubread" )
AR_speciesD_sin <- extraerAR_Species(patients_dir = "~/Daniela/Biota/Muestras/73m", source = "DRAGEN", de_host = "" )
AR_speciesD_sinDH_PD <- extraerAR_Species(patients_dir = "~/Daniela/Biota/Muestras/73m", source = "DRAGEN", de_host = "sinDH_PD" )

df_list <- list(AR_speciesK_Bo, AR_speciesK_BWA, AR_speciesK_Rs, AR_speciesK_sin,
                AR_speciesD_Bo, AR_speciesD_BWA, AR_speciesD_Rs, AR_speciesD_sin,
                AR_speciesD_sinDH_PD)

methodologies <- c("KBo", "KBWA", "KRs", "K",
                   "DBo", "DBWA", "DRs", "DdhD", "D")

df_list <- mapply(function(df, method) {
  df %>%
    pivot_longer(cols = -c(Species, Phylum), names_to = "Sample", values_to = "Abundance") %>%
    mutate(Methodology = method)  # Añadir columna de metodología
}, df_list, methodologies, SIMPLIFY = FALSE)

# Combinar todos los dataframes en uno solo
df_combined <- as.data.frame(bind_rows(df_list))
str(df_combined)
df_combined$Methodology <- as.factor(df_combined$Methodology)
df_combined$Sample <- as.factor(df_combined$Sample)
df_combined$Phylum <- as.factor(df_combined$Phylum)
df_combined$Abundance <- as.numeric(df_combined$Abundance)

#Agrego info de rango y sexo:
MetadataB <- as.data.frame(read_excel("~/Daniela/Biota/PipelineBiota/data/Metadata-soloColumnasUsables.xlsx"))
all(df_combined$Sample  %in% MetadataB$ID)
colnames(MetadataB)[1] <- "Sample"
df_combined <- merge(df_combined, MetadataB[, c("Sample", "Rango etario", "Sexo")], by = "Sample")
df_combined$Sexo <- as.factor(df_combined$Sexo)
df_combined$Rango <- factor(df_combined$`Rango etario`, levels = c("18-35", "35-55", ">55"))
str(df_combined)

library(vegan)

  #  índices de diversidad
diversity_indices <- df_combined %>%
  group_by(Methodology, Sample, Sexo, Rango) %>%
  summarise(
    Shannon = diversity(Abundance, index = "shannon"),
    Simpson = diversity(Abundance, index = "simpson"),
    Richness = specnumber(Abundance),
    #ACE = estimateR(as.integer(Abundance))
  )

df_long <- diversity_indices %>%
  pivot_longer(cols = c(Shannon, Simpson, Richness, ACE), names_to = "Index", values_to = "Abundance")
df_long$Source <- ifelse(grepl("K", df_long$Methodology), "Kraken 2", "DRAGEN")

ggplot(df_long, aes(x = Methodology, y = Abundance, fill = Source)) +
  geom_boxplot() +
  facet_wrap(~ Index, scales = "free_y") +
  theme_minimal() +
  labs(title = "Alpha Diversity by Methodology",
       x = "Methodology",
       y = "Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Generar todas las comparaciones por pares de metodologías
metodologias <- unique(df_long$Methodology)
comparaciones <- combn(metodologias, 2, simplify = FALSE)

# Calcular p-valores manualmente con wilcox.test() para cada índice
p_values <- df_long %>%
  group_by(Index) %>%
  summarise(
    Comparacion = list(sapply(comparaciones, function(pair) paste0(pair[1], " vs ", pair[2]))),
    p_value = list(sapply(comparaciones, function(pair) {
      grupo1 <- df_long %>% filter(Methodology == pair[1]) %>% pull(Abundance)
      grupo2 <- df_long %>% filter(Methodology == pair[2]) %>% pull(Abundance)
      wilcox.test(grupo1, grupo2)$p.value
    }))
  ) %>%
  unnest(cols = c(Comparacion, p_value))  # Expandir las listas en filas separadas



p_values <- df_long %>%
  group_by(Methodology, Index) %>%
  summarise(
    p_value = wilcox.test(Abundance ~ Sexo)$p.value,
    .groups = "drop"  # Para evitar advertencias sobre el agrupamiento
  )


p_values <- df_long %>%
  group_by(Methodology, Index) %>%
  summarise(
    p_value = kruskal.test(Abundance ~ Rango)$p.value,
    .groups = "drop"  # Para evitar advertencias sobre el agrupamiento
  )

df_young_old <- df_long[df_long$Rango %in% c("18-35", ">55"),]
p_values <- df_young_old %>%
  group_by(Methodology, Index) %>%
  summarise(
    p_value = wilcox.test(Abundance ~ Rango)$p.value,
    .groups = "drop"  # Para evitar advertencias sobre el agrupamiento
  )


#DIVERSIDAD BETA:
library(vegan)
abundance_matrices <- list()

# Crear la matriz de abundancia por cada metodología
df_combined$`Rango etario`

abundance_matrices <- df_combined %>%
  group_by(Methodology) %>%
  group_split() %>%
  lapply(function(group_data) {
    group_data %>%
      pivot_wider(names_from = "Species", values_from = "Abundance", values_fill = list(Abundance = 0)) %>%
      select(-c(Methodology, Sample, Sexo, Rango,`Rango etario`, Phylum))  # Eliminar columnas no relevantes
  })

m <- abundance_matrices[[1]]
bray_curtis_matrices <- lapply(abundance_matrices, function(matrix) {
  vegdist(matrix, method = "bray")
})

# Convertir cada matriz de disimilitud a un dataframe
bray_curtis_df <- lapply(bray_curtis_matrices, function(matrix) {
  as.data.frame(as.matrix(matrix))
})

# Realizar un análisis MDS para Bray-Curtis
mds_bray <- lapply(bray_curtis_matrices, metaMDS)

# Graficar el resultado de MDS para Bray-Curtis por cada metodología
par(mfrow = c(1, length(mds_bray)))  # Configura los gráficos para que se muestren en una fila
lapply(mds_bray, function(mds) {
  plot(mds, main = "Diversidad Beta - Bray-Curtis")
})
