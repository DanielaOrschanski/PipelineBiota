library(PipelineBiota)
library(dplyr)
library(readr)
library(tidyr)
library(openxlsx)
patients_dir <- "~/Daniela/Biota/Muestras/73m"
setwd("~/Daniela/Biota/Muestras/73m")
length(list.dirs("~/Daniela/Biota/Muestras/73m", recursive = FALSE, full.names = FALSE))
Bo_vias <- generatePathwaysTable(patients_dir = "~/Daniela/Biota/Muestras/73m", de_host = "Bowtie")
Bo_vias_CPM <- Bo_vias[[1]]
Bo_vias_AR <- Bo_vias[[2]]
Bo_vias_org <- Bo_vias[[3]]

RunKRAKEN(patients_dir = patients_dir, de_host = "BWA")
tabla_otus <- generateOTUsTableGrupal(patients_dir = patients_dir, source = "KRAKEN", conEukaryota = TRUE, de_host = "Bowtie")
out <- group_TaxonomicLevels(patients_dir = patients_dir, tabla_otus = tabla_otus, source = "KRAKEN", de_host = "BWA", conE = FALSE)
list_AR <- out[[1]]
AR_filos <- list_AR[[3]]
AR_generos <- list_AR[[7]]
AR_species <- list_AR[[8]]
colnames(AR_species)[-1] <- gsub(sprintf("_%s", "KRAKEN"), "", colnames(AR_species)[-1])

colSums(AR_species[,-1])


counts <- counts_Tax(patients_dir = patients_dir, source = "KRAKEN", de_host = "Bowtie", conEukaryota = TRUE)
counts <- counts_Tax(patients_dir = patients_dir, source = "KRAKEN", de_host = "BWA", conEukaryota = TRUE)

dehosters <- c("Bowtie", "BWA", "RSubread", "")
for(de_host in dehosters) {
  print(de_host)
  tabla_otusK <- generateOTUsTableGrupal(patients_dir = patients_dir, source = "KRAKEN", conEukaryota = TRUE, de_host = de_host)
  group_TaxonomicLevels(patients_dir = patients_dir, tabla_otus = tabla_otusK, source = "KRAKEN", de_host = de_host, conE = FALSE)

  tabla_otusD <- generateOTUsTableGrupal(patients_dir = patients_dir, source = "DRAGEN", conEukaryota = TRUE, de_host = de_host)
  group_TaxonomicLevels(patients_dir = patients_dir, tabla_otus = tabla_otusD, source = "DRAGEN", de_host = de_host, conE = FALSE)

}

library(vegan)
#K:
AR_speciesKBo <- read_excel("~/Daniela/Biota/Muestras/73m/AR_Bo__Species_KRAKEN.xlsx")
diversidad_KBo <- calculateAlphaDiversity(AR_species = AR_speciesKBo)

AR_speciesKBWA <- read_excel("~/Daniela/Biota/Muestras/73m/AR_bwa__Species_KRAKEN.xlsx")
diversidad_KBWA <- calculateAlphaDiversity(AR_species = AR_speciesKBWA)

AR_speciesKRs <- read_excel("~/Daniela/Biota/Muestras/73m/AR_Rs__Species_KRAKEN.xlsx")
diversidad_KRs <- calculateAlphaDiversity(AR_species = AR_speciesKRs)

AR_speciesKsin <- read_excel("~/Daniela/Biota/Muestras/73m/AR_sin__Species_KRAKEN.xlsx")
diversidad_Ksin <- calculateAlphaDiversity(AR_species = AR_speciesKsin)

#D:
AR_speciesDBo <- read_excel("~/Daniela/Biota/Muestras/73m/AR_Bo__Species_DRAGEN.xlsx")
diversidad_DBo <- calculateAlphaDiversity(AR_species = AR_speciesDBo)

AR_speciesDBWA <- read_excel("~/Daniela/Biota/Muestras/73m/AR_bwa__Species_DRAGEN.xlsx")
diversidad_DBWA <- calculateAlphaDiversity(AR_species = AR_speciesDBWA)

AR_speciesDRs <- read_excel("~/Daniela/Biota/Muestras/73m/AR_Rs__Species_DRAGEN.xlsx")
diversidad_DRs <- calculateAlphaDiversity(AR_species = AR_speciesDRs)

AR_speciesDsin <- read_excel("~/Daniela/Biota/Muestras/73m/AR_sin__Species_DRAGEN.xlsx")
diversidad_Dsin <- calculateAlphaDiversity(AR_species = AR_speciesDsin)

colnames(diversidad_Dsin)


library(dplyr)
library(ggplot2)
library(tidyr)

# Crear columna 'ID' y 'Method' para cada tabla
diversidad_KBo <- diversidad_KBo %>%
  mutate(ID = gsub("_KRAKEN", "", Sample), Method = "KRAKEN_Bo")
diversidad_KBWA <- diversidad_KBWA %>%
  mutate(ID = gsub("_KRAKEN", "", Sample), Method = "KRAKEN_BWA")
diversidad_KRs <- diversidad_KRs %>%
  mutate(ID = gsub("_KRAKEN", "", Sample), Method = "KRAKEN_Rs")
diversidad_Ksin <- diversidad_Ksin %>%
  mutate(ID = gsub("_KRAKEN", "", Sample), Method = "KRAKEN_sin")

diversidad_DBo <- diversidad_DBo %>%
  mutate(ID = gsub("_DRAGEN", "", Sample), Method = "DRAGEN_Bo")
diversidad_DBWA <- diversidad_DBWA %>%
  mutate(ID = gsub("_DRAGEN", "", Sample), Method = "DRAGEN_BWA")
diversidad_DRs <- diversidad_DRs %>%
  mutate(ID = gsub("_DRAGEN", "", Sample), Method = "DRAGEN_Rs")
diversidad_Dsin <- diversidad_Dsin %>%
  mutate(ID = gsub("_DRAGEN", "", Sample), Method = "DRAGEN_sin")

# Combinar todas las tablas
diversidad_combined <- bind_rows(
  diversidad_KBo, diversidad_KBWA, diversidad_KRs, diversidad_Ksin,
  diversidad_DBo, diversidad_DBWA, diversidad_DRs, diversidad_Dsin
)

# Transformar a formato largo
diversidad_long <- diversidad_combined %>%
  pivot_longer(
    cols = c(Shannon, Simpson, Chao1),
    names_to = "Index",
    values_to = "Value"
  )

# Crear el gráfico
ggplot(diversidad_long, aes(x = Method, y = Value, fill = Method)) +
  geom_boxplot(alpha = 0.7, outlier.color = "red") +
  geom_line(aes(group = interaction(ID, Index)), color = "gray50", linetype = "dashed", alpha = 0.5) +
  geom_point(position = position_dodge(width = 0.75), alpha = 0.9) +
  labs(
    title = "Comparación de índices entre métodos de KRAKEN y DRAGEN",
    x = "Método", y = "Valor"
  ) +
  facet_wrap(~Index, scales = "free_y") +  # Ejes Y independientes para cada índice
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top") +
  scale_fill_brewer(palette = "Set3")

# Comparar los valores de cada índice entre los métodos (Wilcoxon test)
p_values <- diversidad_combined %>%
  pivot_longer(cols = c(Shannon, Simpson, Chao1), names_to = "Index", values_to = "Value") %>%
  group_by(Index, Method) %>%
  summarize(
    median = median(Value, na.rm = TRUE),
    p_value = wilcox.test(Value ~ Method, paired = FALSE)$p.value
  )




##########################################################################
library(dplyr)
library(tidyr)
library(ggplot2)

# Crear columna 'ID' para extraer el identificador común
diversidad_KRAKEN <- diversidad_KRAKEN %>%
  mutate(ID = gsub("_KRAKEN", "", Sample),
         Method = "KRAKEN")

diversidad_DRAGEN <- diversidad_DRAGEN %>%
  mutate(ID = gsub("_DRAGEN", "", Sample),
         Method = "DRAGEN")

# Combinar ambas tablas según 'ID'
diversidad_combined <- full_join(
  diversidad_KRAKEN,
  diversidad_DRAGEN,
  by = "ID",
  suffix = c("_KRAKEN", "_DRAGEN")
)

# Transformar la tabla a formato largo
diversidad_long <- diversidad_combined %>%
  pivot_longer(
    cols = c(Shannon_KRAKEN, Simpson_KRAKEN, Chao1_KRAKEN,
             Shannon_DRAGEN, Simpson_DRAGEN, Chao1_DRAGEN),
    names_to = c("Index", "Method"),
    names_sep = "_",
    values_to = "Value"
  )


ggplot(diversidad_long, aes(x = Method, y = Value, fill = Method)) +
  geom_boxplot(alpha = 0.7, outlier.color = "red") +
  geom_line(aes(group = interaction(ID, Index)), color = "gray50", linetype = "dashed", alpha = 0.5) +
  geom_point(position = position_dodge(width = 0.75), alpha = 0.9) +
  labs(
    title = "Comparación de índices entre KRAKEN y DRAGEN",
    x = "Método", y = "Valor"
  ) +
  facet_wrap(~Index, scales = "free_y") +  # Ejes Y independientes para cada índice
  theme_minimal() +
  theme(legend.position = "top")

# Comparar los valores de cada índice entre KRAKEN y DRAGEN (Wilcoxon test)
p_values <- diversidad_combined %>%
  select(ID, Shannon_KRAKEN, Shannon_DRAGEN, Simpson_KRAKEN, Simpson_DRAGEN, Chao1_KRAKEN, Chao1_DRAGEN) %>%
  pivot_longer(cols = -ID, names_to = c("Index", "Method"), names_sep = "_", values_to = "Value") %>%
  spread(key = "Method", value = "Value") %>%
  group_by(Index) %>%
  summarize(
    p_value = wilcox.test(KRAKEN, DRAGEN, paired = FALSE)$p.value
  )

calculateAlphaDiversity <- function(AR_species) {
  AR <- as.matrix(t(AR_species[,-1]))
  shannon <- diversity(AR, index = "shannon")
  #Un índice de Shannon alto para una muestra indica que la comunidad es rica en especies y/o tiene una distribución más uniforme.
  simpson <- diversity(AR, index = "simpson")

  # Índice Chao1
  # Convertir a abundancias absolutas (necesario para Chao1)
  abund_abs <- round(AR * 100000)  # Ajusta 10000 según el rango de tus valores
  chao1 <- estimateR(abund_abs)["S.chao1", ]  # Chao1 extraído

  diversity_indices <- data.frame(
    Sample = colnames(AR_species[,-1]),
    Shannon = shannon,
    Simpson = simpson
    ,Chao1 = chao1
  )

  return(diversity_indices)
}

