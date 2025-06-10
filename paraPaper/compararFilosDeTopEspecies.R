##########################################################################################################
#Cuáles son los filos de las especies más abundantes según cada metodología y dentro de cada categoría?
########################################################################################################

MetadataB <- as.data.frame(read_excel("~/Daniela/Biota/PipelineBiota/data/Metadata-soloColumnasUsables.xlsx"))

# Distribucion de frecuencia de filo proteobacteria en top 10 especies por rango etario ------------------

library(dplyr)
library(tidyr)

top_n_species <- 10
library(PipelineBiota)

extraerAR_Species_Filo <- function(patients_dir, source, de_host ) {
  tabla_otus <- generateOTUsTableGrupal(patients_dir = patients_dir, source = source, conEukaryota = FALSE, de_host = de_host)

  if(any(tabla_otus$Species == "Cutibacterium modestum")) {
    tabla_otus$Species[which(tabla_otus$Species == "Cutibacterium modestum")] <- "Propionibacterium sp. oral taxon 193"
    tabla_otus$Genus[which(tabla_otus$Species ==  "Propionibacterium sp. oral taxon 193")] <- "Propionibacterium"
  }

  out <- group_TaxonomicLevels(patients_dir = patients_dir, tabla_otus = tabla_otus, source = source, de_host = de_host, conE = FALSE)
  list_AR <- out[[1]]
  AR_filos <- list_AR[[3]]
  AR_generos <- list_AR[[7]]
  AR_species <- list_AR[[8]]
  colnames(AR_species)[-1] <- gsub(sprintf("_%s", source), "", colnames(AR_species)[-1])

  otus_unique <- tabla_otus[!duplicated(tabla_otus$Species), c("Species", "Phylum")]
  AR_species_Filo <- merge(AR_species, otus_unique, by = "Species", all.x = TRUE, all.y = FALSE)
  AR_species_Filo <- AR_species_Filo[-which(AR_species_Filo$Phylum == "-"),]

  return(AR_species_Filo)
}

#Quedarme con todas las muestras: sin separar por rango ni sexo
#K:
AR_speciesK_Bo <- extraerAR_Species_Filo(patients_dir = "~/Daniela/Biota/Muestras/73m", source = "KRAKEN", de_host = "Bowtie")
AR_speciesK_BWA <- extraerAR_Species_Filo(patients_dir = "~/Daniela/Biota/Muestras/73m", source = "KRAKEN", de_host = "BWA")
AR_speciesK_Rs <- extraerAR_Species_Filo(patients_dir = "~/Daniela/Biota/Muestras/73m", source = "KRAKEN", de_host = "RSubread" )
AR_speciesK_sin <- extraerAR_Species_Filo(patients_dir = "~/Daniela/Biota/Muestras/73m", source = "KRAKEN", de_host = "" )

#D:
AR_speciesD_BWA <- extraerAR_Species_Filo(patients_dir = "~/Daniela/Biota/Muestras/73m", source = "DRAGEN", de_host = "BWA" )
AR_speciesD_Bo <- extraerAR_Species_Filo(patients_dir = "~/Daniela/Biota/Muestras/73m", source = "DRAGEN", de_host = "Bowtie" )
AR_speciesD_Rs <- extraerAR_Species_Filo(patients_dir = "~/Daniela/Biota/Muestras/73m", source = "DRAGEN", de_host = "RSubread" )
AR_speciesD_sin <- extraerAR_Species_Filo(patients_dir = "~/Daniela/Biota/Muestras/73m", source = "DRAGEN", de_host = "" )
AR_speciesD_sinDH_PD <- extraerAR_Species_Filo(patients_dir = "~/Daniela/Biota/Muestras/73m", source = "DRAGEN", de_host = "sinDH_PD" )


#Todas las muestras:
df_list <- list(AR_speciesK_Bo, AR_speciesK_BWA, AR_speciesK_Rs, AR_speciesK_sin,
                AR_speciesD_Bo, AR_speciesD_BWA, AR_speciesD_Rs, AR_speciesD_sin,
                AR_speciesD_sinDH_PD)

# Asignar nombre de la metodología a cada dataframe de la lista
methodologies <- c("BoK", "bwaK", "RsK", "-K",
                   "BoD", "bwaD", "RsD", "DdhD", "-D")

df_list <- mapply(function(df, method) {
  df %>%
    pivot_longer(cols = -c(Species, Phylum), names_to = "Sample", values_to = "Abundance") %>%
    group_by(Sample) %>%
    slice_max(order_by = Abundance, n = top_n_species) %>%
    ungroup() %>%
    mutate(Methodology = method)
}, df_list, methodologies, SIMPLIFY = FALSE)

df1 <- df_list[[1]]
table(df1$Phylum)

df5 <- df_list[[5]]
table(df5$Phylum)

# Combinar todos los dataframes en uno solo
df_combined <- as.data.frame(bind_rows(df_list))
str(df_combined)
df_combined$Methodology <- as.factor(df_combined$Methodology)
df_combined$Sample <- as.factor(df_combined$Sample)
df_combined$Phylum <- as.factor(df_combined$Phylum)
df_combined$Abundance <- as.numeric(df_combined$Abundance)

# control:
#df <- df[order(df$'1', decreasing = TRUE),]
#all(df$Species[1:10] == df_top_phyla$Species[df_top_phyla$Sample == "1"])
# Definir los filos de interés

df_phylum_freq <- df_combined %>%
  filter(Phylum %in% c("Proteobacteria", "Firmicutes", "Actinobacteria")) %>%
  group_by(Methodology, Sample, Phylum) %>%
  summarise(Freq = n(), .groups = "drop")

df_phylum_freq <- df_phylum_freq %>%
  complete(Phylum = c("Proteobacteria", "Firmicutes", "Actinobacteria"),
           Methodology = unique(df_combined$Methodology),
           Sample = unique(df_combined$Sample),
           fill = list(Freq = 0))
df_phylum_freq <- df_phylum_freq %>%
  group_by(Methodology, Sample) %>%
  mutate(RelFreq = Freq / sum(Freq)) %>%
  ungroup()


# Este se usa:
df_phylum_freq <- df_combined %>%
  filter(Phylum %in% c("Proteobacteria", "Firmicutes", "Actinobacteria")) %>%  # Filtrar los filos de interés
  count(Methodology, Sample, Phylum, name = "Freq") %>%  # Contar la frecuencia de cada Phylum por metodología y muestra
  # Crear todas las combinaciones posibles de Phylum, Methodology y Sample
  complete(Phylum = c("Proteobacteria", "Firmicutes", "Actinobacteria"),
           Methodology = unique(df_combined$Methodology),
           Sample = unique(df_combined$Sample),
           fill = list(Freq = 0)) %>% # Completar con 0 cuando no haya datos
  group_by(Methodology, Sample) %>%  # Agrupar por metodología y muestra
  mutate(RelFreq = Freq / sum(Freq)) %>%  # Calcular la frecuencia relativa dentro de cada muestra
  ungroup()

#Agrego info de rango y sexo:
all(df_phylum_freq$Sample  %in% MetadataB$ID)
colnames(MetadataB)[1] <- "Sample"
df_phylum_freq <- merge(df_phylum_freq, MetadataB[, c("Sample", "Rango etario", "Sexo")], by = "Sample")


# Un boxplot por filo:
library(ggplot2)
library(RColorBrewer)

# Definir colores personalizados
custom_colors <- c(
  "Proteobacteria" = "#8470FF",  # Lila (de Set2)
  "Firmicutes" = "#EE9A00",  # Naranjita (de Set2)
  "Actinobacteria" = "#20B2AA"  # Azul verdoso (de Set2)
)

df_phylum_freq$`Rango etario` <- factor(df_phylum_freq$`Rango etario` , levels = c("18-35", "35-55",">55"))
ggplot(df_phylum_freq[df_phylum_freq$Phylum %in% c("Proteobacteria", "Firmicutes", "Actinobacteria"),],
       aes(x = Methodology, y = RelFreq, fill = Phylum)) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = "Metodología", y = "Frecuencia",
       #title = paste0("Frecuencia Relativa Filos en el top 10 esp - ", rangodet)) +
       title = "Frecuencia Relativa Filos en el top 10 esp ") +
  scale_y_continuous(limits = c(0, 1)) +  # Establecer límites del eje Y de 0 a 1
  scale_fill_manual(values = c("Proteobacteria" = "#8470FF", "Firmicutes" = "#FFA500", "Actinobacteria" = "#66CDAA")) +  # Colores para cada filo
  #facet_wrap(~ Phylum) +  # Crear facetas por filo
  #facet_wrap(~ `Rango etario`) +  # Crear facetas por filo
  #facet_wrap(~ Phylum + `Rango etario`, scales = "free_x")
  facet_wrap(~ Phylum + Sexo, scales = "free_x") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    plot.title = element_text(size = 12)
  )

# >>>>>> FIGURE 5 A <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#BOXPLOT DE FREC PROTEOBACTERIA POR GRUPO ETARIO:
colnames(df_phylum_freq)
df_proteobacteria <- df_phylum_freq[df_phylum_freq$Phylum == "Proteobacteria",]
df_proteobacteria$Rango <- factor(df_proteobacteria$`Rango etario`, levels = c("18-35", "35-55", ">55") )

kruskal_results <- df_proteobacteria %>%
  group_by(Methodology) %>%
  summarise(p_value = kruskal.test(RelFreq ~ Rango)$p.value) %>%
  mutate(p_text = paste0("p = ", signif(p_value, 2)))  # Formatear el p-valor

# Agregar columna TaxMethod según la presencia de "D" o "K"
df_proteobacteria <- df_proteobacteria %>%
  mutate(`Taxonomic Clasiffier` = ifelse(grepl("D", Methodology), "DRAGEN",
                            ifelse(grepl("K", Methodology), "Kraken", "Other")))

unique(df_proteobacteria$Methodology)
df_proteobacteria$Methodology <- factor(df_proteobacteria$Methodology,
                                        levels = c("-D", "BoD", "RsD", "bwaD","DdhD",
                                                   "-K", "BoK", "RsK", "bwaK"))

df_proteobacteria$`Taxonomic Clasiffier` <- as.factor(df_proteobacteria$`Taxonomic Clasiffier`)

colores <- c("DRAGEN" = "#FFA500", "Kraken" = "#8470FF")  # Naranja y violeta
colores <- c("DRAGEN" = "black", "Kraken" = "grey50")  # Naranja y violeta

ggplot(df_proteobacteria, aes(x = Rango, y = RelFreq, pattern = `Taxonomic Clasiffier`, fill = `Taxonomic Clasiffier`)) +
  geom_violin(alpha = 0.3, color = NA) +  # Violin plot con transparencia
  geom_jitter(aes(color = `Taxonomic Clasiffier`), width = 0.2, alpha = 0.5, size = 1) +  # Puntos con transparencia
  geom_boxplot(width = 0.4, outlier.shape = NA, color = "black") +  # Boxplot sin puntos outliers
  scale_fill_manual(values = colores) +
  scale_color_manual(values = colores) +
  scale_pattern_manual(values = patrones) +
  theme_minimal() +
  labs(x = "Methodology", y = "Frequency", title = "Proteobacteria in top 10 species by age group") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    plot.title = element_text(size = 12, face = "bold")
  ) +
  facet_wrap(~ Methodology, scales = "free_x", ncol = 5) +
  geom_text(data = kruskal_results,
            aes(x = 2, y = max(df_proteobacteria$RelFreq, na.rm = TRUE) + 0.05, label = p_text),
            inherit.aes = FALSE, size = 3)


#Texturas diferentes:
library(ggpattern)

# >>>>>>>>>>>>>>>> figure 5 A <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# Definimos los patrones para cada metodología
patrones <- c("DRAGEN" = "stripe", "Kraken" = "none")
df_proteobacteria$Pattern <- ifelse(df_proteobacteria$`Taxonomic Clasiffier` == "Kraken", "none", "stripe")
df_proteobacteria$`Taxonomic Clasiffier`

ggplot(df_proteobacteria, aes(x = Rango, y = RelFreq)) +
  geom_violin(alpha = 0.3, color = NA) +
  geom_jitter(aes(color = "black"), width = 0.2, alpha = 0.5, size = 1) +
  geom_boxplot_pattern(
    #aes(pattern = `Taxonomic Clasiffier`),
    #width = 0.5,
    aes(pattern = `Taxonomic Clasiffier`),
    position = position_dodge(width = 0.25),
    outlier.shape = NA,
    pattern_fill = "black",
    fill = "grey70",
    pattern_spacing = 0.05,
    color = "black"
  ) +
  scale_pattern_manual(values = patrones) +
  scale_color_manual(values = c("DRAGEN" = "#000000", "Kraken" = "#333333")) +
  theme_minimal() +
  labs(
    x = "Methodology",
    y = "Proteobacteria's Frequency",
    title = ""
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    plot.title = element_text(size = 12, face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    legend.box = "horizontal"
  ) +
  facet_wrap(~ Methodology, scales = "free_x", ncol = 5) +
  geom_text(
    data = kruskal_results,
    aes(x = 2, y = max(df_proteobacteria$RelFreq, na.rm = TRUE) + 0.05, label = p_text),
    inherit.aes = FALSE,
    size = 3
  )


#########################
wilcox_results <- df_phylum_freq %>%
  group_by(Phylum, Methodology) %>%
  summarise(p_value = wilcox.test(RelFreq ~ Sexo)$p.value) %>%
  mutate(p_text = paste0("p = ", signif(p_value, 3)))  # Formatear el p-valor

kruskal_results <- df_phylum_freq %>%
  group_by(Phylum, Methodology) %>%
  summarise(p_value = kruskal.test(RelFreq ~ `Rango etario`)$p.value) %>%
  mutate(p_text = paste0("p = ", signif(p_value, 3)))  # Formatear el p-valor



df_firm <- df_phylum_freq[df_phylum_freq$Phylum == "Firmicutes",]
df_firm$Rango <- factor(df_firm$`Rango etario`, levels = c("18-35", "35-55", ">55") )
ggplot(df_firm, aes(x = Sexo, y = RelFreq, fill = Sexo)) +
  #geom_boxplot() +
  geom_boxplot(fill ="#EE9A00") +  # Asignar color fijo a todos los boxplots
  theme_minimal() +
  labs(x = "Metodologia", y = "Frecuencia ", title = "Frec Firmicutes") +
  #scale_fill_manual(values = "#8470FF", aesthetics = "fill", na.value = "grey") +  # Colores personalizados
  scale_y_continuous(limits = c(0, 1)) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    plot.title = element_text(size = 10)
  ) +
  facet_wrap(~ Methodology, scales = "free_x", ncol= 5)

wilcox_results <- df_firm %>%
  group_by(Methodology) %>%
  summarise(p_value = wilcox.test(RelFreq ~ Sexo)$p.value) %>%
  mutate(p_text = paste0("p = ", signif(p_value, 3)))  # Formatear el p-valor


df_act <- df_phylum_freq[df_phylum_freq$Phylum == "Actinobacteria",]
df_act$Rango <- factor(df_act$`Rango etario`, levels = c("18-35", "35-55", ">55") )
ggplot(df_act, aes(x = Sexo, y = RelFreq, fill = Sexo)) +
  #geom_boxplot() +
  geom_boxplot(fill  = "#20B2AA") +  # Asignar color fijo a todos los boxplots
  theme_minimal() +
  labs(x = "Metodologia", y = "Frecuencia ", title = "Frec Actinobacteria") +
  #scale_fill_manual(values = "#8470FF", aesthetics = "fill", na.value = "grey") +  # Colores personalizados
  scale_y_continuous(limits = c(0, 1)) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    plot.title = element_text(size = 10)
  ) +
  facet_wrap(~ Methodology, scales = "free_x", ncol = 5)

# Calcular p-valores por cada metodología
p_values <- df_act %>%
  group_by(Methodology) %>%
  summarise(P_Value = wilcox.test(RelFreq ~ Sexo)$p.value)

# Calcular resumen estadístico por metodología y sexo
summary_stats <- df_act %>%
  group_by(Methodology, Sexo) %>%
  summarise(
    Q1 = quantile(RelFreq, 0.25, na.rm = TRUE),
    Median = median(RelFreq, na.rm = TRUE),
    Q3 = quantile(RelFreq, 0.75, na.rm = TRUE),
    .groups = "drop"
  )
summary_stats[, c("Q1", "Median", "Q3")] <- round(summary_stats[, c("Q1", "Median", "Q3")], 2)

# Quiero ver las especies más abundantes por sexo por software: ---------------------------------------------------

compararFilosDeTopEspecies(categoria = rangodet, AR_speciesKfem_Bo, AR_speciesKfem_BWA, AR_speciesKfem_Rs, AR_speciesKfem_sin,
                           AR_speciesDfem_Bo, AR_speciesDfem_BWA, AR_speciesDfem_Rs, AR_speciesDfem_sin,
                           AR_speciesDfem_sinDHPD)

compararFilosDeTopEspecies <- function(categoria, AR_speciesKfem_Bo, AR_speciesKfem_BWA, AR_speciesKfem_Rs, AR_speciesKfem_sin,
                                       AR_speciesDfem_Bo, AR_speciesDfem_BWA, AR_speciesDfem_Rs, AR_speciesDfem_sin, AR_speciesDfem_sinDHPD) {
  freq_KBo <- table(AR_speciesKfem_Bo$Phylum[1:10])
  freq_KBWA <- table(AR_speciesKfem_BWA$Phylum[1:10])
  freq_KRs <- table(AR_speciesKfem_Rs$Phylum[1:10])
  freq_Ksin <- table(AR_speciesKfem_sin$Phylum[1:10])

  freq_DBo <- table(AR_speciesDfem_Bo$Phylum[1:10])
  freq_DBWA <- table(AR_speciesDfem_BWA$Phylum[1:10])
  freq_DRs <- table(AR_speciesDfem_Rs$Phylum[1:10])
  freq_Dsin <- table(AR_speciesDfem_sin$Phylum[1:10])
  freq_DsinDHPD <- table(AR_speciesDfem_sinDHPD$Phylum[1:10])

  # Combinar las frecuencias en un dataframe

  freq_list <- list(KBo = freq_KBo, KBWA = freq_KBWA, KRs = freq_KRs, Ksin = freq_Ksin,
                    DBo = freq_DBo, DBWA = freq_DBWA, DRs = freq_DRs, Dsin = freq_Dsin,
                    DsinDH_PD = freq_DsinDHPD)

  # Obtener todos los filos únicos de todas las tablas
  all_phyla <- unique(unlist(lapply(freq_list, names)))
  df_frequencies <- data.frame(Phylum = all_phyla)

  # Llenar las columnas con las frecuencias de cada tabla
  for (name in names(freq_list)) {
    df_frequencies[[name]] <- as.numeric(freq_list[[name]][df_frequencies$Phylum])
  }


  # Reemplazar los NA (que aparecen cuando un phylum no está en una tabla) con 0
  df_frequencies[is.na(df_frequencies)] <- 0

  library(reshape2)
  phylum_counts_long <- melt(df_frequencies, id.vars = "Phylum",
                             variable.name = "Source", value.name = "Count")

  # Agregar una columna para agrupar Kraken vs Dragen
  phylum_counts_long$Group <- ifelse(grepl("K", phylum_counts_long$Source), "Kraken", "Dragen")

  phylum_counts_long$Proportion <- phylum_counts_long$Count /
    ave(phylum_counts_long$Count, phylum_counts_long$Source, FUN = sum)


  # Crear el gráfico
  # Crear el gráfico de barras apiladas
  library(ggplot2)
  phylum_counts_long$Source <- factor(phylum_counts_long$Source,
                                      levels = c(sort(unique(phylum_counts_long$Source[phylum_counts_long$Group == "Kraken"])),
                                                 sort(unique(phylum_counts_long$Source[phylum_counts_long$Group == "Dragen"]))))

  ggplot(phylum_counts_long, aes(x = Source, y = Proportion, fill = Phylum)) +
    geom_bar(stat = "identity") +
    labs(title = paste0("Filos de top 10 Especies -", categoria),
         x = "Metodología", y = "Proporción del Total") +
    scale_fill_brewer(palette = "Set2") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_continuous(labels = scales::percent)  +# Mostrar las proporciones en porcentaje
    facet_grid(~ Group, scales = "free_x", space = "free_x")

}

# Test de comparacion de proporciones: -----------------------

proteobacteria_data <- phylum_counts_long %>%
  filter(Phylum == "Proteobacteria")

# Crear tabla de contingencia: filas = Group, columnas = conteo de Proteobacteria
contingency_table <- proteobacteria_data %>%
  group_by(Group) %>%
  summarise(Total = sum(Count)) %>%
  pivot_wider(names_from = Group, values_from = Total) %>%
  as.matrix()

# Test de chi-cuadrado
chi_result <- chisq.test(contingency_table)
print(chi_result)

# Test de Fisher (si hay valores bajos en la tabla) : NO ME ANDA
contingency_table <- proteobacteria_data %>%
  mutate(Presence = ifelse(Count > 0, "Present", "Absent")) %>%
  group_by(Source, Presence) %>%
  summarise(Total = n(), .groups = "drop") %>%
  pivot_wider(names_from = Presence, values_from = Total, values_fill = list(Total = 0)) %>%
  as.matrix()

# Verificar la tabla de contingencia
print(contingency_table)

# Aplicar test de Fisher
fisher_result <- fisher.test(contingency_table)
print(fisher_result)

