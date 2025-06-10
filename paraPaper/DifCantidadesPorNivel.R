
#patients_dir = "~/Daniela/Biota/Muestras/73m"
#source = "KRAKEN"
#de_host = "Bowtie"


extraerCantidadesPorNivel <- function(patients_dir, source, de_host ) {

  cant_nivel <- counts_Tax(patients_dir = patients_dir, source = source, de_host = de_host, conEukaryota = FALSE)
  met <- ifelse(source == "KRAKEN", "K", "D")

  if(de_host  == "Bowtie") {
    de_host_file <- "Bo"
  } else if( de_host == "BWA") {
    de_host_file <- "bwa"
  } else if(de_host == "RSubread") {
    de_host_file <- "Rs"
  } else if(de_host == "") {
    de_host_file <- "sin"
  } else if(de_host == "sinDH_PD") { # sin de host previo ni dehost dragen
    de_host_file <- "sinDH_PD"
  } else {
    stop("de_host must be Bowtie, BWA, RSubread, sinDH_PD or empty string")
  }

  metodologia <- paste0(met, de_host_file, sep="")
  cant_nivel$Metodologia <- metodologia
  #cant_nivelT <- as.data.frame(t(cant_nivel))
  #colnames(cant_nivelT) <- cant_nivelT[1,]
  #cant_nivelT <- cant_nivelT[-1,]
  return(cant_nivel)

}

library(data.table)
library(readxl)
library(openxlsx)
cantKBo <- extraerCantidadesPorNivel(patients_dir = "~/Daniela/Biota/Muestras/73m", source = "KRAKEN", de_host = "Bowtie")
cantKBWA <- extraerCantidadesPorNivel(patients_dir = "~/Daniela/Biota/Muestras/73m", source = "KRAKEN", de_host = "BWA")
cantKRs <- extraerCantidadesPorNivel(patients_dir = "~/Daniela/Biota/Muestras/73m", source = "KRAKEN", de_host = "RSubread")
cantKsin <- extraerCantidadesPorNivel(patients_dir = "~/Daniela/Biota/Muestras/73m", source = "KRAKEN", de_host = "")

cantDBo <- extraerCantidadesPorNivel(patients_dir = "~/Daniela/Biota/Muestras/73m", source = "DRAGEN", de_host = "Bowtie")
cantDBWA <- extraerCantidadesPorNivel(patients_dir = "~/Daniela/Biota/Muestras/73m", source = "DRAGEN", de_host = "BWA")
cantDRs <- extraerCantidadesPorNivel(patients_dir = "~/Daniela/Biota/Muestras/73m", source = "DRAGEN", de_host = "RSubread")
cantDsin <- extraerCantidadesPorNivel(patients_dir = "~/Daniela/Biota/Muestras/73m", source = "DRAGEN", de_host = "")
cantDsinDH_PD <- extraerCantidadesPorNivel(patients_dir = "~/Daniela/Biota/Muestras/73m", source = "DRAGEN", de_host = "sinDH_PD")


ConteosKD <- rbind(cantKBo, cantKBWA, cantKRs, cantKsin,
                   cantDBo, cantDBWA, cantDRs, cantDsin, cantDsinDH_PD)

ConteosKD_largo <- ConteosKD %>%
  pivot_longer(cols = c(Phylum, Genus, Species),
               names_to = "Niveles_Taxonomicos",
               values_to = "Conteos")

# Crear el boxplot
str(ConteosKD_largo)
ConteosKD_largo$Conteos <- as.numeric(ConteosKD_largo$Conteos)
orden_niveles <- colnames(ConteosKD)[2:10]
ConteosKD_largo$Niveles_Taxonomicos <- factor(ConteosKD_largo$Niveles_Taxonomicos, levels = orden_niveles)
unique(ConteosKD_largo$Metodologia)
ConteosKD_largo$Metodologia[ConteosKD_largo$Metodologia == "Ksin"] <- "-K"
ConteosKD_largo$Metodologia[ConteosKD_largo$Metodologia == "Dsin"] <- "DdhD"
ConteosKD_largo$Metodologia[ConteosKD_largo$Metodologia == "DsinDH_PD"] <- "-D"
ConteosKD_largo$Metodologia[ConteosKD_largo$Metodologia == "DsinDH_PD"] <- "-D"
ConteosKD_largo$Metodologia[ConteosKD_largo$Metodologia == "KBo"] <- "BoK"
ConteosKD_largo$Metodologia[ConteosKD_largo$Metodologia == "Kbwa"] <- "bwaK"
ConteosKD_largo$Metodologia[ConteosKD_largo$Metodologia == "KRs"] <- "RsK"

ConteosKD_largo$Metodologia[ConteosKD_largo$Metodologia == "DBo"] <- "BoD"
ConteosKD_largo$Metodologia[ConteosKD_largo$Metodologia == "DRs"] <- "RsD"
ConteosKD_largo$Metodologia[ConteosKD_largo$Metodologia == "Dbwa"] <- "bwaD"

ConteosKD_largo$Metodologia <- factor(ConteosKD_largo$Metodologia, levels = c("-D", "DdhD", "BoD", "RsD", "bwaD",
                                                                              "-K", "BoK", "RsK", "bwaK"))
ConteosKD_largo <- ConteosKD_largo %>%
  mutate(GrupoColor = ifelse(grepl("D", Metodologia), "DRAGEN", "Kraken"))

colores_metodologia <- c("Kraken" = "#FF7F00",  # Azul
                         "DRAGEN" = "#9370DB")  # Verde

# Crear el boxplot con líneas que conecten la misma muestra en distintas metodologías

#mas estetico
ggplot(ConteosKD_largo, aes(x = Metodologia, y = Conteos, fill = GrupoColor)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +  # Boxplots agrupados
  geom_line(aes(group = interaction(ID, Niveles_Taxonomicos)),
            colour = "black", linetype = "11", size = 0.3) +  # Líneas finas conectando puntos

  theme_minimal() +
  labs(x = "Taxonomic Level", y = "Number of Microorganisms",
       title = "Microorganisms Identified by Different Methodologies") +

  theme(
    text = element_text(size = 15),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    axis.text.x = element_text(angle = 90, hjust = 0.5),
    axis.title.x = element_text(margin = margin(t = 15)),  # Separar título eje x de los valores
    axis.title.y = element_text(margin = margin(r = 15)),  # Separar título eje y de los valores
    legend.position = "bottom",  # **Ubica la leyenda abajo**
    legend.box = "horizontal"  # **Hace que los elementos de la leyenda se distribuyan horizontalmente**
  ) +

  scale_fill_manual(values = colores_metodologia, name = "Taxonomic Classifier") +
  facet_wrap(~ Niveles_Taxonomicos, scales = "free_y")  # Crear gráficos por cada nivel taxonómico con escalas independientes


# >>>>>>>>>>>>>>>> FIGURE 3A <<<<<<<<<<<<<<<<<<<
#Patrones en vez de colores:
library(ggpattern)
patrones_metodologia <- c("DRAGEN" = "stripe", "KRAKEN" = "fill")

ggplot(ConteosKD_largo, aes(x = Metodologia, y = Conteos, pattern = GrupoColor)) +
  geom_boxplot_pattern(
    aes(group = interaction(Metodologia, GrupoColor)),
    pattern_fill = "black",
    pattern_density = 0.4,
    pattern_spacing = 0.05,
    pattern_angle = 45,
    color = "black",
    alpha = 0.7,
    position = position_dodge(width = 0.8)
  ) +
  geom_line(
    aes(group = interaction(ID, Niveles_Taxonomicos)),
    colour = "black", linetype = "11", size = 0.2
  ) +

  theme_minimal() +
  labs(x = "Taxonomic Level", y = "Number of Microorganisms",
       title = "") +

  theme(
    text = element_text(size = 15),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    axis.text.x = element_text(angle = 90, hjust = 0.5),
    axis.title.x = element_text(margin = margin(t = 15)),
    axis.title.y = element_text(margin = margin(r = 15)),
    legend.position = "bottom",
    legend.box = "horizontal"
  ) +

  scale_pattern_manual(values = patrones_metodologia, name = "Taxonomic Classifier") +
  facet_wrap(~ Niveles_Taxonomicos, scales = "free_y")




str(ConteosKD_largo)




pvalores <- ConteosKD_largo %>%
  group_by(Niveles_Taxonomicos) %>%
  summarise(
    p_value = kruskal.test(Conteos ~ Metodologia)$p.value
  )
print(pvalores)

library(dplyr)
library(tidyr)

#  p-valores de Wilcoxon por pares
calcular_pvalores_wilcoxon <- function(data) {
  metodologias <- unique(data$Metodologia)
  combinaciones <- combn(metodologias, 2, simplify = FALSE)

  resultados <- lapply(combinaciones, function(par) {
    test <- wilcox.test(
      Conteos ~ Metodologia,
      data = data %>% filter(Metodologia %in% par)
    )

    tibble(
      Niveles_Taxonomicos = unique(data$Niveles_Taxonomicos),
      Metodologia1 = par[1],
      Metodologia2 = par[2],
      p_value = test$p.value
    )
  })

  bind_rows(resultados)
}

# Aplicar la función a cada nivel taxonómico
pvalores_pares <- ConteosKD_largo %>%
  group_split(Niveles_Taxonomicos) %>%
  lapply(calcular_pvalores_wilcoxon) %>%
  bind_rows()

print(pvalores_pares)

write.xlsx(pvalores_pares, file = "~/Daniela/Biota/PipelineBiota/paraPaper/Pvalores_Metodologias_NumberMicroog.xlsx")

# PARA VIAS----------------------------
library(PipelineBiota)
source("~/Daniela/Biota/PipelineBiota/R/RunHuman.R", echo=TRUE)
library(readr)
cantBo <- count_Vias(patients_dir = "~/Daniela/Biota/Muestras/73m", de_host = "Bowtie")
cantBWA <- count_Vias(patients_dir = "~/Daniela/Biota/Muestras/73m", de_host = "BWA")
cantRs <- count_Vias(patients_dir = "~/Daniela/Biota/Muestras/73m",  de_host = "RSubread")
cantsin <- count_Vias(patients_dir = "~/Daniela/Biota/Muestras/73m", de_host = "")

Cantidades_vias <- rbind(cantBo, cantBWA, cantRs, cantsin)
unique(Cantidades_vias$de_host)
Cantidades_vias$de_host[Cantidades_vias$de_host == "Bo"] <- "BoH"
Cantidades_vias$de_host[Cantidades_vias$de_host == "Rs"] <- "RsH"
Cantidades_vias$de_host[Cantidades_vias$de_host == "T"] <- "-H"
Cantidades_vias$de_host[Cantidades_vias$de_host == "bwa"] <- "bwaH"


colores_de_host <- c("BoH"= "#E69F00",
                     "bwaH" = "#A52A2A",
                     "RsH" = "#009E73",
                     "-H" = "#9467BD" )

ggplot(Cantidades_vias, aes(x = de_host, y = Vias, fill = de_host)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = NA) +  # Violin plot con transparencia
  geom_jitter(aes(color = de_host), width = 0.2, size = 1.2, alpha = 0.7) +  # Jitter con mismo color que boxplot
  geom_boxplot(width = 0.6, outlier.shape = NA, color = "black", alpha = 0.8) +  # Boxplot sin outliers visibles
  geom_line(aes(group = ID), colour = "black", linetype = "11", alpha = 0.8) +  # Conectar mismo ID

  theme_minimal(base_size = 15) +
  labs(x = "",
       y = "Number of Pathways",
       title = "Number of Identified Pathways by each Methodology") +  # Título en negrita
  theme(
    plot.title = element_text(size = 14, face = "bold"),  # Título en negrita
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "none",       # Mueve la leyenda abajo
    legend.direction = "vertical"   # Ordena la leyenda en una fila
  ) +
  scale_fill_manual(values = colores_de_host, name = "De-Host") +
  scale_color_manual(values = colores_de_host, guide = "none")  # Usa los mismos colores para jitter y oculta la leyenda extra

pairwise.wilcox.test(Cantidades_vias$Vias, Cantidades_vias$de_host)

#Con patrones:
Cantidades_vias$de_host
ggplot(Cantidades_vias, aes(x = de_host, y = Vias, pattern = de_host)) +
  geom_violin(trim = FALSE, alpha = 0.4, color = "grey15") +  # Gráfico de violín sin patrones
  #geom_jitter(color = "grey20", width = 0.2, size = 1.2, alpha = 0.7) +  # Jitter para dispersión
  geom_boxplot_pattern(
    width = 0.6, outlier.shape = NA,
    color = "black", alpha = 0.8,
    pattern_density = 0.1
  ) +  # Boxplot con patrones
  geom_line(aes(group = ID), colour = "black", linetype = "11", alpha = 0.3) +
  theme_minimal() +
  labs(x = "", y = "% Mapped Reads", title = "") +
  scale_pattern_manual(values = c(
    "BoH" = "stripe",  # Patrones de rayas
    "bwaH" = "crosshatch",  # Patrones de círculos
    "RsH" = "circle",  # Patrones de líneas cruzadas
    "-H" = "none"  # Patrones de cuadrados
  )) +
  labs(x = "",
       y = "Number of Pathways",
       title = "") +  # Título en negrita
  theme(
    plot.title = element_text(size = 1, face = "bold"),  # Título en negrita
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "none"
  ) +
  scale_fill_manual(values = colores_de_host, name = "De-Host") +
  scale_color_manual(values = colores_de_host, guide = "none")  # Usa los mismos colores para jitter y oculta la leyenda extra


#DIAGRAMA DE VENN ---------------------------------
T_vias <- generatePathwaysTable(patients_dir = "~/Daniela/Biota/Muestras/73m", de_host = "")
T_vias_CPM <- T_vias[[1]]
T_vias_AR <- T_vias[[2]]

Bo_vias <- generatePathwaysTable(patients_dir = "~/Daniela/Biota/Muestras/73m", de_host = "Bowtie")
Bo_vias_CPM <- Bo_vias[[1]]
Bo_vias_AR <- Bo_vias[[2]]
Bo_vias_AR$SumaVia <- rowSums(Bo_vias_AR[,-c(1:3)])
Bo_vias_AR <- Bo_vias_AR[order(Bo_vias_AR$SumaVia, decreasing = FALSE),]
Bo_vias_AR[1,]
all(rowSums(Bo_vias_AR[,-c(1:3)])>0)

Rs_vias <- generatePathwaysTable(patients_dir = "~/Daniela/Biota/Muestras/73m", de_host = "RSubread")
Rs_vias_CPM <- Rs_vias[[1]]
Rs_vias_AR <- Rs_vias[[2]]

Bwa_vias <- generatePathwaysTable(patients_dir = "~/Daniela/Biota/Muestras/73m", de_host = "BWA")
Bwa_vias_CPM <- Bwa_vias[[1]]
Bwa_vias_AR <- Bwa_vias[[2]]

library(ggvenn)


data_list <- list(
  '-H' = unique(T_vias_CPM$Pathway),
  BoH = unique(Bo_vias_CPM$Pathway),
  bwaH = unique(Bwa_vias_CPM$Pathway),
  RsH = unique(Rs_vias_CPM$Pathway)
)

# Crear el diagrama de Venn
ggvenn(data_list,
       #fill_color = c("red", "blue", "pink", "green", "purple", "black"),
       fill_color = c(  "#9467BD", "#E69F00","#A52A2A", "#009E73" ),
       stroke_size = 0,
       set_name_size = 4,
       text_size = 3)

# Convertir listas a un tibble con presencia/ausencia

vias_data <- tibble(
  Vias = unique(unlist(data_list)),
  Bo = as.integer(Vias %in% data_list$BoH),
  BWA = as.integer(Vias %in% data_list$bwaH),
  Rs = as.integer(Vias %in% data_list$RsH),
  SinDH = as.integer(Vias %in% data_list$`-H`)
)

# Agrupar por combinaciones y contar especies en cada intersección
vias_summary <- vias_data %>%
  group_by(Bo, BWA, Rs, SinDH) %>%
  summarise(Vias_List = list(Vias), Count = n(), .groups = "drop")

#Las 34 que se pierden con BWA:
bwa_pierde <- vias_summary$Vias_List[which(vias_summary$BWA == 0)][[2]]
Tbwa_pierde_df <- T_vias_AR[which(T_vias_AR$Pathway %in% bwa_pierde),  ]
Tbwa_pierde_df$Min_N <- apply(Tbwa_pierde_df[,-c(1:3)], 1, min)
Tbwa_pierde_df$Max_N <- apply(Tbwa_pierde_df[,-c(1:3)], 1, max)
Tbwa_pierde_df$NSamples_N <- apply(Tbwa_pierde_df[,-c(1:3)], 1, function(x) sum(x > 0))
Tbwa_pierde_df$PercentageSamples_N <- (Tbwa_pierde_df$NSamples/83)*100
Tbwa_pierde_df$Mean_N <- apply(Tbwa_pierde_df[,-c(1:3)], 1, mean)
colnames(Tbwa_pierde_df)
str(Tbwa_pierde_df)
Tbwa_pierde <- Tbwa_pierde_df[,c("Pathway","Clase", "Description","Min_N", "Max_N", "NSamples_N","PercentageSamples_N",  "Mean_N")]

bwa_pierde <- vias_summary$Vias_List[which(vias_summary$BWA == 0)][[2]]
Bobwa_pierde_df <- Bo_vias_AR[which(Bo_vias_AR$Pathway %in% bwa_pierde),  ]
Bobwa_pierde_df$Min_Bo <- apply(Bobwa_pierde_df[,-c(1:3)], 1, min)
Bobwa_pierde_df$Max_Bo <- apply(Bobwa_pierde_df[,-c(1:3)], 1, max)
Bobwa_pierde_df$NSamples_Bo <- apply(Bobwa_pierde_df[,-c(1:3)], 1, function(x) sum(x > 0))
Bobwa_pierde_df$PercentageSamples_Bo <- (Bobwa_pierde_df$NSamples/83)*100
Bobwa_pierde_df$Mean_Bo <- apply(Bobwa_pierde_df[,-c(1:3)], 1, mean)
colnames(Bobwa_pierde_df)
str(Bobwa_pierde_df)
Bobwa_pierde_df <- Bobwa_pierde_df[,c("Pathway","Clase", "Description","Min_Bo", "Max_Bo", "NSamples_Bo","PercentageSamples_Bo",  "Mean_Bo")]

bwa_pierde <- vias_summary$Vias_List[which(vias_summary$BWA == 0)][[2]]
Rsbwa_pierde_df <- Rs_vias_AR[which(Rs_vias_AR$Pathway %in% bwa_pierde),  ]
Rsbwa_pierde_df$Min_Rs <- apply(Rsbwa_pierde_df[,-c(1:3)], 1, min)
Rsbwa_pierde_df$Max_Rs <- apply(Rsbwa_pierde_df[,-c(1:3)], 1, max)
Rsbwa_pierde_df$NSamples_Rs <- apply(Rsbwa_pierde_df[,-c(1:3)], 1, function(x) sum(x > 0))
Rsbwa_pierde_df$PercentageSamples_Rs <- (Rsbwa_pierde_df$NSamples/83)*100
Rsbwa_pierde_df$Mean_Rs <- apply(Rsbwa_pierde_df[,-c(1:3)], 1, mean)
Rsbwa_pierde_df <- Rsbwa_pierde_df[,c("Pathway","Clase", "Description","Min_Rs", "Max_Rs", "NSamples_Rs","PercentageSamples_Rs","Mean_Rs")]

Todos_pierde_bwa_df <- merge(Tbwa_pierde, Bobwa_pierde_df, by = c("Pathway", "Clase", "Description"))
Todos_pierde_bwa_df <- merge(Todos_pierde_bwa_df, Rsbwa_pierde_df, by =  c("Pathway", "Clase", "Description"))
Todos_pierde_bwa_df[,-c(1:3)] <- round(Todos_pierde_bwa_df[,-c(1:3)], 2)
Todos_pierde_bwa_df <- Todos_pierde_bwa_df[order(Todos_pierde_bwa_df$NSamples_N, decreasing = TRUE),]

write.xlsx(Todos_pierde_bwa_df, file = "~/Daniela/Biota/PipelineBiota/paraPaper/Todos_Pierde_BWA_Vias.xlsx")

todos_menossinDH_pierden <- vias_summary$Vias_List[which(vias_summary$SinDH == 1)][[1]]
todos_menossinDH_pierden_df <- T_vias_AR[which(T_vias_AR$Pathway %in% todos_menossinDH_pierden), c(1,2,3) ]

vs <- Bo_vias_AR[Bo_vias_AR$Pathway %in% c("PWY-7242"),]
vs <- Rs_vias_AR[Rs_vias_AR$Pathway %in% c("PWY-7242"),]
vs <- Bwa_vias_AR[Bwa_vias_AR$Pathway %in% c("PWY-7242"),]
vs <- T_vias_AR[T_vias_AR$Pathway %in% c("PWY-7242"),]
vs <- T_vias_AR[T_vias_AR$Pathway %in% bwa_pierde,]
colnames(vs)
vs$Cant_muestras <- rowSums(vs[,-c(1,2,3,87)] > 0.000000)
vs <- vs[order(vs$Cant_muestras, decreasing = TRUE),]
vs[1,]

vs$MeanAbundance <- apply(vs[,-c(1,2,3,87)], 1, function(row) {
  values <- row[row > 0]
  if (length(values) > 0) {
    mean(values)  # Calcular el promedio solo de esos valores
  }
})
vs_lim <- vs[,c(1,2,3,87,88)]

write.xlsx(vs_lim, file = "~/Daniela/Biota/PipelineBiota/paraPaper/Vias_BWA_pierde.xlsx")

bosi_bwano <- Bo_vias_CPM[which(Bo_vias_CPM$Pathway %in% bwa_pierde),]


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




