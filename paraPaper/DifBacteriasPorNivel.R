library(PipelineBiota)
library(openxlsx)
library(dplyr)

# Especies:
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

# Generos:
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



#Diagrama de venn: -----------------------------------------------------------------------------------
library(ggvenn)
# Convertir los datos a una lista
data_list <- list(
  BoK = unique(AR_speciesK_Bo$Species),
  bwaK = unique(AR_speciesK_BWA$Species),
  RsK = unique(AR_speciesK_Rs$Species),
  '-K' = unique(AR_speciesK_sin$Species),

  BoD = unique(AR_speciesD_Bo$Species),
  bwaD = unique(AR_speciesD_BWA$Species),
  RsD = unique(AR_speciesD_Rs$Species),
  DdhD = unique(AR_speciesD_sin$Species),
  '-D' = unique(AR_speciesD_sinDH_PD$Species)
)


data_list <- list(
  BoK = unique(AR_genK_Bo$Genus),
  bwaK = unique(AR_genK_BWA$Genus),
  RsK = unique(AR_genK_Rs$Genus),
  '-K' = unique(AR_genK_sin$Genus),

  BoD = unique(AR_genD_Bo$Genus),
  bwaD = unique(AR_genD_BWA$Genus),
  RsD = unique(AR_genD_Rs$Genus),
  DdhD = unique(AR_genD_sin$Genus),
  '-D' = unique(ARgenD_sinDH_PD$Genus)
)



data_list <- list(
  FilosK = unique(AR_filosK$Phylum),
  FilosD = unique(AR_filosD$Phylum)
)

# Crear el diagrama de Venn

library(UpSetR)
library(ComplexHeatmap)  # Necesario para fromList()
data_matrix <- fromList(data_list)

#>>>>>>>>>>>>>>>> FIGURE 3B <<<<<<<<<<<<<<<<<<
# Definir los colores: negro para la mayoría, colores solo en las primeras 4 barras
bar_colors <- rep("black", 25)
bar_colors[c(2,3,9)] <- "grey40"  # Colores para destacar

bar_colors <- rep("grey15", 25)
bar_colors[c(2,3,9)] <- c("#9370DB", "#FF7F00", "#56B4E9")  # Colores para destacar

upset_plot <- upset(
  data_matrix,
  sets = names(data_list),
  order.by = "freq",
  mainbar.y.label = "Number of Species",
  #mainbar.y.label = "Number of Genus",
  sets.x.label = "Nº Detected Species",
  text.scale = c(1.3, 1.3, 1, 1, 1.3, 1.3),  # Aumenta la legibilidad del texto
  point.size = 2,  # Puntos más grandes en la matriz
  line.size = 0,  # Líneas más gruesas para mayor visibilidad
  keep.order = TRUE,  # Mantiene el orden original de los conjuntos
  sets.bar.color = "grey50",  # Todas las barras de conjuntos en negro
  matrix.color = "black",  # Matriz en negro
  main.bar.color = bar_colors,
  shade.color = "gray90",  # Fondo de sombreado más claro para contraste
  number.angles = 0,  # Mejor legibilidad de números
  #mainbar.y.max = max(rowSums(data_matrix)) * 0.8,  # Reduce la altura de la barra principal
  #mb.ratio = c(0.6, 0.4)  # Ajusta la proporción entre matriz y barra principal
)

upset_plot

upset_gen <- upset_plot
upset_esp <- upset_plot

# Convertir a data.frame para manipulación
generos <- unique(unlist(data_list))
df_presence <- as.data.frame(data_matrix)
df_presence$Genero <- generos  # Agregar nombres de géneros
df_presence$Interseccion <- apply(df_presence[,-ncol(df_presence)], 1, paste, collapse = "-")

# Agrupar géneros por intersección
intersecciones <- split(df_presence$Genero, df_presence$Interseccion)
colnames(df_presence)
g_enK_noD <- df_presence$Genero[df_presence$Interseccion == "1-1-1-1-0-0-0-0-0"]
g_enD_noK <- df_presence$Genero[df_presence$Interseccion == "0-0-0-0-1-1-1-1-1"]
g_noBWA <- df_presence$Genero[df_presence$Interseccion == "1-0-1-1-1-0-1-1-1"]

e_enK_noD <- df_presence$Genero[df_presence$Interseccion == "1-1-1-1-0-0-0-0-0"]
e_enD_noK <- df_presence$Genero[df_presence$Interseccion == "0-0-0-0-1-1-1-1-1"]
e_noBWA <- df_presence$Genero[df_presence$Interseccion == "1-0-1-1-1-0-1-1-1"]

#Ver en cuantas muestras aparecen estas especies que se pierden: ----------------

  #las que no aparecen en Dragen:

KBo_enoBWA <- AR_speciesK_Bo[AR_speciesK_Bo$Species %in% e_enK_noD,]
KBo_enoBWA <- AR_speciesK_Rs[AR_speciesK_Rs$Species %in% e_enK_noD,]
KBo_enoBWA <- AR_speciesK_BWA[AR_speciesK_BWA$Species %in% e_enK_noD,]
KBo_enoBWA <- AR_speciesK_sin[AR_speciesK_sin$Species %in% e_enK_noD,]

KBo_enoBWA <- AR_speciesK_sin[AR_speciesK_sin$Species %in% e_enK_noD,]
KBo_enoBWA <- AR_speciesK_Rs[AR_speciesK_Rs$Species %in% e_enK_noD,]
KBo_enoBWA <- AR_speciesK_BWA[AR_speciesK_BWA$Species %in% e_enK_noD,]

 #las que no aparecen en bwa:
KBo_enoBWA <- AR_speciesK_Bo[AR_speciesK_Bo$Species %in% e_noBWA,]
KBo_enoBWA <- AR_speciesD_Bo[AR_speciesD_Bo$Species %in% e_noBWA,]
KBo_enoBWA <- AR_speciesD_sin[AR_speciesD_sin$Species %in% e_noBWA,]
KBo_enoBWA <- AR_speciesK_sin[AR_speciesK_sin$Species %in% e_noBWA,]
KBo_enoBWA <- AR_speciesK_Rs[AR_speciesK_Rs$Species %in% e_noBWA,]
KBo_enoBWA <- AR_speciesD_Rs[AR_speciesD_Rs$Species %in% e_noBWA,]

#Veo en cuantas muestras aparecen (rango):
colnames(KBo_enoBWA)
KBo_enoBWA$CountSamples <- rowSums(KBo_enoBWA[,-c(1,85)] > 0)
summary(KBo_enoBWA$CountSamples)

# Veo cual es la especie que aparece en mayor cantidad de muestras:
KBo_enoBWA <- KBo_enoBWA[order(KBo_enoBWA$CountSamples, decreasing = TRUE),]
KBo_enoBWA[1,] #"Sulfurihydrogenibium sp."
#cuando es con DRAGEN es "Aquifex aeolicus"

str(KBo_enoBWA)
KBo_enoBWA$MeanAbundance <- apply(KBo_enoBWA[,-c(1,85, 86)], 1, function(row) {
  values <- row[row > 0]
  if (length(values) > 0) {
    mean(values)  # Calcular el promedio solo de esos valores
  }
})


str(AR_speciesK_Bo)
rowMeans(AR_speciesK_Bo[AR_speciesK_Bo$Species == "Aquifex aeolicus", -c(1,85)])
rowMeans(AR_speciesD_Bo[AR_speciesD_Bo$Species == "Aquifex aeolicus", -c(1,85)])
rowSums(AR_speciesK_Bo[AR_speciesK_Bo$Species == "Aquifex aeolicus",-c(1,85)] > 0)
rowSums(AR_speciesD_Bo[AR_speciesD_Bo$Species == "Aquifex aeolicus",-c(1,85)] > 0)

#---------------------------------------------------------------------

# Encuentra la longitud máxima de los vectores
max_length <- max(length(g_enK_noD), length(g_enD_noK), length(g_noBWA))
max_length <- max(length(e_enK_noD), length(e_enD_noK), length(e_noBWA))

# Rellena cada vector con NA hasta la misma longitud
parabuscar <- data.frame(
  g_enK_noD = c(g_enK_noD, rep(NA, max_length - length(g_enK_noD))),
  g_enD_noK = c(g_enD_noK, rep(NA, max_length - length(g_enD_noK))),
  g_noBWA   = c(g_noBWA, rep(NA, max_length - length(g_noBWA)))
)

parabuscar <- data.frame(
  e_enK_noD = c(e_enK_noD, rep(NA, max_length - length(e_enK_noD))),
  e_enD_noK = c(e_enD_noK, rep(NA, max_length - length(e_enD_noK))),
  e_noBWA   = c(e_noBWA, rep(NA, max_length - length(e_noBWA)))
)

write.xlsx(parabuscar, file = "~/Daniela/Biota/PipelineBiota/paraPaper/generos_parabuscar.xlsx")
write.xlsx(parabuscar, file = "~/Daniela/Biota/PipelineBiota/paraPaper/especies_parabuscar.xlsx")

names(data_list)
ggvenn(data_list[c(1,2,5,6,9)],
       #fill_color = c("red", "blue", "green", "black"),
       #fill_color = c("red", "blue", "green", "black", "grey"),
       fill_color = c("darkred", "darkgreen", "red",  "green", "black"),
       stroke_size = 0.3,
       set_name_size = 3.5,
       text_size = 3.5)

library(VennDiagram)
length(data_list)
venn.plot <- venn.diagram(
  x = data_list,
  filename = NULL,
  fill = rainbow(length(data_list)),
  alpha = 0.5,  # Transparencia
  cex = 1.5,  # Tamaño del texto
  cat.cex = 1.2,  # Tamaño de nombres de conjuntos
  lwd = 1  # Grosor de las líneas
)

grid::grid.draw(venn.plot)




library(VennDiagram)
venn.plot <- venn.diagram(
  x = list(SpeciesK = unique(AR_speciesK$Species), SpeciesD = unique(AR_speciesD$Species)),
  category.names = c("Species KRAKEN", "Species DRAGEN"),
  filename = NULL, # Para no guardar en un archivo, puedes cambiarlo si lo deseas
  output = TRUE,
  #fill = c("lightpink", "lightgreen"), # Colores de las áreas
  fill = c(alpha("#440154ff", 0.3), alpha('#21908dff', 0.3)),
  col = c("#440154ff", '#21908dff'),
  alpha = 0.7, # Transparencia
  cex = 1.2, # Tamaño de texto
  fontface = "bold", # Estilo de fuente
  fontfamily = "serif",
  cat.cex = 1.2, # Tamaño del texto de las categorías
  cat.fontface = "bold",
  lwd = 1.5, # Ancho de las líneas
  euler.d = TRUE, # Ajusta el diagrama a un estilo más compacto
  scaled = TRUE, # Escala el diagrama
  cat.pos = c(-20, 20), # Posición de los textos de las categorías
  cat.dist = c(0.05, 0.05), # Distancia de las categorías
  cat.col = c("#440154ff", '#21908dff')
)
# Mostrar el diagrama
grid.draw(venn.plot)

species_data <- tibble(
  Species = unique(unlist(data_list)),
  Kraken = as.integer(Species %in% data_list$SpeciesK),
  Dragen = as.integer(Species %in% data_list$SpeciesD),
)
species_data <- tibble(
  Generos = unique(unlist(data_list)),
  Kraken = as.integer(Generos %in% data_list$GenerosK),
  Dragen = as.integer(Generos %in% data_list$GenerosD),
)
species_data <- tibble(
  Filos = unique(unlist(data_list)),
  Kraken = as.integer(Filos %in% data_list$FilosK),
  Dragen = as.integer(Filos %in% data_list$FilosD),
)

# Agrupar por combinaciones y contar especies en cada intersección
species_summary <- species_data %>%
  group_by(Kraken, Dragen) %>%
  summarise(Species_List = list(Species), Count = n(), .groups = "drop")
#summarise(Filos_List = list(Filos), Count = n(), .groups = "drop")
#summarise(Generos_List = list(Generos), Count = n(), .groups = "drop")
