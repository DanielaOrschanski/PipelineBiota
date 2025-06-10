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

# >>>>>>>>>>>> FIGURE S1 and 5b <<<<<<<<<<<<

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

dfs <- lista_p_gen
nombres <- names(dfs)
# Función para renombrar columnas
library(dplyr)
library(purrr)

colnames(dfs[[1]])
renombrar_columnas <- function(df, nombre) {
  df %>%
    rename_with(~ paste0(.x, "_", nombre), c("Wilcoxon", "GrupoDominante"))
  }

# Renombrar y unir por "Categoria" y "Via"
str(dfs)
dfs_renombrados <- map2(dfs, nombres, renombrar_columnas)

df_final <- dfs_renombrados[[1]]

for (i in 2:length(dfs_renombrados)) {
  df_final <- merge(df_final, dfs_renombrados[[i]], by = c("Categoria", "Via"), all = TRUE)
}

df_final <- na.omit(df_final)

df_significativo <- df_final %>%
  filter(
    rowSums(select(., starts_with("Wilcoxon_")) < 0.05, na.rm = TRUE) > 0 &
      rowSums(select(., starts_with("Wilcoxon_")) > 0.05, na.rm = TRUE) > 0
  )

df_sign_rango <- df_significativo[df_significativo$Categoria == "Rango etario",]
df_sign_sexo <- df_significativo[df_significativo$Categoria == "Sexo",]

file = "~/Daniela/Biota/PipelineBiota/paraPaper/DifSign_Genero_RangoEtario.xlsx"

#Fig 5b ---------------------------------------------------------------
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

bar_colors <- rep("grey40", 45)
bar_colors[c(9, 10, 16)] <- "black"
bar_colors[c(1,3)] <- "grey15"


upset(
  data_matrix,
  sets = names(data_list),
  order.by = "freq",
  mainbar.y.label = "Nº Genera signif dif between sexes",
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



# Fig S1: -------------------------------------------------------------------
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

library(UpSetR)
library(ComplexHeatmap)
data_matrix <- fromList(data_list)

bar_colors <- rep("grey40", 47)
bar_colors[c(5,11,14)] <- "black"  # Colores para destacar
bar_colors[c(1,2)] <- "grey15"  # Colores para destacar

upset(
  data_matrix,
  sets = names(data_list),
  order.by = "freq",
  mainbar.y.label = "Nº Genera signif dif between age ranges",
  sets.x.label = "Nº Signif Genera Detected",
  text.scale = c(1.5, 1.3, 1.2, 1, 1.3, 1.3),
  point.size = 2,
  line.size = 0,
  keep.order = TRUE,  # Mantiene el orden original de los conjuntos
  sets.bar.color = "grey50",  # Todas las barras de conjuntos en negro
  matrix.color = "black",  # Matriz en negro
  main.bar.color = bar_colors,
  shade.color = "gray90",  # Fondo de sombreado más claro para contraste
  number.angles = 0,  # Mejor legibilidad de números
  mb.ratio = c(0.7, 0.3)
)


# >>>>>>>>>>>>>>>>>>> FIGURE 5C <<<<<<<<<<<<<<<<<<<<<<<

df_long <- df_gen_sign_sexo %>%
  pivot_longer(cols = -Via,
               names_to = c(".value", "Method"),
               names_sep = "_")


df_long$pValor <- ifelse(df_long$Wilcoxon<0.05, -log10(df_long$Wilcoxon), NA)

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
    breaks = c(-log10(0.049), -log10(0.01), -log10(0.0015)),
    labels = c("0.05", "0.01", "0.001"),
    range = c(2, 12)
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
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13)
  ) + guides(
    fill = guide_legend(
      override.aes = list(size = 5)
    )
  )


