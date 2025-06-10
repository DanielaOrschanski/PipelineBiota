patients_dir <- "~/Daniela/Biota/Muestras/73m"
de_host = "Bowtie"

QCcontrol(patients_dir = "~/Daniela/Biota/Muestras/73m", de_host ="BWA", generate_QCReport_Individual = FALSE, generate_QCReport_Grupal = FALSE, FastQC_trimmed = TRUE)
QCcontrol(patients_dir = "~/Daniela/Biota/Muestras/73m", de_host ="RSubread", generate_QCReport_Individual = FALSE, generate_QCReport_Grupal = FALSE, FastQC_trimmed = TRUE)
QCcontrol(patients_dir = "~/Daniela/Biota/Muestras/73m", de_host ="", generate_QCReport_Individual = FALSE, generate_QCReport_Grupal = FALSE, FastQC_trimmed = TRUE)

GC_Bo <- plot_GCcontent(patients_dir= patients_dir, de_host = "Bowtie")
Tabla_Bo <- GC_Bo[[2]]

GC_BWA <- plot_GCcontent(patients_dir= patients_dir, de_host = "BWA")
Tabla_BWA <- GC_BWA[[2]]

GC_Rs <- plot_GCcontent(patients_dir= patients_dir, de_host = "RSubread")
Tabla_Rs <- GC_Rs[[2]]

GC_sin <- plot_GCcontent(patients_dir= patients_dir, de_host = "")
Tabla_sin <- GC_sin[[2]]

# Unimos todas las tablas
Tabla_todas <- bind_rows(Tabla_Bo, Tabla_BWA, Tabla_Rs, Tabla_sin)
str(Tabla_todas)

library(dplyr)

# Ordenar los datos por Method, ID y GC content
Tabla_todas_ordenada <- Tabla_todas %>%
  arrange(Method, ID, `GC Content`)
unique(Tabla_todas_ordenada$Method)
Tabla_todas_ordenada$Method[Tabla_todas_ordenada$Method == ""] <- "No-dehost"
Tabla_todas_ordenada$Method <- factor(Tabla_todas_ordenada$Method, levels = c("No-dehost", "Bowtie", "BWA", "RSubread"))

write.xlsx(Tabla_todas_ordenada, file = "~/Daniela/Biota/PipelineBiota/paraPaper/Scripts Reproducir Paper/Fig2b-GC.xlsx")

GC_content_AllMethods <- read_excel("Daniela/Biota/PipelineBiota/paraPaper/Scripts Reproducir Paper/Fig2b-GC.xlsx")
GC_content_AllMethods$Method <- factor(GC_content_AllMethods$Method, levels = c("No-dehost", "Bowtie", "BWA", "RSubread"))

# Crear gráfico: group por ID + Method
GC_content_AllMethods <- Tabla_todas_ordenada
Tabla_todas_ordenada <- GC_content_AllMethods
ggplot(GC_content_AllMethods, aes(x = `GC Content`, y = Count, group = interaction(ID, Method), color = Method)) +
  #geom_point(alpha = 0.4, size = 1) +
  geom_line(alpha = 0.7) +
  facet_wrap(~ Method, ncol = 4) +
  labs(
    title = "",
    x = "GC content (%)",
    y = "Number of sequences"
  ) +
  scale_color_manual(values = c(
    "Bowtie" = "#E69F00",
    "BWA" = "#56B4E9",
    "RSubread" = "#009E73",
    "No-dehost" = "#EE6363"
  )) +
  scale_x_continuous(breaks = seq(0, 100, by = 20)) +
  scale_y_continuous(breaks = seq(0, 200000, by = 25000)) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "none",
    strip.text = element_text( size = 15)
  )

max(GC_content_AllMethods$Count)
#Comparacion entre picos alrededor de 40%:
library(dplyr)

# Filtrar las regiones de GC
Tabla_40 <- Tabla_todas_ordenada %>%
  filter(`GC Content` >= 30, `GC Content` <= 50)

Tabla_60 <- Tabla_todas_ordenada %>%
  filter(`GC Content` >= 50, `GC Content` <= 70)

maximos40 <- Tabla_40 %>%
  group_by(ID, Method) %>%
  summarise(MaxCount = max(Count), .groups = "drop") %>%
  mutate(GC_peak = "GC ~40%")

maximos60 <- Tabla_60 %>%
  group_by(ID, Method) %>%
  summarise(MaxCount = max(Count), .groups = "drop") %>%
  mutate(GC_peak = "GC ~60%")

maximos <- bind_rows(maximos40, maximos60)

maximos$Method <- factor(maximos$Method, levels = c("No-dehost", "Bowtie", "BWA", "RSubread"))
maximos$GC_peak <- factor(maximos$GC_peak, levels = c("GC ~40%", "GC ~60%"))

unique(maximos$Method)
str(maximos)
maximos <- as.data.frame(maximos)
maximos$ID <- as.factor(maximos$ID)
#Mas eficientes:

wilcox_results_peak <- maximos %>%
  filter(Method %in% c("BWA", "RSubread", "Bowtie")) %>%
  group_by(GC_peak) %>%
  nest() %>%
  mutate(
    data_wide = map(data, ~ pivot_wider(.x, names_from = Method, values_from = MaxCount)),

    # Wilcoxon BWA vs RSubread
    test_BWA_Rs = map(data_wide, ~ wilcox.test(.x$BWA, .x$RSubread, paired = TRUE)),
    p_BWA_Rs = map_dbl(test_BWA_Rs, ~ .x$p.value),

    # Wilcoxon BWA vs Bowtie
    test_BWA_Bo = map(data_wide, ~ wilcox.test(.x$BWA, .x$Bowtie, paired = TRUE)),
    p_BWA_Bo = map_dbl(test_BWA_Bo, ~ .x$p.value),

    # Wilcoxon RSubread vs Bowtie
    test_Rs_Bo = map(data_wide, ~ wilcox.test(.x$RSubread, .x$Bowtie, paired = TRUE)),
    p_Rs_Bo = map_dbl(test_Rs_Bo, ~ .x$p.value),

    # Diferencias porcentuales promedio
    pct_diff_BWA_Rs = map_dbl(data_wide, ~ mean((.x$BWA - .x$RSubread) / .x$BWA * 100, na.rm = TRUE)),
    pct_diff_BWA_Bo = map_dbl(data_wide, ~ mean((.x$BWA - .x$Bowtie) / .x$BWA * 100, na.rm = TRUE)),
    pct_diff_Bo_Rs  = map_dbl(data_wide, ~ mean((.x$Bowtie - .x$RSubread) / .x$Bowtie * 100, na.rm = TRUE))
)

wide <- wilcox_results_peak$data_wide[[1]]
wide_60 <- wilcox_results_peak$data_wide[[which(wilcox_results_peak$GC_peak == "GC ~60%")]]
wide_40 <- wilcox_results_peak$data_wide[[which(wilcox_results_peak$GC_peak == "GC ~40%")]]

wilcox.test(wide_60$BWA, wide_60$RSubread, paired = TRUE)
wilcox.test(wide_60$Bowtie, wide_60$RSubread, paired = TRUE)

wide_60$Bo_BWA <- (wide_60$Bowtie - wide_60$BWA) / wide_60$Bowtie *100
wide_40$Bo_BWA <- (wide_40$Bowtie - wide_40$BWA) / wide_40$Bowtie *100

#------


maximos_wide <- maximos %>%
  pivot_wider(names_from = Method, values_from = MaxCount) %>%
  mutate(Diferencia_RS_Bo = ((RSubread - Bowtie)/Bowtie)*100) %>%
  mutate(Diferencia_N_Rs = ((`No-dehost` - RSubread)/RSubread)*100) %>%
  mutate(Diferencia_N_Bo = ((`No-dehost` - Bowtie)/Bowtie)*100) %>%
  mutate(Diferencia_Bo_BWA = ((Bowtie - BWA)/BWA)*100) %>%
  mutate(Diferencia_Rs_BWA = ((RSubread - BWA)/BWA)*100)

maximos_wide$`No-dehost`

maximos_long <- maximos_wide %>%
  select(ID, GC_peak, Diferencia_RS_Bo, Diferencia_N_Rs, Diferencia_N_Bo, Diferencia_Bo_BWA, Diferencia_Rs_BWA) %>%
  pivot_longer(
    cols = starts_with("Diferencia"),
    names_to = "Comparacion",
    values_to = "Porcentaje"
  )


#Boxplot de diferencias entre metodos:
ggplot(maximos_long, aes(x = GC_peak, y = Porcentaje, fill = Comparacion)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_violin(alpha = 0.6, color = NA) +
  geom_boxplot(width = 0.4, outlier.shape = NA, alpha = 0.5, color = "black") +
  #geom_jitter(width = 0.1, size = 1, color = "black") +
  facet_wrap(~ Comparacion, scales = "fixed", nrow=1) +
  scale_fill_manual(values = c(
    "Diferencia_RS_Bo" = "#009E73",
    "Diferencia_N_Rs" = "#E69F00",
    "Diferencia_N_Bo" = "#56B4E9"
  )) +
  theme_minimal() +
  labs(
    title = "Porcentaje de diferencia entre métodos en picos de GC",
    x = "Pico de GC",
    y = "Diferencia porcentual (%)"
  ) +
  theme(
    text = element_text(size = 12),
    strip.text = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 12, face = "bold"),
    legend.position = "none"
  )


wilcox.test(maximos_wide$`No-dehost`, maximos_wide$Bowtie, paired = TRUE)
wilcox.test(maximos_wide$`No-dehost`, maximos_wide$BWA, paired = TRUE)
wilcox.test(maximos_wide$`No-dehost`, maximos_wide$RSubread, paired = TRUE)

wilcox.test(maximos_wide$BWA, maximos_wide$Bowtie, paired = TRUE)
wilcox.test(maximos_wide$BWA, maximos_wide$RSubread, paired = TRUE)
wilcox.test(maximos_wide$Bowtie, maximos_wide$RSubread, paired = TRUE)

wilcox.test(maximos_wide$`No-dehost`, maximos_wide$Bowtie, paired = FALSE)
wilcox.test(maximos_wide$`No-dehost`, maximos_wide$BWA, paired = FALSE)
wilcox.test(maximos_wide$`No-dehost`, maximos_wide$RSubread, paired = FALSE)

wilcox.test(maximos_wide$BWA, maximos_wide$Bowtie, paired = FALSE)
wilcox.test(maximos_wide$BWA, maximos_wide$RSubread, paired = FALSE)
wilcox.test(maximos_wide$Bowtie, maximos_wide$RSubread, paired = FALSE)


library(ggplot2)

library(ggplot2)
library(dplyr)

# Gráfico final
ggplot(maximos[maximos$Method != "No-dehost", ], aes(x = Method, y = MaxCount, fill = Method, color = Method)) +
  geom_violin(alpha = 0.6, width = 1, color = NA) +
  geom_jitter(position = position_jitter(width = 0.15), size = 2, alpha = 0.7) +
  geom_line(aes(group=ID), colour="black", linetype="11", alpha = 0.8) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.5, color = "black") +

  facet_wrap(~ GC_peak, scales = "fixed") +
  scale_fill_manual(values = c(
    "Bowtie" = "#E69F00",
    "BWA" = "#56B4E9",
    "RSubread" = "#009E73"
  )) +
  scale_color_manual(values = c(
    "Bowtie" = "#E69F00",
    "BWA" = "#56B4E9",
    "RSubread" = "#009E73"
  )) +
  theme_minimal() +
  labs(
    title = "Peak sequence count at GC content ~40% and ~60%",
    x = "Method",
    y = "Maximum sequence count"
  ) +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(face = "bold", size = 14),
    legend.position = "none"
  )

# Subset para GC ~40%
df_40 <- maximos40
p_values <- pairwise.wilcox.test(
  x = df_40$MaxCount,
  g = df_40$Method,
  paired = FALSE, # o FALSE, según tu diseño
  p.adjust.method = "BH"
)
print(p_values)

# Subset para GC ~60%
df_60 <- maximos60
p_values <- pairwise.wilcox.test(
  x = df_60$MaxCount,
  g = df_60$Method,
  paired = FALSE, # o FALSE, según tu diseño
  p.adjust.method = "BH"
)
print(p_values)

plot_GCcontent <- function(patients_dir, de_host) {

  if(de_host  == "Bowtie") {
    de_host_file <- "DHBo"
  } else if( de_host == "BWA") {
    de_host_file <- "DHbwa"
  } else if(de_host == "RSubread") {
    de_host_file <- "DHRs"
  } else if(de_host == "") {
    de_host_file <- "T"
  }  else {
    stop("de_host must be Bowtie, BWA, RSubread or empty string")
  }

  #Separación de los modulos del fastqc para generar graficos
  dir_list <- list.dirs(path = patients_dir, full.names = TRUE, recursive = FALSE)
  cant_patients <- length(dir_list)

  list_modulo_GC <- list() #Va a almacenar los dataframes para el grafico Per Base Sequence Quality
  Table_GC <- data.frame()
  p=1
  #Recorro cada paciente para guardar las tablas para el grafico y sacar el promedio de las medias.
  for (p in 1:cant_patients) {

    id <- basename(dir_list[p])
    print(id)
    txt_R1 <- sprintf("%s/trimmed/%s%s_S04_L001_R1_001_fastqc/fastqc_data.txt", dir_list[p], id, de_host_file)
    txt_R2 <- sprintf("%s/trimmed/%s%s_S04_L001_R2_001_fastqc/fastqc_data.txt", dir_list[p], id, de_host_file)

    if(!file.exists(txt_R1)) {
      fastqc_R1_zip <- sprintf("%s/trimmed/%s%s_S04_L001_R1_001_fastqc.zip", dir_list[p], id, de_host_file)
      if(file.exists(fastqc_R1_zip)){
        unzip(fastqc_R1_zip, exdir = sprintf("%s/trimmed", dir_list[p]))
      }

    }
    if( !file.exists(txt_R2)) {
      fastqc_R2_zip <- sprintf("%s/trimmed/%s%s_S04_L001_R2_001_fastqc.zip", dir_list[p], id, de_host_file)
      unzip(fastqc_R2_zip, exdir = sprintf("%s/trimmed", dir_list[p]))
    }

    library(readr)
    report_R1 <- read_file(txt_R1)
    library(stringr)
    module_R1 <- str_split(report_R1, ">>")

    report_R2 <- read_file(txt_R2)
    module_R2 <- str_split(report_R2, ">>")

    #Basic stats:
    GCcontent_R1 <- create_FQCdata(module_R1[[1]][12])
    GCcontent_R2 <- create_FQCdata(module_R2[[1]][12])
    GCcontent_R1$ID <- id
    Table_GC <- rbind(Table_GC, GCcontent_R1)

    #Plot principal: -----------------------------------------------------------------
    str(GCcontent_R1)
    list_modulo_GC[[p]] <- GCcontent_R1
    message(sprintf("The module of patient%s for GC content graph has been loaded", id))

  }


  # Crear el gráfico
  curva_GC <- ggplot(Table_GC, aes(x = `GC Content`, y = Count, color = as.factor(ID))) +
    geom_point() +  # Agrega puntos
    geom_line() +   # Agrega la curva
    labs(x = "GC Content", y = "Count", color = "ID") +
    theme_minimal()

  Table_GC$Method <- de_host
  return(list(curva_GC, Table_GC))

}

#' @title create FastQC data
#' @description Prepare the data for plotting
#' @param dat data from module of FastQC output
create_FQCdata <- function(dat) {
  DBname <- read_lines(dat)
  len <- length(DBname)
  DBname <- DBname[2:len]
  DBname <- read.table(text = DBname, sep = "\t")
  colnames(DBname) <- strsplit(read_lines(dat)[2],"\t", fixed=TRUE)[[1]]
  colnames(DBname)[1] <- str_replace(colnames(DBname)[1], "#", "")
  return(DBname)
}

