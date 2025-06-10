######################################################################
#   cantidad de lecturas finales ############################################################
############################################################################################
patients_dir <- "~/Daniela/Biota/Muestras/73m"

#cant_lecturasT <- calculate_cant_lecturas(patients_dir = patients_dir, de_host = "")
cant_lecturasT <- as.data.frame(read_excel("~/Daniela/Biota/PipelineBiota/paraPaper/cant_lecturas_T.xlsx"))
cant_lecturasBo <- as.data.frame(read_excel("~/Daniela/Biota/PipelineBiota/paraPaper/cant_lecturas_DHBo.xlsx"))
#cant_lecturasBo <- calculate_cant_lecturas(patients_dir = patients_dir, de_host = "Bowtie")
cant_lecturasBWA <- as.data.frame(read_excel("~/Daniela/Biota/PipelineBiota/paraPaper/cant_lecturas_DHbwa.xlsx"))
#cant_lecturasBWA <- calculate_cant_lecturas(patients_dir = patients_dir, de_host = "BWA")
cant_lecturasRs <- as.data.frame(read_excel("~/Daniela/Biota/PipelineBiota/paraPaper/cant_lecturas_DHRs.xlsx"))
#cant_lecturasRs <- calculate_cant_lecturas(patients_dir = patients_dir, de_host = "RSubread")

#cant_lecturasDRAGEN <- calculate_cant_lecturas(patients_dir = patients_dir, de_host = "DRAGEN")
cant_lecturasDRAGEN <- as.data.frame(read_excel("~/Daniela/Biota/PipelineBiota/paraPaper/cant_lecturas_DRAGEN.xlsx"))


cant_lecturas_totales <- rbind(cant_lecturasT, cant_lecturasBo, cant_lecturasBWA, cant_lecturasRs, cant_lecturasDRAGEN)
cant_lecturas_totales$Method[cant_lecturas_totales$Method == "DHBo"] <- "Bowtie"
cant_lecturas_totales$Method[cant_lecturas_totales$Method == "DHRs"] <- "RSubread"
cant_lecturas_totales$Method[cant_lecturas_totales$Method == "DHbwa"] <- "BWA"
cant_lecturas_totales$Method[cant_lecturas_totales$Method == "T"] <- "Sin-DeHost"

ggplot(cant_lecturas_totales, aes(x = Method, y = R1, fill = Method)) +
  geom_boxplot() +
  geom_line(aes(group=Sample), colour="black", linetype="11") +
  theme_minimal() +
  labs(x = "", y = "Total Reads", title = "Number of reads") +
  #scale_fill_manual(values = c("Bowtie" = "red", "BWA" = "green", "RSubread" = "blue", "Sin-DeHost" = "black")) +  # Colores personalizados
  scale_fill_manual(values = c(
    "Bowtie" = scales::alpha("red", 0.5),
    "BWA" = scales::alpha("green", 0.5),
    "RSubread" = scales::alpha("blue", 0.5),
    "Sin-DeHost" = scales::alpha("black", 0.5),
    "DRAGEN" = scales::alpha("orange", 0.5)
  ))  +
  theme(
    text = element_text(size = 15),         # Cambia el tamaño de todo el texto
    axis.text = element_text(size = 12),    # Cambia el tamaño del texto en los ejes
    axis.title = element_text(size = 14),   # Cambia el tamaño de los títulos de los ejes
    plot.title = element_text(size = 16),   # Cambia el tamaño del título del gráfico
    legend.text = element_text(size = 12),  # Cambia el tamaño del texto de la leyenda
    legend.title = element_text(size = 14),  # Cambia el tamaño del título de la leyenda
    legend.position = "none"  # Ocultar la leyenda
  ) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))   # Rotar etiquetas del eje x

#solo dehosters:
unique(cant_lecturas_totales$Method)
ggplot(cant_lecturas_totales[-which(cant_lecturas_totales$Method == "Sin-DeHost"),], aes(x = Method, y = R1, fill = Method)) +
  geom_boxplot() +
  geom_line(aes(group=Sample), colour="black", linetype="11") +
  theme_minimal() +
  labs(x = "", y = "Total Reads", title = "Number of reads") +
  #scale_fill_manual(values = c("Bowtie" = "red", "BWA" = "green", "RSubread" = "blue", "Sin-DeHost" = "black")) +  # Colores personalizados
  scale_fill_manual(values = c(
    "Bowtie" = scales::alpha("red", 0.5),
    "BWA" = scales::alpha("green", 0.5),
    "RSubread" = scales::alpha("blue", 0.5),
    "DRAGEN" = scales::alpha("orange", 0.5)
  ))  +
  theme(
    text = element_text(size = 12),         # Cambia el tamaño de todo el texto
    axis.text = element_text(size = 10),    # Cambia el tamaño del texto en los ejes
    axis.title = element_text(size = 12),   # Cambia el tamaño de los títulos de los ejes
    plot.title = element_text(size = 14),   # Cambia el tamaño del título del gráfico
    legend.text = element_text(size = 10),  # Cambia el tamaño del texto de la leyenda
    legend.title = element_text(size = 12),  # Cambia el tamaño del título de la leyenda
    legend.position = "none"  # Ocultar la leyenda
  ) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))   # Rotar etiquetas del eje x


#P-valores:
library(dplyr)
cant_lecturas_totales <- cant_lecturas_totales[,-which(colnames(cant_lecturas_totales) == "R2")]
data_Bo_BWA <- cant_lecturas_totales %>% filter(Method %in% c("Bowtie", "BWA"))
data_Bo_Rs  <- cant_lecturas_totales %>% filter(Method %in% c("Bowtie", "RSubread"))
data_BWA_Rs  <- cant_lecturas_totales %>% filter(Method %in% c("BWA", "RSubread"))
data_BWA_D  <- cant_lecturas_totales %>% filter(Method %in% c("BWA", "DRAGEN"))
data_Bo_D  <- cant_lecturas_totales %>% filter(Method %in% c("Bowtie", "DRAGEN"))
data_D_Rs  <- cant_lecturas_totales %>% filter(Method %in% c("DRAGEN", "RSubread"))

# Asegurar que los datos están en el mismo orden por muestra
data_Bo_BWA <- data_Bo_BWA %>% arrange(Sample, Method) %>% pivot_wider(names_from = Method, values_from = R1)
data_Bo_Rs  <- data_Bo_Rs %>% arrange(Sample, Method) %>% pivot_wider(names_from = Method, values_from = R1)
data_BWA_Rs  <- data_BWA_Rs %>% arrange(Sample, Method) %>% pivot_wider(names_from = Method, values_from = R1)
data_BWA_D  <- data_BWA_D %>% arrange(Sample, Method) %>% pivot_wider(names_from = Method, values_from = R1)
data_Bo_D  <- data_Bo_D %>% arrange(Sample, Method) %>% pivot_wider(names_from = Method, values_from = R1)
data_D_Rs  <- data_D_Rs %>% arrange(Sample, Method) %>% pivot_wider(names_from = Method, values_from = R1)

wilcox.test(data_Bo_BWA$Bowtie, data_Bo_BWA$BWA, paired = FALSE)
wilcox.test(data_Bo_Rs$Bowtie, data_Bo_Rs$RSubread, paired = FALSE)
wilcox.test(data_BWA_Rs$BWA, data_BWA_Rs$RSubread, paired = FALSE)
wilcox.test(data_BWA_D$BWA, data_BWA_D$DRAGEN, paired = FALSE)
wilcox.test(data_D_Rs$RSubread, data_D_Rs$DRAGEN, paired = FALSE)
wilcox.test(data_Bo_D$DRAGEN, data_Bo_D$Bowtie, paired = FALSE)


calculate_cant_lecturas <- function(patients_dir, de_host) {
  cant_lecturas <- data.frame("Sample"= c(), "R1" = c(), "R2"= c(), "Method" = c())

  if(de_host  == "Bowtie") {
    de_host_file <- "DHBo"
  } else if( de_host == "BWA") {
    de_host_file <- "DHbwa"
  } else if(de_host == "RSubread") {
    de_host_file <- "DHRs"
  } else if(de_host == "") {
    de_host_file <- "T"
  } else if(de_host == "DRAGEN") {
    de_host_file <- "DRAGEN"
  }

  i=1
  dir_list <- list.dirs(path = patients_dir, full.names = TRUE, recursive = FALSE)
  cant_patients <- length(dir_list)

  for (p in dir_list) {
    p <- dir_list[[i]]
    print(p)

    if(de_host == "DRAGEN") {
      id <- basename(p)
      patient_dir <- paste(p, "/trimmed", sep="")
      nombre_p <- sprintf("%sT", id)

      folder_report <- list.dirs( sprintf("%s/DRAGEN_Reports", patient_dir), full.names = TRUE, recursive = FALSE)
      metrics_report <- sprintf("%s/%s.mapping_metrics.csv", folder_report, nombre_p)
      metrics_report <- read_csv(metrics_report)

      if(nrow(metrics_report) == 0) { #Tengo que descargar los reportes

        bs_path = path.expand("~/Daniela/Biota/PipelineBiota/data/bs")
        project_id = "436799364"
        carpeta_bs = "MuestrasTrimmed"
        command <- sprintf("%s list appsession --project-id %s", bs_path, project_id)

        output <- capture.output(system(command, intern = TRUE))
        df_output <- as.data.frame(do.call(rbind, strsplit(output, " +", perl = TRUE)))

        df_output <- df_output[which(df_output$V3 == "DRAGENHG38_sin" | df_output$V4 == "DRAGENHG38_sin"),]
        df_output <- df_output[, c(3,4,5,6,7)]

        appsession_id <- unique(df_output[which(df_output[,2] == id), 4])
        if(length(nchar(appsession_id)) == 0) {
          appsession_id <- unique(df_output[which(df_output[, 3] == id), 5])
        }

        if(length(nchar(appsession_id)) == 0) {
          df_output <- df_output[, c(1,2,4)]
          colnames(df_output) <- c("DRAGENHG38_sin", "SAMPLE", "ID")
          appsession_id <- df_output$ID[which(df_output$SAMPLE == id)]
        }

        if(length(nchar(appsession_id)) == 0) {
          output <- system(command, intern = TRUE)
          output_text <- paste(output, collapse = "\n")
          output_lines <- unlist(strsplit(output_text, "\n"))
          data_lines <- output_lines[-which(grepl("-", output_lines))]
          data_lines <- data.frame(data_lines)
          data_lines <- as.data.frame(lapply(data_lines, function(x) gsub("^\\||\\|$", "", x)))
          data_lines <- data_lines %>%
            separate(data_lines, into = c("BioSampleName", "Id", "ContainerName", "ContainerPosition", "Status"), sep = "\\|")
          appsession_id <- data_lines$Id[which(grepl(id, data_lines$BioSampleName))]
          appsession_id <- gsub(" ", "", appsession_id)

        }

        dir.create(sprintf("%s/DRAGEN_Reports", patient_dir))
        command <- sprintf("%s download appsession -i %s -o %s",
                           bs_path,
                           appsession_id,
                           sprintf("%s/DRAGEN_Reports", patient_dir))
        system(command)

      }

      folder_report <- list.dirs( sprintf("%s/DRAGEN_Reports", patient_dir), full.names = TRUE, recursive = FALSE)
      metrics_report <- sprintf("%s/%s.mapping_metrics.csv", folder_report, nombre_p)
      metrics_report <- read_csv(metrics_report)
      metrics_report <- rbind(colnames(metrics_report), metrics_report)
      colnames(metrics_report) <- c("Mapping", "Sample", "Metric", "Value", "Percentage")
      porcentage_mapeo <- as.numeric(metrics_report$Percentage[metrics_report$Metric == "Mapped reads" & grepl(id, metrics_report$Sample)])
      num_readsR1 <- round(as.numeric(metrics_report$Value[metrics_report$Metric == "Unmapped reads" & grepl(id, metrics_report$Sample)])/2)
      num_readsR2 <- num_readsR1

    } else { #Para bowtie, Rs y bwa:

      file_list <- list.files(sprintf("%s/trimmed/", p))
      gzip <- ifelse(length(nchar(file_list[endsWith(file_list, sprintf("%s_S04_L001_R1_001.fastq.gz", de_host_file))])) == 0, "", ".gz")
      fileR1 <- paste0(p, "/trimmed/", file_list[endsWith(file_list, sprintf("%s_S04_L001_R1_001.fastq%s", de_host_file, gzip))], sep="")
      fileR2 <- paste0(p, "/trimmed/", file_list[endsWith(file_list, sprintf("%s_S04_L001_R2_001.fastq%s", de_host_file, gzip))], sep="")

      library(ShortRead)
      fastq_dataR1 <- readFastq(fileR1)
      num_readsR1 <- length(fastq_dataR1)

      fastq_dataR2 <- readFastq(fileR2)
      num_readsR2 <- length(fastq_dataR2)
    }

    cant_lecturas[i, "Sample"] <- basename(p)
    cant_lecturas[i, "R1"] <- num_readsR1
    cant_lecturas[i, "R2"] <- num_readsR2
    cant_lecturas[i, "Method"] <- de_host_file

    i= i+1
  }

  write.xlsx(cant_lecturas, file = sprintf("~/Daniela/Biota/PipelineBiota/paraPaper/cant_lecturas_%s.xlsx", de_host_file))
  return(cant_lecturas)
}



#problema con 182:
p
deHosting(patient_dir = "/home/daniela/Daniela/Biota/Muestras/73m/182", de_host = "Bowtie")

#################################################################
#% mapeo a humano  ###############################################
##############################################################
patient_dir <- "~/Daniela/Biota/Muestras/73m/140"

calculate_BWA_mapeo <- function(patient_dir) {

  id <- basename(patient_dir)
  patient_dir_trim <- paste0(patient_dir, "/trimmed", sep="")
  file_list_trimmed <- list.files(patient_dir_trim, full.names = TRUE, recursive = FALSE)
  aligned_bam <- path.expand(sprintf("%s/%s_bwa.bam", patient_dir_trim, id))
  #bwa_summary <- path.expand(sprintf("%s/%s_bwa.summary", patient_dir_trim, id))

  if(!file.exists(aligned_bam)) {

    r1_trim <- sprintf("%s/%sT_S04_L001_R1_001.fastq.gz", patient_dir_trim, id)
    r2_trim <- sprintf("%s/%sT_S04_L001_R2_001.fastq.gz", patient_dir_trim, id)
    BWA <- downloadBWA()
    indexBWA <- "/home/daniela/R/x86_64-pc-linux-gnu-library/4.1/PipelineBiota-Softwares/HG38/indexBWA"
    # Guardar el BAM intermedio antes de bedtools
    system2(
      BWA,
      args = c(
        "mem", "-t", "10",
        file.path(indexBWA, "Homo_sapiens.GRCh38.dna_sm.primary_assembly"),
        r1_trim, r2_trim,
        sprintf("| /home/daniela/R/x86_64-pc-linux-gnu-library/4.1/PipelineBiota-Softwares/Samtools/samtools-1.16.1/samtools view -bh -o %s", aligned_bam)
      ))
  }

  # Función para ejecutar comandos y leer salida
  run_cmd <- function(cmd) {
    con <- pipe(cmd)
    output <- readLines(con)
    close(con)
    return(output)
  }
  flagstat_output <- run_cmd(sprintf("/home/daniela/R/x86_64-pc-linux-gnu-library/4.1/PipelineBiota-Softwares/Samtools/samtools-1.16.1/samtools flagstat %s", aligned_bam))
  # Extraer valores clave
  total_fragments <- as.numeric(strsplit(flagstat_output[1], " ")[[1]][1])
  mapped_fragments <- as.numeric(strsplit(flagstat_output[7], " ")[[1]][1])
  porcentaje_mapeo <- round(mapped_fragments/total_fragments*100,2)

  return(porcentaje_mapeo)
}

calculate_mapeo <- function(patients_dir, de_host) {
  mapeo <- data.frame("Sample"= c(), "PorcentajeMapeo" = c(), "Method" = c())

  if(de_host  == "Bowtie") {
    de_host_file <- "DHBo"
  } else if( de_host == "BWA") {
    de_host_file <- "DHbwa"
  } else if(de_host == "RSubread") {
    de_host_file <- "DHRs"
  } else if(de_host == "") {
    de_host_file <- "T"
  } else if(de_host == "DRAGEN") {
    de_host_file <- "DRAGEN"
  }

  i=1
  dir_list <- list.dirs(path = patients_dir, full.names = TRUE, recursive = FALSE)
  cant_patients <- length(dir_list)

  for (p in dir_list) {
    p <- dir_list[[i]]
    print(p)
    id <- basename(p)

    if(de_host == "BWA"){
      porcentaje_mapeo <- calculate_BWA_mapeo(patient_dir = p)

    } else if (de_host == "RSubread") {
      Rs_summary <- as.data.frame(read_delim(sprintf("%s/trimmed/%s_R1R2_Rsubread.bam.summary", p, id)))
      Rs_summary <- rbind(colnames(Rs_summary) , Rs_summary)
      colnames(Rs_summary) <- c("Metrica", "Valor")
      total_fragments <- as.numeric(Rs_summary$Valor[Rs_summary$Metrica == "Total_fragments"])
      mapped_fragments <- as.numeric(Rs_summary$Valor[Rs_summary$Metrica == "Mapped_fragments"])
      porcentaje_mapeo <- round(mapped_fragments/total_fragments*100,2)

    } else if (de_host == "Bowtie") {
      error_file <- sprintf("%s/trimmed/%s_Bo_Error.txt", p, id)
      if(!file.exists(error_file)) {
        dhbo_file_R1 <- sprintf("%s/trimmed/%sDHBo_S04_L001_R1_001.fastq.gz", p, id)
        dhbo_file_R2 <- sprintf("%s/trimmed/%sDHBo_S04_L001_R2_001.fastq.gz", p, id)

        viejo_dhbo_file_R1 <- sprintf("%s/trimmed/viejo_%sDHBo_S04_L001_R1_001.fastq.gz", p, id)
        viejo_dhbo_file_R2 <- sprintf("%s/trimmed/viejo_%sDHBo_S04_L001_R2_001.fastq.gz", p, id)

        file.rename(from = dhbo_file_R1, to = viejo_dhbo_file_R1)
        file.rename(from = dhbo_file_R2, to = viejo_dhbo_file_R2)

        deHosting(patient_dir = p, de_host = "Bowtie")
      }

      Bo_error <- as.data.frame(read_delim(error_file))
      porcentaje_mapeo <- Bo_error[nrow(Bo_error), 1]
      porcentaje_mapeo <- round(as.numeric(gsub("%", "", porcentaje_mapeo)),2)

    } else if (de_host == "DRAGEN") {

      patient_dir <- paste(p, "/trimmed", sep="")
      nombre_p <- sprintf("%sT", id)

      folder_report <- list.dirs( sprintf("%s/DRAGEN_Reports", patient_dir), full.names = TRUE, recursive = FALSE)
      metrics_report <- sprintf("%s/%s.mapping_metrics.csv", folder_report, nombre_p)
      metrics_report <- read_csv(metrics_report)
      metrics_report <- rbind(colnames(metrics_report), metrics_report)
      colnames(metrics_report) <- c("Mapping", "Sample", "Metric", "Value", "Percentage")
      porcentaje_mapeo <- as.numeric(metrics_report$Percentage[metrics_report$Metric == "Mapped reads" & grepl(id, metrics_report$Sample)])
      #num_readsR1 <- round(as.numeric(metrics_report$Value[metrics_report$Metric == "Unmapped reads" & grepl(id, metrics_report$Sample)])/2)
      #num_readsR2 <- num_readsR1
    }

    mapeo[i, "Sample"] <- basename(p)
    mapeo[i, "PorcentajeMapeo"] <- porcentaje_mapeo
    mapeo[i, "Method"] <- de_host_file

    i= i+1
  }

  write.xlsx(mapeo, file = sprintf("~/Daniela/Biota/PipelineBiota/paraPaper/Porcentaje_mapeo_%s.xlsx", de_host_file))
  return(mapeo)
}


patients_dir = "~/Daniela/Biota/Muestras/73m"
mapeo_Rs <- calculate_mapeo(patients_dir = "~/Daniela/Biota/Muestras/73m", de_host = "RSubread")
mapeo_Bo <- calculate_mapeo(patients_dir = "~/Daniela/Biota/Muestras/73m", de_host = "Bowtie")
mapeo_BWA <- calculate_mapeo(patients_dir = "~/Daniela/Biota/Muestras/73m", de_host = "BWA")
mapeo_D <- calculate_mapeo(patients_dir = "~/Daniela/Biota/Muestras/73m", de_host = "DRAGEN")

dif_BWA_Rs <- merge(mapeo_Rs, mapeo_BWA, by = "Sample")
dif_BWA_Rs$Diferencia <- dif_BWA_Rs$PorcentajeMapeo.y - dif_BWA_Rs$PorcentajeMapeo.x
dif_BWA_Rs$DiferenciaPorcentaje <- ((dif_BWA_Rs$PorcentajeMapeo.y - dif_BWA_Rs$PorcentajeMapeo.x) / ((dif_BWA_Rs$PorcentajeMapeo.y + dif_BWA_Rs$PorcentajeMapeo.x) / 2)) * 100

any(is.na(mapeo_Rs))
any(is.na(mapeo_Bo))
any(is.na(mapeo_BWA))
any(is.na(mapeo_D))

#problema con 130:
p
deHosting(patient_dir = "/home/daniela/Daniela/Biota/Muestras/73m/185", de_host = "Bowtie")

mapeo_totales <- rbind( mapeo_Bo, mapeo_BWA, mapeo_Rs, mapeo_D)
mapeo_totales$Method[mapeo_totales$Method == "DHBo"] <- "Bowtie"
mapeo_totales$Method[mapeo_totales$Method == "DHRs"] <- "RSubread"
mapeo_totales$Method[mapeo_totales$Method == "DHbwa"] <- "BWA"

colnames(mapeo_totales)[2] <- "MappedReads"

write.xlsx(mapeo_totales, file = "~/Daniela/Biota/PipelineBiota/paraPaper/Scripts Reproducir Paper/Fig2a-MappedReads.xlsx")

MappedReads_AllMethods <- mapeo_totales

#mas estetico:
unique(MappedReads_AllMethods$Method )
MappedReads_AllMethods$Method <- factor(MappedReads_AllMethods$Method, levels = c("Bowtie","BWA", "RSubread", "DRAGEN"))

ggplot(MappedReads_AllMethods, aes(x = Method, y = MappedReads, fill = Method)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = NA) +  # Violin plot con transparencia
  geom_jitter(aes(color = Method), width = 0.2, size = 1.2, alpha = 0.7) +  # Jitter para dispersión
  geom_boxplot(width = 0.6, outlier.shape = NA, color = "black", alpha = 0.8) +  # Boxplot sin outliers visibles
  geom_line(aes(group=Sample), colour="black", linetype="11", alpha = 0.3) +
  theme_minimal() +
  labs(x = "", y = "% Mapped Reads", title = "") +
  scale_y_continuous(breaks = seq(0, 100, by =20))+
  scale_fill_manual(values = c(
    "Bowtie" = "#E69F00",  # Naranja suave
    "BWA" = "#56B4E9",  # Azul suave
    "RSubread" = "#009E73",  # Verde suave
    "DRAGEN" = "#9467BD"   # Violeta suave
  )) +
  scale_color_manual(values = c(
    "Bowtie" = "#E69F00",
    "BWA" = "#56B4E9",
    "RSubread" = "#009E73",
    "DRAGEN" = "#9467BD"
  )) +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 12),
    legend.position = "none"
  ) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))  # Mantiene alineación del eje X

#Figure 2a - osea la nterior pero con àtterns:
# Instalar y cargar ggpattern si no lo tienes
install.packages("ggpattern")
library(ggpattern)

# Instalar y cargar ggpattern si no lo tienes
install.packages("ggpattern")
library(ggpattern)

# Crear el gráfico con patrones solo en el boxplot
mapeo_totales$Method
ggplot(MappedReads_AllMethods, aes(x = Method, y = MappedReads, pattern = Method)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = "grey15") +  # Gráfico de violín sin patrones
  #geom_jitter(color = "grey20", width = 0.2, size = 1.2, alpha = 0.7) +  # Jitter para dispersión
  geom_boxplot_pattern(
    width = 0.6, outlier.shape = NA,
    color = "black", alpha = 0.8,
    pattern_density = 0.1
  ) +  # Boxplot con patrones
  geom_line(aes(group = Sample), colour = "black", linetype = "11", alpha = 0.3) +
  theme_minimal() +
  labs(x = "", y = "% Mapped Reads", title = "") +
  scale_y_continuous(breaks = seq(0, 100, by = 20)) +
  scale_pattern_manual(values = c(
    "Bowtie" = "stripe",  # Patrones de rayas
    "BWA" = "crosshatch",  # Patrones de círculos
    "RSubread" = "circle",  # Patrones de líneas cruzadas
    "DRAGEN" = "none"  # Patrones de cuadrados
  )) +

  theme(
    text = element_text(size = 14),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 1, face = "bold"),
    axis.title.y = element_text(size = 12),
    legend.position = "none"
  ) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))  # Mantiene alineación del eje X

#P-values:

library(dplyr)
library(tidyr)
library(purrr)

df_wide <- MappedReads_AllMethods %>%
  pivot_wider(names_from = Method, values_from = MappedReads)

df_wide$Dif_BWA_Rs <- ((df_wide$BWA - df_wide$RSubread ) / df_wide$BWA)*100

methods <- setdiff(colnames(df_wide), "Sample")
method_pairs <- combn(methods, 2, simplify = FALSE)

wilcox_results <- map(method_pairs, function(pair) {
  test <- wilcox.test(df_wide[[pair[1]]], df_wide[[pair[2]]], paired = FALSE)
  tibble(
    Method1 = pair[1],
    Method2 = pair[2],
    p_value = test$p.value
  )
}) %>% bind_rows()

wilcox_results <- wilcox_results %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))  # o "bonferroni"


###########################################################################
#Abundancia relativa que se clasifica como eukaryota despues del dehost
###############################################################################

patients_dir = "~/Daniela/Biota/Muestras/73m"

homo_KBo <- extraerAR_HomoSapiens(patients_dir = patients_dir, source = "KRAKEN", de_host = "Bowtie")
homo_KBWA <-extraerAR_HomoSapiens(patients_dir = patients_dir, source = "KRAKEN", de_host = "BWA")
homo_KRs <-extraerAR_HomoSapiens(patients_dir = patients_dir, source = "KRAKEN", de_host = "RSubread")
homo_Ksin <-extraerAR_HomoSapiens(patients_dir = patients_dir, source = "KRAKEN", de_host = "")

homo_Dsin <- extraerAR_HomoSapiens(patients_dir = patients_dir, source = "DRAGEN", de_host = "")
homo_DsinDH_PD <-extraerAR_HomoSapiens(patients_dir = patients_dir, source = "DRAGEN", de_host = "sinDH_PD")
homo_DBo <-extraerAR_HomoSapiens(patients_dir = patients_dir, source = "DRAGEN", de_host = "Bowtie")
homo_DBWA <-extraerAR_HomoSapiens(patients_dir = patients_dir, source = "DRAGEN", de_host = "BWA")
homo_DRs <-extraerAR_HomoSapiens(patients_dir = patients_dir, source = "DRAGEN", de_host = "RSubread")


extraerAR_HomoSapiens <- function(patients_dir, source, de_host) {
  tabla_otus <- generateOTUsTableGrupal(patients_dir = patients_dir, source = source, conEukaryota = TRUE, de_host = de_host)
  out <- group_TaxonomicLevels(patients_dir = patients_dir, tabla_otus = tabla_otus, source = source, de_host = de_host, conE = TRUE)
  list_AR <- out[[1]]
  #AR_filos <- list_AR[[3]]
  #AR_generos <- list_AR[[7]]
  AR_species <- list_AR[[8]]
  colnames(AR_species)[-1] <- gsub(sprintf("_%s", source), "", colnames(AR_species)[-1])

  AR_homosapiens <- AR_species[AR_species$Species == "Homo sapiens",]
  s <- ifelse(source == "DRAGEN", "D", "K")
  AR_homosapiens$Method <- paste0(s, de_host, sep ="")

  write.xlsx(AR_homosapiens, file = sprintf("~/Daniela/Biota/PipelineBiota/paraPaper/AR_homosapiens_%s_%s.xlsx", s, de_host))
  return(AR_homosapiens)
}

todos_homosapiens <- rbind(homo_KBo, homo_KBWA, homo_KRs, homo_Ksin,
                          homo_DBo, homo_DBWA, homo_DRs, homo_Dsin, homo_DsinDH_PD)
df_long <- todos_homosapiens %>%
  pivot_longer(cols = -c(Method, Species), names_to = "Sample", values_to = "Abundance")

unique(df_long$Method)

df_long$Method[df_long$Method == "D"] <- "DdhD"
df_long$Method[df_long$Method == "DsinDH_PD"] <- "-D"
df_long$Method[df_long$Method == "KBowtie"] <- "BoK"
df_long$Method[df_long$Method == "KBWA"] <- "bwaK"
df_long$Method[df_long$Method == "KRSubread"] <- "RsK"
df_long$Method[df_long$Method == "K"] <- "-K"
df_long$Method[df_long$Method == "DBowtie"] <- "BoD"
df_long$Method[df_long$Method == "DBWA"] <- "bwaD"
df_long$Method[df_long$Method == "DRSubread"] <- "RsD"

df_long$Method <- factor(df_long$Method, levels = c("BoK", "BoD",
                                                    "bwaK", "bwaD", "RsK", "RsD",
                                                    "DdhD", "-K","-D" ))

df_long$Method <- factor(df_long$Method, levels = c("-K","-D", "DdhD",  "BoK", "BoD",
                                                    "bwaK", "bwaD", "RsK", "RsD"))


# mas estetico:

write.xlsx(df_long, file = "~/Daniela/Biota/PipelineBiota/paraPaper/Scripts Reproducir Paper/Fig2c-HomoSapiens.xlsx")
library(ggplot2)

HomoSapiens_AllMethods <- df_long

ggplot(HomoSapiens_AllMethods, aes(x = Method, y = Abundance, fill = Method)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = NA) +  # Violin con transparencia
  geom_jitter(aes(color = Method), width = 0.2, size = 1.2, alpha = 0.4) +  # Jitter sutil y transparente
  geom_line(aes(group=Sample), colour="grey50", linetype="11", alpha = 0.2) +
  geom_boxplot(width = 0.6, outlier.shape = NA, color = "black", alpha = 0.8) +  # Boxplot clásico
  theme_minimal() +
  labs(x = "", y = "Relative Abundance (log10)", title = "") +
  scale_y_log10() +
  scale_fill_manual(values = c(
    "BoK" = "#E69F00",  # Naranja suave
    "BoD" = "#F0C987",

    "bwaK" = "#56B4E9", # Azul suave
    "bwaD" = "#A0D8F0",

    "RsK" = "#009E73",  # Verde suave
    "RsD" = "#8FCB88",

    "-K" = "#666666",  # Gris oscuro
    "DdhD" = "#9467BD", # Violeta suave
    "-D" = "#999999"
  )) +
  scale_color_manual(values = c(
    "BoK" = "#E69F00", "BoD" = "#F0C987",
    "bwaK" = "#56B4E9", "bwaD" = "#A0D8F0",
    "RsK" = "#009E73", "RsD" = "#8FCB88",
    "-K" = "#666666", "DdhD" = "#9467BD", "-D" = "#999999"
  )) +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 12),
    legend.position = "none"
  ) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))  # Rotar etiquetas suavemente

#Con patrones:
df_long$DeHost <- ifelse(grepl("Bo", df_long$Method), "Bowtie",
                         ifelse(grepl("bwa", df_long$Method), "BWA",
                                ifelse(grepl("-", df_long$Method), "None",
                                       ifelse(grepl("Rs", df_long$Method), "RSubread", "DRAGEN"))))

df_long$DeHost <- factor(df_long$DeHost, levels = c("Bowtie", "BWA", "RSubread", "DRAGEN", "None"))
df_long$Tax <- ifelse(grepl("K", df_long$Method), "Kraken", "DRAGEN")
df_long$Tax  <- as.factor(df_long$Tax )

ggplot(df_long, aes(x = Method, y = Abundance, pattern = Tax)) +
  geom_violin(trim = FALSE, alpha = 0.8, color = "grey15") +  # Violin con transparencia
  geom_jitter(color = "black", width = 0.1, size = 1, alpha = 0.2) +  # Jitter sutil y transparente
  geom_boxplot_pattern(
    aes(pattern = Tax),
    pattern_fill = "black",
    pattern_density = 0.4,
    pattern_spacing = 0.05,
    pattern_angle = 45,
    color = "black",
    alpha = 0.7,
    position = position_dodge(width = 0.8)
  ) +
  theme_minimal() +
  labs(x = "", y = "Relative Abundance (log10)", title = "") +
  scale_y_log10() +
  scale_pattern_manual(
    name = "Tax Classifier",
    values = c( "Kraken" = "none", "DRAGEN" = "stripe"  )
    )+
  theme(
    text = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 1, face = "bold"),
    axis.title.y = element_text(size = 11),
    legend.position = "right"
  )
#P-valores:

library(dplyr)
library(tidyr)
library(purrr)

# Todos los pares de métodos
method_pairs <- combn(unique(df_long$Method), 2, simplify = FALSE)

# Función para comparar un par de métodos
dar_pvalor <- function(df, method1, method2) {
  df_pair <- df %>%
    filter(Method %in% c(method1, method2)) %>%
    pivot_wider(names_from = Method, values_from = Abundance)

  test <- wilcox.test(df_pair[[method1]], df_pair[[method2]], paired = TRUE)

  diff_pct <- mean(((df_pair[[method1]] - df_pair[[method2]]) / df_pair[[method1]]) * 100)

  tibble(
    Method1 = method1,
    Method2 = method2,
    p_value = test$p.value,
    diff_pct = diff_pct
  )
}

# Aplicar a todos los pares de métodos
wilcoxon_results <- map_dfr(method_pairs, ~ dar_pvalor(df_long, .x[1], .x[2])) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))


