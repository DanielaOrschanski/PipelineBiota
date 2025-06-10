# Definir rutas y parámetros
library(PipelineBiota)

#instalar bracken:
setwd("/home/daniela/R/x86_64-pc-linux-gnu-library/4.1/PipelineBiota-Softwares")
system("git clone https://github.com/jenniferlu717/Bracken")
setwd("/home/daniela/R/x86_64-pc-linux-gnu-library/4.1/PipelineBiota-Softwares/Bracken")
system("bash install_bracken.sh")
chmod +x ~/R/x86_64-pc-linux-gnu-library/4.1/PipelineBiota-Softwares/Bracken/bracken
chmod +x ~/R/x86_64-pc-linux-gnu-library/4.1/PipelineBiota-Softwares/Bracken/bracken-build
export PATH=$PATH:~/R/x86_64-pc-linux-gnu-library/4.1/PipelineBiota-Softwares/Bracken
export PATH=$PATH:~/R/x86_64-pc-linux-gnu-library/4.1/PipelineBiota-Softwares/Bracken/bracken-build
export PATH=$PATH:~/R/x86_64-pc-linux-gnu-library/4.1/PipelineBiota-Softwares/KRAKEN2/kraken2/kraken2

source ~/.bashrc


#-----------------------

system("cp -r ~/R/x86_64-pc-linux-gnu-library/4.1/PipelineBiota-Softwares/KRAKEN2/minikraken2_v2_8GB_201904_UPDATE ~/Daniela/Biota/PipelineBiota-Softwares/KRAKEN2/backup_kraken2_db")
sudo apt-get remove kraken2

#-------------------------------

path_kraken2 <- "/usr/bin/kraken2"
system("/usr/bin/kraken2 --version")

ruta_kraken2_db <- "~/R/x86_64-pc-linux-gnu-library/4.1/PipelineBiota-Softwares/KRAKEN2/minikraken2_v2_8GB_201904_UPDATE"

ruta_bracken <- "/home/daniela/R/x86_64-pc-linux-gnu-library/4.1/PipelineBiota-Softwares/Bracken/bracken"
ruta_bracken_build <- "/home/daniela/R/x86_64-pc-linux-gnu-library/4.1/PipelineBiota-Softwares/Bracken/bracken-build"
reporte_kraken2 <- " ~/Daniela/Biota/Muestras/73m/1/trimmed/Resultados_KRAKEN/report_Bo.sequences"
reporte_bracken <- " ~/Daniela/Biota/Muestras/73m/1/trimmed/Resultados_KRAKEN/report_Bo.bracken"
nivel_taxonomico <- "S"  # Ajusta según el nivel deseado (S para especie, G para género, etc.)

# Crear índice de Bracken
bracken_build_command <- paste(
  #ruta_bracken,
  #"bracken-build",
  ruta_bracken_build,
  "-x", path_kraken2,
  "-d", ruta_kraken2_db,
  "-t", 4,
  "-k", 35,
  "-l", 100
)

bracken_build_command <- paste(
  ruta_bracken_build,
  "-x", "/home/daniela/Daniela/Biota/Kraken/kraken2/kraken2",
  "-d", ruta_kraken2_db,
  "-t", 4,
  "-k", 35,
  "-l", 100
)


system(bracken_build_command)


# Ejecutar Bracken en el reporte de Kraken2
runBracken <- function() {
  ruta_kraken2_db <- "~/R/x86_64-pc-linux-gnu-library/4.1/PipelineBiota-Softwares/KRAKEN2/minikraken2_v2_8GB_201904_UPDATE"

  ruta_bracken <- "/home/daniela/R/x86_64-pc-linux-gnu-library/4.1/PipelineBiota-Softwares/Bracken/bracken"
  ruta_bracken_build <- "/home/daniela/R/x86_64-pc-linux-gnu-library/4.1/PipelineBiota-Softwares/Bracken/bracken-build"
  reporte_kraken2 <- " ~/Daniela/Biota/Muestras/73m/1/trimmed/Resultados_KRAKEN/report_Bo.sequences"
  reporte_bracken <- " ~/Daniela/Biota/Muestras/73m/1/trimmed/Resultados_KRAKEN/report_Bo.bracken"
  nivel_taxonomico <- "S"  # Ajusta según el nivel deseado (S para especie, G para género, etc.)

  system(paste(ruta_bracken, "-d", ruta_kraken2_db, "-i", reporte_kraken2, "-o", reporte_bracken,  "-l", nivel_taxonomico))

}

setwd("~/Daniela/Biota/Muestras/73m")

# Comparación Kraken vs Bracken: ver el impacto del ajuste

patient_dir <- "~/Daniela/Biota/Muestras/73m/1/trimmed"
de_host_file <- "Bo"
library(data.table)
Species_Bracken <- fread(paste(patient_dir, sprintf("/Resultados_KRAKEN/report_%s.bracken", de_host_file), sep=""),  header = FALSE,sep = "\t")
Species_Bracken <- as.data.frame(Species_Bracken)
colnames(Species_Bracken) <- Species_Bracken[1,]
Species_Bracken <- Species_Bracken[-1,]
Species_Bracken$AR_Bracken <- as.numeric(Species_Bracken$fraction_total_reads)*100
sum(Species_Bracken$AR_Bracken)
Species_Bracken$AR_Kraken <- as.numeric(Species_Bracken$kraken_assigned_reads)/sum(as.numeric(Species_Bracken$kraken_assigned_reads))*100

#BLAND ALTMAN  -----------------------
library(tidyr)
library(ggplot2)
library(dplyr)


library(ggplot2)

Species_Bracken$kraken_assigned_reads <- as.numeric(Species_Bracken$kraken_assigned_reads )
Species_Bracken$new_est_reads <- as.numeric(Species_Bracken$new_est_reads)

# Calcular medias y diferencias
Species_Bracken$mean <- (Species_Bracken$AR_Bracken + Species_Bracken$AR_Kraken) / 2
Species_Bracken$diff <- Species_Bracken$AR_Bracken - Species_Bracken$AR_Kraken

Species_Bracken$mean <- (Species_Bracken$new_est_reads + Species_Bracken$kraken_assigned_reads) / 2
Species_Bracken$diff <- Species_Bracken$new_est_reads - Species_Bracken$kraken_assigned_reads

# Calcular media y límites de acuerdo
mean_diff <- mean(Species_Bracken$diff, na.rm = TRUE)
sd_diff <- sd(Species_Bracken$diff, na.rm = TRUE)
upper_limit <- mean_diff + 1.96 * sd_diff
lower_limit <- mean_diff - 1.96 * sd_diff

# Graficar Bland-Altman
ggplot(Species_Bracken, aes(x = mean, y = diff)) +
  geom_point(alpha = 0.5, color = "blue") +  # Puntos de datos
  geom_hline(yintercept = mean_diff, linetype = "solid", color = "red") +  # Línea de media
  geom_hline(yintercept = upper_limit, linetype = "dashed", color = "black") +  # Límite superior
  geom_hline(yintercept = lower_limit, linetype = "dashed", color = "black") +  # Límite inferior
  labs(title = "Bland-Altman Plot: CONTEOS Bracken vs Kraken - Muestra 1 - 150 especies",
       x = "Media de Bracken y Kraken",
       y = "Bracken - Kraken") +
  theme_minimal()



as.numeric(Species_Bracken$new_est_reads[1])/sum(as.numeric(Species_Bracken$new_est_reads))*100
Species_Bracken$fraction_total_reads[1]



