.onLoad <- function(libname, pkgname) {
  
  message("ESTOY EN EL ONLOAD")
  
  libPath <- dirname(system.file(package = "PipelineBiota"))
  print(libPath)
  
  #libPath <- Sys.getenv('R_LIBS_USER')
  
  #Folder where the softwares will be saved
  if(!(file.exists(sprintf("%s/PipelineBiota-Softwares", libPath)))) {
    dir.create(sprintf("%s/PipelineBiota-Softwares", libPath))
  }
  
  if (!(file.exists(sprintf("%s/PipelineBiota-Softwares/path_to_soft.txt", libPath)))) { #Solo en instalacion
    write("",file = sprintf("%s/PipelineBiota-Softwares/path_to_soft.txt", libPath))
  }
  
  omicsdo_sof <<- sprintf("%s/PipelineBiota-Softwares", dirname(system.file(package = "PipelineBiota")))
  
  #check_packages()
  
  #downloadFastQC()
  #downloadTrimGalore()
  #downloadSamtools()
  #downloadBWA()
  #downloadHG38()
  
}

check_packages <- function() {
  
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    remotes::install_cran("openxlsx")
  }
  library(openxlsx)
  
  if (!requireNamespace("readr", quietly = TRUE)) {
    remotes::install_cran("readr")
  }
  library(readr)
  
  if (!requireNamespace("stringr", quietly = TRUE)) {
    remotes::install_cran("stringr")
  }
  library(stringr)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    remotes::install_cran("ggplot2")
  }
  library(ggplot2)
  
  if (!requireNamespace("tibble", quietly = TRUE)) {
    remotes::install_cran("tibble")
  }
  library(tibble)
  
  if (!requireNamespace("tidyverse", quietly = TRUE)) {
    remotes::install_cran("tidyverse")
  }
  library(tidyverse)
  
  if (!requireNamespace("showtext", quietly = TRUE)) {
    remotes::install_cran("showtext")
  }
  
  library(showtext)
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    remotes::install_cran("dplyr")
  }
  library(dplyr)
  
  if (!requireNamespace("ggtext", quietly = TRUE)) {
    remotes::install_cran("ggtext")
  }
  library(ggtext)
  
  if (!requireNamespace("readxl", quietly = TRUE)) {
    remotes::install_cran("readxl")
  }
  library(readxl)
  
  if (!requireNamespace("magrittr", quietly = TRUE)) {
    remotes::install_cran("magrittr")
  }
  library(magrittr)
  
  if (!requireNamespace("httr", quietly = TRUE)) {
    remotes::install_cran("httr")
  }
  library(httr)
  
  if (!requireNamespace("Rsamtools", quietly = TRUE)) {
    BiocManager::install("Rsamtools")
  }
  library(Rsamtools)
  
  if (!requireNamespace("viridis", quietly = TRUE)) {
    install.packages("viridis")
  }
  library(viridis)
  
  if (!requireNamespace("reshape2", quietly = TRUE)) {
    install.packages("reshape2")
  }
  library(reshape2)
  
  #Para usar MIXTURE
  if (!requireNamespace("nnls", quietly = TRUE)) {
    install.packages("nnls")
  }
  library(nnls)
  
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    BiocManager::install("ComplexHeatmap", force = TRUE)
  }
  library(ComplexHeatmap)
  
  #if (!requireNamespace("MIXTURE", quietly = TRUE)) {
  #  install_github("elmerfer/MIXTURE")
  #}
  #library(MIXTURE)
  
  if (!requireNamespace("Rsubread", quietly = TRUE)) {
    BiocManager::install("Rsubread")
  }
  library(Rsubread)
  
  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    BiocManager::install("GEOquery")
  }
  library(GEOquery)
  
  message("The R packages required have been successfully installed")
}
