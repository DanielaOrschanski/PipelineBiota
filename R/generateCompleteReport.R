#' @title Generate Complete Results' Report
#' @description It generates the report in pdf format that register all the predictions and the analizes of the microbiome composition of a determined patient.
#' @param patients_dir indicated the patient directory that will be analyzed.
#' @export
#' @import rmarkdown
#' @import ggplot2
#' @import kableExtra
#' @import gridExtra
#' @import readxl
#' @import data.table
#' @import openxlsx

generateCompleteReport <- function(patient_dir) {

  library(qpdf)
  id <- basename(patient_dir)
  id_mr <- str_pad(id, width = 4, side = "left", pad = "0")

  caratulaComplete(patient_dir = patient_dir)
  indice(patient_dir = patient_dir, path_metadata = "/media/4tb2/Daniela/Biota/PipelineBiota-master/data/Metadata-Completa-SinLimpiar.xlsx")
  contenidoComplete(patient_dir = patient_dir)

  contenido_file <- sprintf("%s/CompleteContenido_ID%s.pdf", patient_dir, id_mr)
  caratula_file <- sprintf("%s/CompleteCaratula_ID%s.pdf", patient_dir, id_mr)
  indice_file <- sprintf("%s/indice_MR_%s.pdf", patient_dir, id_mr)

  output_file <- sprintf("%s/Doctor_MR_%s.pdf", patient_dir, id_mr)

  pdf_combine(c(caratula_file, indice_file, contenido_file), output_file)


  #Elimino todo lo que se genera y no me sirve:
  file.remove(sprintf("%s/cantidades_barra.png", patient_dir))
  file.remove(sprintf("%s/score_barra.png", patient_dir))
  file.remove(sprintf("%s/rango_etario_barra.png", patient_dir))
  file.remove(sprintf("%s/fagos_barra.png", patient_dir))
  file.remove(sprintf("%s/virus_barra.png", patient_dir))
  file.remove(sprintf("%s/tipodepiel_barra.png", patient_dir))

  file.remove(sprintf("%s/Tabla_Genus_masAbundantes.xlsx", patient_dir))
  file.remove(sprintf("%s/Tabla_Phylum_masAbundantes.xlsx", patient_dir))
  file.remove(sprintf("%s/Tabla_Species_masAbundantes.xlsx", patient_dir))
  file.remove(sprintf("%s/Tabla_SubSpecies_masAbundantes.xlsx", patient_dir))
  file.remove(sprintf("%s/Tabla_SubEspecies_RangoEtario.xlsx", patient_dir))
  file.remove(sprintf("%s/Tabla_CantidadesTax.xlsx", patient_dir))

  file.remove(sprintf("%s/CompleteContenido_ID%s.Rmd", patient_dir, id_mr))
  file.remove(sprintf("%s/CompleteCaratula_ID%s.Rmd", patient_dir, id_mr))
  file.remove(sprintf("%s/CompleteCaratula_ID%s.log", patient_dir, id_mr))
  file.remove(sprintf("%s/CompleteContenido_ID%s.pdf", patient_dir, id_mr))
  file.remove(sprintf("%s/CompleteCaratula_ID%s.pdf", patient_dir, id_mr))
  file.remove(sprintf("%s/indice_MR_%s.log", patient_dir, id_mr))
  file.remove(sprintf("%s/indice_MR_%s.pdf", patient_dir, id_mr))

  #Esto despues se va:

  file.remove(sprintf("%s/habitos_barra.png", patient_dir))

  #-----------------------

  file.remove(sprintf("%s/Doctor_MR_%s.log", patient_dir, id_mr))


}

generateCompleteReportV <- function(patient_dir) {

  id <- basename(patient_dir)
  patients_dir <- dirname(patient_dir)
  MetadataB <- as.data.frame(read_excel("~/Daniela/Biota/PipelineBiota/data/Metadata-soloColumnasUsables.xlsx"))
  BIOTALIFE_SKIN_Respuestas_ <- read_excel("~/Daniela/Biota/BIOTALIFE SKIN  (Respuestas).xlsx", sheet = 4)
  #BIOTALIFE_SKIN_Respuestas_ <- read_excel("~/Daniela/Biota/BIOTALIFE SKIN- ActualizadoOct.xlsx")

  #Obtener info gral
  BIOTALIFE_SKIN_Respuestas_$ID  <- gsub("\\.0$", "", as.character(BIOTALIFE_SKIN_Respuestas_$ID))
  BIOTALIFE_SKIN_Respuestas_$ID[which(BIOTALIFE_SKIN_Respuestas_$ID == "118")] <- "118-1"
  BIOTALIFE_SKIN_Respuestas_$ID[which(BIOTALIFE_SKIN_Respuestas_$ID == "184")] <- "184-1"

  info <- BIOTALIFE_SKIN_Respuestas_[which(BIOTALIFE_SKIN_Respuestas_$ID == id),]
  colnames(info)
  info <- info[, which(colnames(info) %in% c("ID", "1. Nombre y Apellido","3. Sexo", "2. Fecha de Nacimiento", "Edad", "Fecha Cita" ))]
  info <- as.data.frame(t(info))
  info <- data.frame("Característica" = rownames(info), "Información" = info$V1)
  info$Característica <- c("ID", "Nombre y Apellido", "Fecha de nacimiento", "Sexo", "Fecha Cita", "Edad")

  nombre <- info$Información[which(info$Característica == "Nombre y Apellido")]
  capitalizar_apellido <- function(nombre) {
    partes <- strsplit(nombre, " ")[[1]]  # Divide el nombre en partes
    partes[2] <- paste0(toupper(substring(partes[2], 1, 1)), substring(partes[2], 2))  # Capitaliza la primera letra del apellido
    return(paste(partes, collapse = " "))  # Vuelve a unir las partes
  }
  nombre_modificado <- capitalizar_apellido(nombre)

  info$Información[which(info$Característica == "Nombre y Apellido")] <- nombre_modificado

  #Agrego info del paciente si es que NO esta en MetadabaB pero SI en BiotaLife respuestas: --------------------------
  MetadataB <- MetadataB[, -which(colnames(MetadataB) == "TipoPielViejo")]
  MetadataB <- MetadataB[, -which(colnames(MetadataB) == "Secuenciado")]
  colnames(BIOTALIFE_SKIN_Respuestas_)[c(2,5,6,7,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,25,26,27,28,29,30,31,33,34,57, 60, 61, 32,35)] <- c("ID", "FechadeNacimiento", "Sexo", "Email",
                                                                                                                                                 "ColorCabello", "ColorOjos", "FacilidadBroncearse", "Peso",
                                                                                                                                                 "Altura", "AntecedenteEnfermedad", "CualEnfermedad", "AfeccionesPiel",
                                                                                                                                                 "OtraAfeccionPiel", "TratamientoEstético3meses", "MétodoAnticonceptivo", "EmbarazadaoAmantando",
                                                                                                                                                 "TratamientoMédico", "CualTratamientoMedico", "FrecuenciaTratamientoMedico", "Tabaco",
                                                                                                                                                 "CigarrillosporSemana", "Alcohol", "AlcoholporSemana", "ActividadFísica",
                                                                                                                                                 "CuándoLimpiaCara", "ConquéLimpiaCara", "AplicacionProtectorSolar", "CuándoMaquillaje",
                                                                                                                                                 "Fecha Cita", "Edad", "Rango etario",
                                                                                                                                                 "Tipodepiel", "Maquillaje_Base")

  if(!(id %in% MetadataB$ID) & id %in% BIOTALIFE_SKIN_Respuestas_$ID) {
    MetadataB <- rbind(MetadataB, BIOTALIFE_SKIN_Respuestas_[which(BIOTALIFE_SKIN_Respuestas_$ID == id), which(colnames(BIOTALIFE_SKIN_Respuestas_) %in% colnames(MetadataB))])
    #modifico lo de maquillaje base y lo de facilidad de bronceado
    MetadataB$FacilidadBroncearse <- sub("^(\\d+).*", "\\1", MetadataB$FacilidadBroncearse)
    MetadataB$FacilidadBroncearse[which(MetadataB$FacilidadBroncearse == "1")] <- "2"
    MetadataB$CuándoMaquillaje[which(MetadataB$CuándoMaquillaje == "1 o 2 veces por semana")] <- "Diariamente"
    MetadataB$AplicacionProtectorSolar[which(MetadataB$AplicacionProtectorSolar == "1 o 2 veces por semana")] <- "Solo ante exposición en verano"

    if(any(MetadataB$Maquillaje_Base != "0" | MetadataB$Maquillaje_Base != "1")) {
      maq_base <-  MetadataB$Maquillaje_Base[-which(MetadataB$Maquillaje_Base == "0" | MetadataB$Maquillaje_Base == "1")]
      MetadataB$Maquillaje_Base[which(MetadataB$Maquillaje_Base == maq_base)] <- ifelse(grepl("base", maq_base) | grepl("corrector", maq_base) | grepl("polvo", maq_base), "1", "0")
    }
    MetadataB$Maquillaje_Base[which(is.na(MetadataB$Maquillaje_Base))] <- "0"
  }
  #----------------------------------------------------------------------------------------------------------


  # Obtener datos de las tablas
  counts_tax <- counts_Tax(patients_dir = patient_dir, source = "KRAKEN", de_host = "Biota", conEukaryota = FALSE)

  reportCountsTax(id=id, MetadataB = MetadataB, patients_dir = patients_dir, de_host = "Biota", conEukaryota = FALSE)
  reportMasAbundantes(id = id, MetadataB = MetadataB)
  indicadorRangoEtario(id = id, MetadataB = MetadataB)
  titulo_titulos <- "Qué podés encontrar en este informe:"
  titulo1 <- "1. ¿Qué es el microbioma y por qué es importante para la salud de tu piel?"
  texto1 <- "El microbioma de la piel es una comunidad compleja de microorganismos, como bacterias, hongos y virus, que viven en su superficie. Estos microorganismos cumplen funciones importantes como: proteger la piel de agentes dañinos, mantener la hidratación y el equilibrio del pH, y fortalecer el sistema inmunológico. Factores como la edad, el sexo, la dieta y el estilo de vida pueden influir en la composición del microbioma. Mantenerlo equilibrado es esencial para una piel sana, ya que su desequilibrio puede contribuir a problemas como el acné, la dermatitis atópica, la rosácea y el envejecimiento prematuro."

  titulo2 <- " 2. ¿Cómo interpretar tus resultados?"
  texto2 <- "En este informe, se analiza la cantidad de microorganismos que componen el microbioma de tu piel. Se evalúan en diferentes niveles taxonómicos (filo, clase, orden, familia, género y especie) y se clasifica según su presencia en comparación con una base de datos de referencia de piel saludable de Biotalife Skin. Es importante destacar que un valor fuera del rango no significa necesariamente un problema de salud, ya que el microbioma puede variar por factores temporales como cambios en el entorno, uso de productos cosméticos o consumo de medicamentos."

  titulo3 <- "3. ¿En qué consiste nuestro estudio?"
  texto3 <- "Este estudio utiliza tecnología de secuenciación de nueva generación (NGS) de Illumina para analizar el ADN de los microorganismos presentes en una muestra de la piel obtenida mediante un hisopado estéril. La secuenciación permite identificar con precisión los diferentes microorganismos, brindando una visión detallada de la composición del microbioma de la piel."

  titulo4 <- "4. Descargo de responsabilidad"
  texto4 <- "Cada parámetro evaluado se clasifica como ‘Bajo’, ‘Equilibrado’ o ‘Alto’ utilizando un algoritmo desarrollado por nuestro equipo. Además, identificamos ciertas especies que pueden actuar como indicadores de la edad de la piel, lo que nos permite predecir si el microbioma corresponde a tu edad biológica. Es importante tener en cuenta que estos resultados son informativos, no es un diagnóstico médico ni reemplaza a una consulta médica. "

  titulo5 <- "5. Descripción de la composición de tu microbioma"
  texto5 <- "Tu microbioma es único, como tu huella digital. En esta sección, detallaremos la variabilidad y cantidad de microorganismos, destacando los grupos más abundantes.  Se mostrarán índices que van del 0 al 1, presentando resultados más altos en microbiomas más equilibrados, lo que contribuye a mantener una piel saludable."
  titulo_score <- "Índice de microorganismos más abundantes"
  texto_final_score <- "A continuación, se detallan los microorganismos más abundantes a nivel de Filo, Género y Especie (ver gráfico):"

  subtitulo_filos <- "Filos"
  texto_filos <- "Los FILOS bacterianos más comunes en la piel son Actinobacteria,  Firmicutes, Proteobacteria y Bacteroidetes."
  Filos_masAbundantes <- read_excel(sprintf("%s/Tabla_Phylum_masAbundantes.xlsx", patient_dir))

  subtitulo_generos <- "Géneros"
  Genus_masAbundantes <- read_excel(sprintf("%s/Tabla_Genus_masAbundantes.xlsx", patient_dir))
  generos_top3 <- paste0(Genus_masAbundantes$Géneros[1:3], collapse = ", ")
  #if(any(grepl("virus", Genus_masAbundantes$Géneros))) {
  #  warning_virus_generos <- "La presencia de virus dentro de los géneros más abundantes de la piel es un factor de riesgo para la salud cutánea."
  #} else {
  #  warning_virus_generos <- ""
  #}
  texto_generos <- sprintf("Los GÉNEROS más abundantes de tu piel son %s.", generos_top3)


  subtitulo_especies <- "Especies"
  esp_masAbundantes <- read_excel(sprintf("%s/Tabla_Species_masAbundantes.xlsx", patient_dir))
  especies_top3 <- paste0(esp_masAbundantes$Especies[1:3], collapse = ", ")
  #if(any(grepl("virus", esp_masAbundantes$Especies))) {
  #  warning_virus_especies <- "La presencia de virus dentro de las especies más abundantes de la piel es un factor de riesgo para la salud cutánea."
  #} else {
  #  warning_virus_especies <- ""
  #}
  texto_especies <- sprintf("Las ESPECIES más abundantes de tu piel son %s.", especies_top3)

  #subtitulo_subespecies <- "5.4. Subespecies más abundantes"
  #subesp_masAbundantes <- read_excel(sprintf("%s/Tabla_SubSpecies_masAbundantes.xlsx", patient_dir))
  #if(any(grepl("virus", subesp_masAbundantes$SubEspecies))) {
  #  warning_virus <- "La presencia de virus dentro de las subespecies más abundantes de la piel es un factor de riesgo para la salud cutánea. Se recomienda solicitar una consulta con un dermatólogo."
  #} else {
  #  warning_virus <- ""
  #}
  #texto_subespecies <- sprintf("En algunos casos, se pueden identificar subespecies, lo que proporciona una visión aún más detallada de tu microbioma. Las subespecies pueden variar entre personas y ofrecer información específica sobre la composición y variaciones en tu microbioma cutáneo. %s", warning_virus)
  perfil_taxonomico <- "5.3. Perfil taxonómico"
  subtitulo_biomarcador_rango <- "6. Envejecimiento de la piel"
  titulotipopiel <- "7. Rasgos de la piel"
  titulo_riesgo <- "8. Factores de riesgo"
  titulo_habitos <- "9. Factores del medio ambiente y estilo de vida que afectan tu piel"
  titulo_fagos <- "9. Factores beneficiosos"

  titulo_devolucion <- "10. Devolución personalizada"

  texto_pre_recomendacion <- "A continuación te ofrecemos algunas sugerencias personalizadas para mantener y mejorar el equilibrio de tu microbioma. Estas recomendaciones están diseñadas para promover un microbioma equilibrado y, por lo tanto, una piel más saludable."
  titulo_recomendacion <- "11. Recomendaciones"
  texto_recomendacion <- "Estamos trabajando en el desarrollo de cremas personalizadas para ayudarte a equilibrar tu microbioma. Mientras tanto para mejorar tu piel, a continuación, te dejamos algunas sugerencias:"

  #Armado de recomendaciones personalizadas:
  #Armar recomendaciones: ------------------------------------------------------------------------
  out <- reportCountsTax(id = id, MetadataB = MetadataB, patients_dir = patients_dir, de_host = "Biota", conEukaryota = FALSE)
  dev_cantidades <- out[[4]]
  out <- generateScore(patient_dir = patient_dir)
  dev_score_gral <- out[[2]]
  dev_scores <- paste(dev_cantidades, dev_score_gral, sep ="")

  out <- indicadorTipoPiel(id = id, MetadataB = MetadataB)
  dev_piel <- out[[4]]

  out <- indicadorRangoEtario(id = id, MetadataB = MetadataB)
  dev_edad <- out[[7]]

  Recomendaciones <- as.data.frame(read_excel("~/Daniela/Biota/PipelineBiota/data/topTax-Recomendaciones.xlsx", sheet=3))

  mensaje_recomendaciones1 <- ""
  mensaje_recomendaciones3 <- ""

  if(grepl("desequilibrado", dev_scores) |  grepl("menor", dev_scores)) {
    rec_deseq <- Recomendaciones$Recomendación[which(grepl("Índices desbalanceados", Recomendaciones$Indicador))]
    rec_subtitulo <- strsplit(rec_deseq, split= ":")[[1]][1]
    rec_texto <- strsplit(rec_deseq, split= ":")[[1]][2]
    mensaje_recomendaciones1 <- sprintf("__Descripción de la composición de tu microbioma - Índices desbalanceados:__

  - __%s__: %s", rec_subtitulo, rec_texto)
  }

  tipos_piel <- c("Piel Mixta", "Piel Seca", "Piel Grasa")
  for(tipo in tipos_piel) {
    print(tipo)
    if(grepl(tipo, dev_piel)) {
      rec_deseq <- Recomendaciones$Recomendación[which(grepl(tipo, Recomendaciones$Indicador))]
      rec_subtitulo1 <- strsplit(rec_deseq, split= ":")[[1]][1]
      rec_texto1 <- strsplit(rec_deseq, split= ":")[[1]][2]
      rec_subtitulo2 <- strsplit(rec_deseq, split= ":")[[2]][1]
      rec_texto2 <- strsplit(rec_deseq, split= ":")[[2]][2]
      mensaje_recomendaciones2 <- sprintf("__Rasgos de la piel - %s:__

- __%s:__ %s

- __%s:__ %s", tipo, rec_subtitulo1, rec_texto1, rec_subtitulo2, rec_texto2)
    }
  }

  if(grepl("envejecida", dev_edad)) {
    rec_deseq <- Recomendaciones$Recomendación[which(grepl("Envejecida", Recomendaciones$Indicador))]
    rec_subtitulo1 <- strsplit(rec_deseq, split= ":")[[1]][1]
    rec_texto1 <- strsplit(rec_deseq, split= ":")[[1]][2]
    rec_subtitulo2 <- strsplit(rec_deseq, split= ":")[[2]][1]
    rec_texto2 <- strsplit(rec_deseq, split= ":")[[2]][2]
    mensaje_recomendaciones3 <- sprintf("__Envejecimiento de la piel - Envejecida:__
- __%s__:%s

- __%s__:%s", rec_subtitulo1, rec_texto1, rec_subtitulo2, rec_texto2)
  }

  recomendaciones <- c(mensaje_recomendaciones1, mensaje_recomendaciones2, mensaje_recomendaciones3)
  #recomendaciones <- paste(Filter(nzchar, mensajes), collapse = "\n\n")


  #---------------------------------------------------------

  titulo_final <- "12. Comentarios finales"
  texto_final <- "¡Gracias por haber participado en el estudio de Biotalife Skin!


  Tu colaboración nos ayuda a comprender el microbioma de la piel.

  Seguimos trabajando en desarrollar soluciones personalizadas para tu bienestar. Mantenete atenta/o a nuestras redes sociales para conocer los avances. Continúa cuidando de tu piel y recordá que cada pequeño cambio puede marcar una gran diferencia en tu salud."

  titulo_contacto <- "13. Contactanos"
  texto_contacto <- "Si querés conocer más sobre los resultados de tu estudio, te recomendamos consultar a un dermatólogo o podés agendar una sesión con alguno de nuestros especialistas para una orientación personalizada.

  contacto@biotalifeskin.com

  www.biotalifeskin.com


  Redes sociales:

  https://www.linkedin.com/company/biotalife/

  https://www.instagram.com/biotalife"

  espacio_vacio <- "    "

  # Contenido del archivo Rmd
  rmd_content <- sprintf('
---
title: "INFORME DE LA COMPOSICIÓN DEL MICROBIOMA DE LA PIEL"
author: "Muestra %s"
output:
  pdf_document:
    includes:
      in_header: "~/Daniela/Biota/PipelineBiota/data/header.tex"
  html_document:
    df_print: paged
    includes:
      in_header: "~/Daniela/Biota/PipelineBiota/data/header.html"
editor_options:
  markdown:
    wrap: 72
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)

library(readxl)
library(kableExtra)
library(dplyr)
library(tidyr)
library(png)

InfoGral <- info

print(reportMasAbundantes(id = id, MetadataB = MetadataB))

patient_dir <- "%s"

out <- reportCountsTax(id = id, MetadataB = MetadataB, patients_dir = patients_dir, de_host = "Biota", conEukaryota = FALSE)
CantidadesTax <- read_excel(sprintf("%s/Tabla_CantidadesTax.xlsx", patient_dir))
png_cantidades <- sprintf("%s/cantidades_barra.png", patient_dir)
mensaje_cantidades <- out[[3]]
dev_cantidades <- out[[4]]

Especies_masAbundantes <- read_excel(sprintf("%s/Tabla_Species_masAbundantes.xlsx", patient_dir))
Genus_masAbundantes <- read_excel(sprintf("%s/Tabla_Genus_masAbundantes.xlsx", patient_dir))
Phylum_masAbundantes <- read_excel(sprintf("%s/Tabla_Phylum_masAbundantes.xlsx", patient_dir))

SubEspecies_RangoEtario <- read_excel(sprintf("%s/Tabla_SubEspecies_RangoEtario.xlsx", patient_dir))

out <- reportMasAbundantes(id = id, MetadataB = MetadataB)
grid_plot <- out[[2]]
especies_paraCompleto <- out[[3]]
out <- out[[1]]
species_plot <- out[[1]]
genus_plot <- out[[2]]
phylum_plot <- out[[3]]

out <- generateScore(patient_dir = patient_dir)
mensaje_score <- out[[1]]
dev_score_gral <- out[[2]]
png_score <- (sprintf("%s/score_barra.png", patient_dir))

dev_scores <- paste(dev_cantidades, dev_score_gral)

png_barra <- sprintf("%s/rango_etario_barra.png", patient_dir)
png_biomarcadores <- (sprintf("%s/PieChart_BiomarcadoresRangoEtario.png", patient_dir))
out <- indicadorRangoEtario(id = id, MetadataB = MetadataB)
rango_predicho <- out[[2]]
plot_rango_real <- out[[3]]
combined_plot <- out[[4]]
mensaje <- out[[5]]
barra_plot <- out[[6]]
dev_edad <- out[[7]]

out <- indicadorTipoPiel(id = id, MetadataB = MetadataB)
mensaje_piel <- out[[2]]
dev_piel <- out[[4]]
png_tipopiel <- sprintf("%s/tipodepiel_barra.png", patient_dir)


out <- indicadorVirus(id = id, MetadataB = MetadataB)
mensaje_virus <- out[[1]]
dev_virus <- out[[2]]
png_virus <- sprintf("%s/virus_barra.png", patient_dir)
virus_paraCompleto <- out[[3]]
tabla_virus <- out[[4]]

out <- indicadorFagos(id = id, MetadataB = MetadataB)
mensaje_fagos <- out[[1]]
png_fagos <- sprintf("%s/fagos_barra.png", patient_dir)

png_niveles <- "~/Daniela/Biota/PipelineBiota/data/Niveles_tax2.png"



```

```{r, echo=FALSE, fig.align="center"}

cat("\n\n")

cat("\n\n")
knitr::kable(InfoGral)
cat("\n\n")
knitr::asis_output(paste0("### ", espacio_vacio, "\n\n"))
knitr::asis_output(paste0("### ", espacio_vacio, "\n\n"))

cat("\n\n")

knitr::asis_output(paste0("# ", titulo_titulos, "\n\n"))
cat("\n\n")
knitr::asis_output(paste0(titulo1, "\n\n"))
knitr::asis_output(paste0(titulo2, "\n\n"))
knitr::asis_output(paste0(titulo3, "\n\n"))
knitr::asis_output(paste0(titulo4, "\n\n"))
knitr::asis_output(paste0(titulo5, "\n\n"))
knitr::asis_output(paste0(subtitulo_biomarcador_rango, "\n\n"))
knitr::asis_output(paste0(titulotipopiel, "\n\n"))

knitr::asis_output(paste0(titulo_riesgo, "\n\n"))
knitr::asis_output(paste0(titulo_fagos, "\n\n"))
knitr::asis_output(paste0(titulo_devolucion, "\n\n"))
knitr::asis_output(paste0(titulo_recomendacion, "\n\n"))
knitr::asis_output(paste0(titulo_final, "\n\n"))
knitr::asis_output(paste0(titulo_contacto, "\n\n"))

cat("\n\n")
knitr::asis_output(paste0("# ", espacio_vacio, "\n\n"))
knitr::asis_output(paste0("# ", espacio_vacio, "\n\n"))
knitr::asis_output(paste0("# ", espacio_vacio, "\n\n"))

cat("\n\n")
knitr::asis_output(paste0(" ", "\n\n"))
knitr::asis_output(paste0("### ", titulo1, "\n\n"))
knitr::asis_output(paste0(texto1, "\n\n"))

knitr::asis_output(paste0("### ", titulo2, "\n\n"))
knitr::asis_output(paste0(texto2, "\n\n"))
cat("\n\n")
knitr::include_graphics(png_niveles)
cat("\n\n")

knitr::asis_output(paste0("### ", titulo3, "\n\n"))
knitr::asis_output(paste0(texto3, "\n\n"))

knitr::asis_output(paste0("### ", titulo4, "\n\n"))
knitr::asis_output(paste0(texto4, "\n\n"))

knitr::asis_output(paste0("### ", titulo5, "\n\n"))
knitr::asis_output(paste0(texto5, "\n\n"))

cat("\n\n")

knitr::asis_output(paste0("### ", "5.1. Índice de variabilidad", "\n\n"))
knitr::include_graphics(png_cantidades)
cat("\n\n")
knitr::asis_output(paste0(mensaje_cantidades, "\n\n"))
cat("\n\n")

knitr::asis_output(paste0("### ", "5.2. Índice de cantidad", "\n\n"))
cat("\n\n")
knitr::include_graphics(png_score)
cat("\n\n")
knitr::asis_output(paste0(mensaje_score, "\n\n"))
cat("\n\n")

knitr::asis_output(paste0("### ", "5.3. Perfil taxonómico", "\n\n"))
knitr::asis_output(paste0(texto_final_score, "\n\n"))
cat("\n\n")

knitr::asis_output(paste0(texto_filos, "\n\n"))
cat("\n\n")

knitr::asis_output(paste0(texto_generos, "\n\n"))
cat("\n\n")

knitr::asis_output(paste0(texto_especies, "\n\n"))
cat("\n\n")

knitr::asis_output(paste0("En la sección de recomendaciones se sugieren cambios en los hábitos que pueden colaborar a mejorar la cantidad y el equilibrio de los microorganismos en tu piel.", "\n\n", "\n\n"))
cat("\n\n")
cat("\n\n")
cat("\n\n")
plot(grid_plot)
cat("\n\n")
knitr::asis_output(paste0("Gráfico que indica la cantidad de Filos, Géneros y Especies más abundantes de tu muestra en comparación con los valores equilibrados de la plataforma de Biotalife Skin.", "\n\n"))
cat("\n\n")
knitr::kable(Filos_masAbundantes)
cat("\n\n")

knitr::kable(Genus_masAbundantes)
cat("\n\n")
knitr::asis_output(paste0("### ", espacio_vacio, "\n\n"))
knitr::asis_output(paste0("### ", espacio_vacio, "\n\n"))
cat("\n\n")
knitr::kable(esp_masAbundantes)
cat("\n\n")
knitr::asis_output(paste0(especies_paraCompleto, "\n\n"))
cat("\n\n")

knitr::asis_output(paste0("### ", subtitulo_biomarcador_rango, "\n\n"))
cat("\n\n")
knitr::include_graphics(png_barra)
cat("\n\n")
knitr::asis_output(paste0(mensaje, "\n\n"))
cat("\n\n")

knitr::asis_output(paste0("### ", titulotipopiel, "\n\n"))
cat("\n\n")
knitr::include_graphics(png_tipopiel)
cat("\n\n")
knitr::asis_output(paste0(mensaje_piel, "\n\n"))
cat("\n\n")


knitr::asis_output(paste0("### ", titulo_riesgo, "\n\n"))
cat("\n\n")
knitr::include_graphics(png_virus)
cat("\n\n")
knitr::asis_output(paste0(mensaje_virus, "\n\n"))
cat("\n\n")
knitr::asis_output(paste0(virus_paraCompleto, "\n\n"))
cat("\n\n")


knitr::asis_output(paste0("### ", titulo_fagos, "\n\n"))
cat("\n\n")
knitr::include_graphics(png_fagos)
cat("\n\n")
knitr::asis_output(paste0(mensaje_fagos, "\n\n"))
cat("\n\n")

knitr::asis_output(paste0("### ", titulo_devolucion, "\n\n"))
cat("\n\n")
knitr::asis_output(paste0(dev_scores, "\n\n"))
cat("\n\n")
knitr::asis_output(paste0(dev_edad, dev_piel, dev_virus, "\n\n"))
cat("\n\n")
cat("\n\n")
cat("\n\n")
knitr::asis_output(paste0(texto_pre_recomendacion, "\n\n"))
cat("\n\n")
knitr::asis_output(paste0("### ", titulo_recomendacion, "\n\n"))
cat("\n\n")
knitr::asis_output(paste0(texto_recomendacion, "\n\n"))
cat("\n\n")
knitr::asis_output(paste0(recomendaciones, "\n\n"))
cat("\n\n")
knitr::asis_output(paste0("### ", titulo_final, "\n\n"))
cat("\n\n")
knitr::asis_output(paste0("#### ",texto_final, "\n\n"))
cat("\n\n")

knitr::asis_output(paste0("### ", titulo_contacto, "\n\n"))
cat("\n\n")
knitr::asis_output(paste0("### ",texto_contacto, "\n\n"))


```
',  id, patient_dir, patient_dir, patient_dir, patient_dir, patient_dir, patient_dir, patient_dir, patient_dir,patient_dir, patient_dir, patient_dir, patient_dir, patient_dir, patient_dir, patient_dir)

  rmd_file <- sprintf("%s/CompleteReportResults_Biota_ID%s.Rmd", patient_dir, id)
  writeLines(rmd_content, con = rmd_file)
  output_file <- sprintf("%s/CompleteReportResults_Biota_ID%s.pdf", patient_dir, id)
  rmarkdown::render(rmd_file, output_file = output_file)
  #output_file <- sprintf("%s/ReportResults_Biota.doc", patient_dir)
  #rmarkdown::render(rmd_file, output_file = output_file)

  #Elimino todo lo que se genera y no me sirve:
  file.remove(sprintf("%s/cantidades_barra.png", patient_dir))
  file.remove(sprintf("%s/score_barra.png", patient_dir))
  file.remove(sprintf("%s/rango_etario_barra.png", patient_dir))
  file.remove(sprintf("%s/fagos_barra.png", patient_dir))
  file.remove(sprintf("%s/virus_barra.png", patient_dir))
  file.remove(sprintf("%s/tipodepiel_barra.png", patient_dir))

  file.remove(sprintf("%s/Genus_masAbundantes.png", patient_dir))
  file.remove(sprintf("%s/Species_masAbundantes.png", patient_dir))
  file.remove(sprintf("%s/Phylum_masAbundantes.png", patient_dir))

  file.remove(sprintf("%s/Tabla_Genus_masAbundantes.xlsx", patient_dir))
  file.remove(sprintf("%s/Tabla_Phylum_masAbundantes.xlsx", patient_dir))
  file.remove(sprintf("%s/Tabla_Species_masAbundantes.xlsx", patient_dir))
  file.remove(sprintf("%s/Tabla_SubSpecies_masAbundantes.xlsx", patient_dir))
  file.remove(sprintf("%s/Tabla_SubEspecies_RangoEtario.xlsx", patient_dir))
  file.remove(sprintf("%s/Tabla_CantidadesTax.xlsx", patient_dir))



  file.remove(sprintf("%s/CompleteReportResults_Biota_ID%s.log", patient_dir, id))
  file.remove(sprintf("%s/CompleteReportResults_Biota_ID%s.Rmd", patient_dir, id))

  return(output_file)

}


#########################################################################################

caratulaComplete <- function(patient_dir) {

  id <- basename(patient_dir)
  patients_dir <- dirname(patient_dir)
  MetadataB <- as.data.frame(read_excel("/media/4tb2/Daniela/Biota/PipelineBiota-master/data/Metadata-soloColumnasUsables.xlsx"))
  #BIOTALIFE_SKIN_Respuestas_ <- read_excel("~/Daniela/Biota/BIOTALIFE SKIN  (Respuestas).xlsx", sheet = 4)
  path_metadata <- "/media/4tb2/Daniela/Biota/PipelineBiota-master/data/Metadata-Completa-SinLimpiar.xlsx"
  BIOTALIFE_SKIN_Respuestas_ <- read_excel(path_metadata)
  
  #Obtener info gral
  BIOTALIFE_SKIN_Respuestas_$ID  <- gsub("\\.0$", "", as.character(BIOTALIFE_SKIN_Respuestas_$ID))
  BIOTALIFE_SKIN_Respuestas_$ID[which(BIOTALIFE_SKIN_Respuestas_$ID == "118")] <- "118-1"
  BIOTALIFE_SKIN_Respuestas_$ID[which(BIOTALIFE_SKIN_Respuestas_$ID == "184")] <- "184-1"

  info <- BIOTALIFE_SKIN_Respuestas_[which(BIOTALIFE_SKIN_Respuestas_$ID == id),]
  colnames(info)
  info <- info[, which(colnames(info) %in% c("ID", "1. Nombre y Apellido","3. Sexo", "2. Fecha de Nacimiento", "Edad", "Fecha Cita" ))]
  info <- as.data.frame(t(info))
  info <- data.frame("Característica" = rownames(info), "Información" = info$V1)
  info$Característica <- c("ID", "Nombre y Apellido", "Fecha de nacimiento", "Sexo", "Fecha Cita", "Edad")

  nombre <- info$Información[which(info$Característica == "Nombre y Apellido")]
  capitalizar_apellido <- function(nombre) {
    partes <- strsplit(nombre, " ")[[1]]  # Divide el nombre en partes
    partes[2] <- paste0(toupper(substring(partes[2], 1, 1)), substring(partes[2], 2))  # Capitaliza la primera letra del apellido
    return(paste(partes, collapse = " "))  # Vuelve a unir las partes
  }
  nombre_modificado <- capitalizar_apellido(nombre)

  info$Información[which(info$Característica == "Nombre y Apellido")] <- nombre_modificado

  #Agrego info del paciente si es que NO esta en MetadabaB pero SI en BiotaLife respuestas: --------------------------
  MetadataB <- MetadataB[, -which(colnames(MetadataB) == "TipoPielViejo")]
  MetadataB <- MetadataB[, -which(colnames(MetadataB) == "Secuenciado")]
  colnames(BIOTALIFE_SKIN_Respuestas_)[c(2,5,6,7,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,25,26,27,28,29,30,31,33,34,57, 60, 61, 32,35)] <- c("ID", "FechadeNacimiento", "Sexo", "Email",
                                                                                                                                                 "ColorCabello", "ColorOjos", "FacilidadBroncearse", "Peso",
                                                                                                                                                 "Altura", "AntecedenteEnfermedad", "CualEnfermedad", "AfeccionesPiel",
                                                                                                                                                 "OtraAfeccionPiel", "TratamientoEstético3meses", "MétodoAnticonceptivo", "EmbarazadaoAmantando",
                                                                                                                                                 "TratamientoMédico", "CualTratamientoMedico", "FrecuenciaTratamientoMedico", "Tabaco",
                                                                                                                                                 "CigarrillosporSemana", "Alcohol", "AlcoholporSemana", "ActividadFísica",
                                                                                                                                                 "CuándoLimpiaCara", "ConquéLimpiaCara", "AplicacionProtectorSolar", "CuándoMaquillaje",
                                                                                                                                                 "Fecha Cita", "Edad", "Rango etario",
                                                                                                                                                 "Tipodepiel", "Maquillaje_Base")

  if(!(id %in% MetadataB$ID) & id %in% BIOTALIFE_SKIN_Respuestas_$ID) {
    MetadataB <- rbind(MetadataB, BIOTALIFE_SKIN_Respuestas_[which(BIOTALIFE_SKIN_Respuestas_$ID == id), which(colnames(BIOTALIFE_SKIN_Respuestas_) %in% colnames(MetadataB))])
    #modifico lo de maquillaje base y lo de facilidad de bronceado
    MetadataB$FacilidadBroncearse <- sub("^(\\d+).*", "\\1", MetadataB$FacilidadBroncearse)
    MetadataB$FacilidadBroncearse[which(MetadataB$FacilidadBroncearse == "1")] <- "2"
    MetadataB$CuándoMaquillaje[which(MetadataB$CuándoMaquillaje == "1 o 2 veces por semana")] <- "Diariamente"
    MetadataB$AplicacionProtectorSolar[which(MetadataB$AplicacionProtectorSolar == "1 o 2 veces por semana")] <- "Solo ante exposición en verano"

    if(any(MetadataB$Maquillaje_Base != "0" | MetadataB$Maquillaje_Base != "1")) {
      maq_base <-  MetadataB$Maquillaje_Base[-which(MetadataB$Maquillaje_Base == "0" | MetadataB$Maquillaje_Base == "1")]
      MetadataB$Maquillaje_Base[which(MetadataB$Maquillaje_Base == maq_base)] <- ifelse(grepl("base", maq_base) | grepl("corrector", maq_base) | grepl("polvo", maq_base), "1", "0")
    }
    MetadataB$Maquillaje_Base[which(is.na(MetadataB$Maquillaje_Base))] <- "0"
  }
  #----------------------------------------------------------------------------------------------------------


  # Contenido del archivo Rmd
  rmd_content <- sprintf('
---
title: "INFORME DE LA COMPOSICIÓN DEL MICROBIOMA DE LA PIEL"
author: "Muestra %s"
output:
  pdf_document:
    includes:
      in_header: "/media/4tb2/Daniela/Biota/PipelineBiota-master/data/header.tex"
  html_document:
    df_print: paged
    includes:
      in_header: "/media/4tb2/Daniela/Biota/PipelineBiota-master/data/header.html"
editor_options:
  markdown:
    wrap: 72
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)

library(readxl)
library(kableExtra)
library(dplyr)
library(tidyr)
library(png)

InfoGral <- info

```

```{r, echo=FALSE, fig.align="center"}

cat("\n\n")

cat("\n\n")
knitr::kable(InfoGral)
cat("\n\n")


```
',  id, patient_dir, patient_dir, patient_dir, patient_dir, patient_dir, patient_dir, patient_dir, patient_dir,patient_dir, patient_dir, patient_dir, patient_dir, patient_dir, patient_dir, patient_dir)

  rmd_file <- sprintf("%s/CompleteCaratula_ID%s.Rmd", patient_dir, id_mr)
  writeLines(rmd_content, con = rmd_file)
  output_file <- sprintf("%s/CompleteCaratula_ID%s.pdf", patient_dir, id_mr)
  rmarkdown::render(rmd_file, output_file = output_file)

  return(output_file)
}

################################################################################

#################################################################################################

contenidoComplete <- function(patient_dir) {

  id <- basename(patient_dir)
  patients_dir <- dirname(patient_dir)
  MetadataB <- as.data.frame(read_excel("/media/4tb2/Daniela/Biota/PipelineBiota-master/data/Metadata-soloColumnasUsables.xlsx"))
  #BIOTALIFE_SKIN_Respuestas_ <- read_excel("~/Daniela/Biota/BIOTALIFE SKIN  (Respuestas).xlsx", sheet = 4)
  path_metadata <- "/media/4tb2/Daniela/Biota/PipelineBiota-master/data/Metadata-Completa-SinLimpiar.xlsx"
  BIOTALIFE_SKIN_Respuestas_ <- read_excel(path_metadata)
  
  #Obtener info gral
  BIOTALIFE_SKIN_Respuestas_$ID  <- gsub("\\.0$", "", as.character(BIOTALIFE_SKIN_Respuestas_$ID))
  BIOTALIFE_SKIN_Respuestas_$ID[which(BIOTALIFE_SKIN_Respuestas_$ID == "118")] <- "118-1"
  BIOTALIFE_SKIN_Respuestas_$ID[which(BIOTALIFE_SKIN_Respuestas_$ID == "184")] <- "184-1"

  info <- BIOTALIFE_SKIN_Respuestas_[which(BIOTALIFE_SKIN_Respuestas_$ID == id),]
  colnames(info)
  info <- info[, which(colnames(info) %in% c("ID", "1. Nombre y Apellido","3. Sexo", "2. Fecha de Nacimiento", "Edad", "Fecha Cita" ))]
  info <- as.data.frame(t(info))
  info <- data.frame("Característica" = rownames(info), "Información" = info$V1)
  info$Característica <- c("ID", "Nombre y Apellido", "Fecha de nacimiento", "Sexo", "Fecha Cita", "Edad")

  nombre <- info$Información[which(info$Característica == "Nombre y Apellido")]
  capitalizar_apellido <- function(nombre) {
    partes <- strsplit(nombre, " ")[[1]]  # Divide el nombre en partes
    partes[2] <- paste0(toupper(substring(partes[2], 1, 1)), substring(partes[2], 2))  # Capitaliza la primera letra del apellido
    return(paste(partes, collapse = " "))  # Vuelve a unir las partes
  }
  nombre_modificado <- capitalizar_apellido(nombre)

  info$Información[which(info$Característica == "Nombre y Apellido")] <- nombre_modificado

  #Agrego info del paciente si es que NO esta en MetadabaB pero SI en BiotaLife respuestas: --------------------------
  #MetadataB <- MetadataB[, -which(colnames(MetadataB) == "TipoPielViejo")]
  #MetadataB <- MetadataB[, -which(colnames(MetadataB) == "Secuenciado")]
  colnames(BIOTALIFE_SKIN_Respuestas_)[c(2,5,6,7,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,25,26,27,28,29,30,31,33,34,57, 60, 61, 32,35)] <- c("ID", "FechadeNacimiento", "Sexo", "Email",
                                                                                                                                                 "ColorCabello", "ColorOjos", "FacilidadBroncearse", "Peso",
                                                                                                                                                 "Altura", "AntecedenteEnfermedad", "CualEnfermedad", "AfeccionesPiel",
                                                                                                                                                 "OtraAfeccionPiel", "TratamientoEstético3meses", "MétodoAnticonceptivo", "EmbarazadaoAmantando",
                                                                                                                                                 "TratamientoMédico", "CualTratamientoMedico", "FrecuenciaTratamientoMedico", "Tabaco",
                                                                                                                                                 "CigarrillosporSemana", "Alcohol", "AlcoholporSemana", "ActividadFísica",
                                                                                                                                                 "CuándoLimpiaCara", "ConquéLimpiaCara", "AplicacionProtectorSolar", "CuándoMaquillaje",
                                                                                                                                                 "Fecha Cita", "Edad", "Rango etario",
                                                                                                                                                 "Tipodepiel", "Maquillaje_Base")

  if(!(id %in% MetadataB$ID) & id %in% BIOTALIFE_SKIN_Respuestas_$ID) {
    MetadataB <- rbind(MetadataB, BIOTALIFE_SKIN_Respuestas_[which(BIOTALIFE_SKIN_Respuestas_$ID == id), which(colnames(BIOTALIFE_SKIN_Respuestas_) %in% colnames(MetadataB))])
    #modifico lo de maquillaje base y lo de facilidad de bronceado
    MetadataB$FacilidadBroncearse <- sub("^(\\d+).*", "\\1", MetadataB$FacilidadBroncearse)
    MetadataB$FacilidadBroncearse[which(MetadataB$FacilidadBroncearse == "1")] <- "2"
    MetadataB$CuándoMaquillaje[which(MetadataB$CuándoMaquillaje == "1 o 2 veces por semana")] <- "Diariamente"
    MetadataB$AplicacionProtectorSolar[which(MetadataB$AplicacionProtectorSolar == "1 o 2 veces por semana")] <- "Solo ante exposición en verano"

    if(any(MetadataB$Maquillaje_Base != "0" | MetadataB$Maquillaje_Base != "1")) {
      maq_base <-  MetadataB$Maquillaje_Base[-which(MetadataB$Maquillaje_Base == "0" | MetadataB$Maquillaje_Base == "1")]
      MetadataB$Maquillaje_Base[which(MetadataB$Maquillaje_Base == maq_base)] <- ifelse(grepl("base", maq_base) | grepl("corrector", maq_base) | grepl("polvo", maq_base), "1", "0")
    }
    MetadataB$Maquillaje_Base[which(is.na(MetadataB$Maquillaje_Base))] <- "0"
  }
  #----------------------------------------------------------------------------------------------------------


  # Obtener datos de las tablas
  counts_tax <- counts_Tax(patients_dir = patient_dir, source = "KRAKEN", de_host = "Bowtie", conEukaryota = FALSE)

  reportCountsTax(id=id, MetadataB = MetadataB, patients_dir = patients_dir, de_host = "Bowtie", conEukaryota = FALSE)
  reportMasAbundantes(patient_dir = patient_dir, MetadataB = MetadataB)
  indicadorRangoEtario(id = id, MetadataB = MetadataB)
  titulo_titulos <- "Qué podés encontrar en este informe:"
  titulo1 <- "1. ¿Qué es el microbioma y por qué es importante para la salud de tu piel?"
  texto1 <- "El microbioma de la piel es una comunidad compleja de microorganismos, como bacterias, hongos y virus, que viven en su superficie. Estos microorganismos cumplen funciones importantes como: proteger la piel de agentes dañinos, mantener la hidratación y el equilibrio del pH, y fortalecer el sistema inmunológico. Factores como la edad, el sexo, la dieta y el estilo de vida pueden influir en la composición del microbioma. Mantenerlo equilibrado es esencial para una piel sana, ya que su desequilibrio puede contribuir a problemas como el acné, la dermatitis atópica, la rosácea y el envejecimiento prematuro."

  titulo2 <- " 2. ¿Cómo interpretar tus resultados?"
  texto2 <- "En este informe, se analiza la cantidad de microorganismos que componen el microbioma de tu piel. Se evalúan en diferentes niveles taxonómicos (filo, clase, orden, familia, género y especie) y se clasifica según su presencia en comparación con una base de datos de referencia de piel saludable de Biotalife Skin. Es importante destacar que un valor fuera del rango no significa necesariamente un problema de salud, ya que el microbioma puede variar por factores temporales como cambios en el entorno, uso de productos cosméticos o consumo de medicamentos."

  titulo3 <- "3. ¿En qué consiste nuestro estudio?"
  texto3 <- "Este estudio utiliza tecnología de secuenciación de nueva generación (NGS) de Illumina para analizar el ADN de los microorganismos presentes en una muestra de la piel obtenida mediante un hisopado estéril. La secuenciación permite identificar con precisión los diferentes microorganismos, brindando una visión detallada de la composición del microbioma de la piel."

  titulo4 <- "4. Descargo de responsabilidad"
  texto4 <- "Cada parámetro evaluado se clasifica como ‘Bajo’, ‘Equilibrado’ o ‘Alto’ utilizando un algoritmo desarrollado por nuestro equipo. Además, identificamos ciertas especies que pueden actuar como indicadores de la edad de la piel, lo que nos permite predecir si el microbioma corresponde a tu edad biológica. Es importante tener en cuenta que estos resultados son informativos, no es un diagnóstico médico ni reemplaza a una consulta médica. "

  titulo5 <- "5. Descripción de la composición de tu microbioma"
  texto5 <- "Tu microbioma es único, como tu huella digital. En esta sección, detallaremos la variabilidad y cantidad de microorganismos, destacando los grupos más abundantes.  Se mostrarán índices que van del 0 al 1, presentando resultados más altos en microbiomas más equilibrados, lo que contribuye a mantener una piel saludable."
  titulo_score <- "Índice de microorganismos más abundantes"
  texto_final_score <- "A continuación, se detallan los microorganismos más abundantes a nivel de Filo, Género y Especie (ver gráfico):"

  subtitulo_filos <- "Filos"
  texto_filos <- "Los FILOS bacterianos más comunes en la piel son Actinobacteria,  Firmicutes, Proteobacteria y Bacteroidetes."
  Filos_masAbundantes <- read_excel(sprintf("%s/Tabla_Phylum_masAbundantes.xlsx", patient_dir))

  subtitulo_generos <- "Géneros"
  Genus_masAbundantes <- read_excel(sprintf("%s/Tabla_Genus_masAbundantes.xlsx", patient_dir))
  generos_top3 <- paste0(Genus_masAbundantes$Géneros[1:3], collapse = ", ")
  texto_generos <- sprintf("Los GÉNEROS más abundantes de tu piel son %s.", generos_top3)


  subtitulo_especies <- "Especies"
  esp_masAbundantes <- read_excel(sprintf("%s/Tabla_Species_masAbundantes.xlsx", patient_dir))
  especies_top3 <- paste0(esp_masAbundantes$Especies[1:3], collapse = ", ")

  texto_especies <- sprintf("Las ESPECIES más abundantes de tu piel son %s.", especies_top3)


  #texto_subespecies <- sprintf("En algunos casos, se pueden identificar subespecies, lo que proporciona una visión aún más detallada de tu microbioma. Las subespecies pueden variar entre personas y ofrecer información específica sobre la composición y variaciones en tu microbioma cutáneo. %s", warning_virus)
  perfil_taxonomico <- "5.3. Perfil taxonómico"
  subtitulo_biomarcador_rango <- "6. Envejecimiento de la piel"
  titulotipopiel <- "7. Rasgos de la piel"
  titulo_riesgo <- "8. Factores de riesgo"
  titulo_habitos <- "9. Factores del medio ambiente y estilo de vida que afectan tu piel"
  titulo_fagos <- "9. Factores beneficiosos"

  titulo_devolucion <- "10. Devolución personalizada"

  texto_pre_recomendacion <- "A continuación te ofrecemos algunas sugerencias personalizadas para mantener y mejorar el equilibrio de tu microbioma. Estas recomendaciones están diseñadas para promover un microbioma equilibrado y, por lo tanto, una piel más saludable."
  titulo_recomendacion <- "11. Recomendaciones"
  texto_recomendacion <- "Estamos trabajando en el desarrollo de cremas personalizadas para ayudarte a equilibrar tu microbioma. Mientras tanto para mejorar tu piel, a continuación, te dejamos algunas sugerencias:"

  #Armado de recomendaciones personalizadas:
  #Armar recomendaciones: ------------------------------------------------------------------------
  out <- reportCountsTax(id = id, MetadataB = MetadataB, patients_dir = patients_dir, de_host = "Bowtie", conEukaryota = FALSE)
  dev_cantidades <- out[[4]]
  out <- generateScore(patient_dir = patient_dir)
  dev_score_gral <- out[[2]]
  dev_scores <- paste(dev_cantidades, dev_score_gral, sep ="")

  out <- indicadorTipoPiel(patient_dir = patient_dir, MetadataB = MetadataB)
  dev_piel <- out[[4]]

  out <- indicadorRangoEtario(id = id, MetadataB = MetadataB)
  dev_edad <- out[[7]]

  Recomendaciones <- as.data.frame(read_excel("/media/4tb2/Daniela/Biota/PipelineBiota-master/data/topTax-Recomendaciones.xlsx", sheet=3))

  mensaje_recomendaciones1 <- ""
  mensaje_recomendaciones3 <- ""

  if(grepl("desequilibrado", dev_scores) |  grepl("menor", dev_scores)) {
    rec_deseq <- Recomendaciones$Recomendación[which(grepl("Índices desbalanceados", Recomendaciones$Indicador))]
    rec_subtitulo <- strsplit(rec_deseq, split= ":")[[1]][1]
    rec_texto <- strsplit(rec_deseq, split= ":")[[1]][2]
    mensaje_recomendaciones1 <- sprintf("__Descripción de la composición de tu microbioma - Índices desbalanceados:__

  - __%s__: %s", rec_subtitulo, rec_texto)
  }

  tipos_piel <- c("Piel Mixta", "Piel Seca", "Piel Grasa")
  for(tipo in tipos_piel) {
    print(tipo)
    if(grepl(tipo, dev_piel)) {
      rec_deseq <- Recomendaciones$Recomendación[which(grepl(tipo, Recomendaciones$Indicador))]
      rec_subtitulo1 <- strsplit(rec_deseq, split= ":")[[1]][1]
      rec_texto1 <- strsplit(rec_deseq, split= ":")[[1]][2]
      rec_subtitulo2 <- strsplit(rec_deseq, split= ":")[[2]][1]
      rec_texto2 <- strsplit(rec_deseq, split= ":")[[2]][2]
      mensaje_recomendaciones2 <- sprintf("__Rasgos de la piel - %s:__

- __%s:__ %s

- __%s:__ %s", tipo, rec_subtitulo1, rec_texto1, rec_subtitulo2, rec_texto2)
    }
  }

  if(grepl("envejecida", dev_edad)) {
    rec_deseq <- Recomendaciones$Recomendación[which(grepl("Envejecida", Recomendaciones$Indicador))]
    rec_subtitulo1 <- strsplit(rec_deseq, split= ":")[[1]][1]
    rec_texto1 <- strsplit(rec_deseq, split= ":")[[1]][2]
    rec_subtitulo2 <- strsplit(rec_deseq, split= ":")[[2]][1]
    rec_texto2 <- strsplit(rec_deseq, split= ":")[[2]][2]
    mensaje_recomendaciones3 <- sprintf("__Envejecimiento de la piel - Envejecida:__
- __%s__:%s

- __%s__:%s", rec_subtitulo1, rec_texto1, rec_subtitulo2, rec_texto2)
  }

  recomendaciones <- c(mensaje_recomendaciones1, mensaje_recomendaciones2, mensaje_recomendaciones3)
  #recomendaciones <- paste(Filter(nzchar, mensajes), collapse = "\n\n")


  #---------------------------------------------------------

  titulo_final <- "12. Comentarios finales"
  texto_final <- "¡Gracias por haber participado en el estudio de Biotalife Skin!


  Tu colaboración nos ayuda a comprender el microbioma de la piel.

  Seguimos trabajando en desarrollar soluciones personalizadas para tu bienestar. Mantenete atenta/o a nuestras redes sociales para conocer los avances. Continúa cuidando de tu piel y recordá que cada pequeño cambio puede marcar una gran diferencia en tu salud."

  titulo_contacto <- "13. Contactanos"
  texto_contacto <- "Si querés conocer más sobre los resultados de tu estudio, te recomendamos consultar a un dermatólogo o podés agendar una sesión con alguno de nuestros especialistas para una orientación personalizada.

  contacto@biotalifeskin.com

  www.biotalifeskin.com


  Redes sociales:

  https://www.linkedin.com/company/biotalife/

  https://www.instagram.com/biotalife"

  espacio_vacio <- "    "

  # Contenido del archivo Rmd
  rmd_content <- sprintf('


```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)

library(readxl)
library(kableExtra)
library(dplyr)
library(tidyr)
library(png)

InfoGral <- info

print(reportMasAbundantes(patient_dir = patient_dir, MetadataB = MetadataB))

patient_dir <- "%s"

out <- reportCountsTax(id = id, MetadataB = MetadataB, patients_dir = patients_dir, de_host = "Bowtie", conEukaryota = FALSE)
CantidadesTax <- read_excel(sprintf("%s/Tabla_CantidadesTax.xlsx", patient_dir))
png_cantidades <- sprintf("%s/cantidades_barra.png", patient_dir)
mensaje_cantidades <- out[[3]]
dev_cantidades <- out[[4]]

Especies_masAbundantes <- read_excel(sprintf("%s/Tabla_Species_masAbundantes.xlsx", patient_dir))
Genus_masAbundantes <- read_excel(sprintf("%s/Tabla_Genus_masAbundantes.xlsx", patient_dir))
Phylum_masAbundantes <- read_excel(sprintf("%s/Tabla_Phylum_masAbundantes.xlsx", patient_dir))

SubEspecies_RangoEtario <- read_excel(sprintf("%s/Tabla_SubEspecies_RangoEtario.xlsx", patient_dir))

out <- reportMasAbundantes(patient_dir = patient_dir, MetadataB = MetadataB)
grid_plot <- out[[2]]
especies_paraCompleto <- out[[3]]
out <- out[[1]]
species_plot <- out[[1]]
genus_plot <- out[[2]]
phylum_plot <- out[[3]]

out <- generateScore(patient_dir = patient_dir)
mensaje_score <- out[[1]]
dev_score_gral <- out[[2]]
png_score <- (sprintf("%s/score_barra.png", patient_dir))

dev_scores <- paste(dev_cantidades, dev_score_gral)

png_barra <- sprintf("%s/rango_etario_barra.png", patient_dir)
png_biomarcadores <- (sprintf("%s/PieChart_BiomarcadoresRangoEtario.png", patient_dir))
out <- indicadorRangoEtario(id = id, MetadataB = MetadataB)
rango_predicho <- out[[2]]
plot_rango_real <- out[[3]]
combined_plot <- out[[4]]
mensaje <- out[[5]]
barra_plot <- out[[6]]
dev_edad <- out[[7]]

out <- indicadorTipoPiel(patient_dir = patient_dir, MetadataB = MetadataB)
mensaje_piel <- out[[2]]
dev_piel <- out[[4]]
png_tipopiel <- sprintf("%s/tipodepiel_barra.png", patient_dir)


out <- indicadorVirus(id = id, MetadataB = MetadataB)
mensaje_virus <- out[[1]]
dev_virus <- out[[2]]
png_virus <- sprintf("%s/virus_barra.png", patient_dir)
virus_paraCompleto <- out[[3]]
tabla_virus <- out[[4]]

out <- indicadorFagos(id = id, MetadataB = MetadataB)
mensaje_fagos <- out[[1]]
png_fagos <- sprintf("%s/fagos_barra.png", patient_dir)

png_niveles <- "/media/4tb2/Daniela/Biota/PipelineBiota-master/data/Niveles_tax2.png"

```

```{r, echo=FALSE, fig.align="center"}

knitr::asis_output(paste0(" ", "\n\n"))
knitr::asis_output(paste0("### ", titulo1, "\n\n"))
knitr::asis_output(paste0(texto1, "\n\n"))

knitr::asis_output(paste0("### ", titulo2, "\n\n"))
knitr::asis_output(paste0(texto2, "\n\n"))
cat("\n\n")
knitr::include_graphics(png_niveles)
cat("\n\n")

knitr::asis_output(paste0("### ", titulo3, "\n\n"))
knitr::asis_output(paste0(texto3, "\n\n"))

knitr::asis_output(paste0("### ", titulo4, "\n\n"))
knitr::asis_output(paste0(texto4, "\n\n"))

knitr::asis_output(paste0("### ", titulo5, "\n\n"))
knitr::asis_output(paste0(texto5, "\n\n"))

cat("\n\n")

knitr::asis_output(paste0("### ", "5.1. Índice de variabilidad", "\n\n"))
knitr::include_graphics(png_cantidades)
cat("\n\n")
knitr::asis_output(paste0(mensaje_cantidades, "\n\n"))
cat("\n\n")

knitr::asis_output(paste0("### ", "5.2. Índice de cantidad", "\n\n"))
cat("\n\n")
knitr::include_graphics(png_score)
cat("\n\n")
knitr::asis_output(paste0(mensaje_score, "\n\n"))
cat("\n\n")

knitr::asis_output(paste0("### ", "5.3. Perfil taxonómico", "\n\n"))
knitr::asis_output(paste0(texto_final_score, "\n\n"))
cat("\n\n")

knitr::asis_output(paste0(texto_filos, "\n\n"))
cat("\n\n")

knitr::asis_output(paste0(texto_generos, "\n\n"))
cat("\n\n")

knitr::asis_output(paste0(texto_especies, "\n\n"))
cat("\n\n")

knitr::asis_output(paste0("En la sección de recomendaciones se sugieren cambios en los hábitos que pueden colaborar a mejorar la cantidad y el equilibrio de los microorganismos en tu piel.", "\n\n", "\n\n"))
cat("\n\n")
cat("\n\n")
cat("\n\n")
plot(grid_plot)
cat("\n\n")
knitr::asis_output(paste0("Gráfico que indica la cantidad de Filos, Géneros y Especies más abundantes de tu muestra en comparación con los valores equilibrados de la plataforma de Biotalife Skin.", "\n\n"))
cat("\n\n")
knitr::kable(Filos_masAbundantes)
cat("\n\n")

knitr::kable(Genus_masAbundantes)
cat("\n\n")
knitr::asis_output(paste0("### ", espacio_vacio, "\n\n"))
knitr::asis_output(paste0("### ", espacio_vacio, "\n\n"))
cat("\n\n")
knitr::kable(esp_masAbundantes)
cat("\n\n")
knitr::asis_output(paste0(especies_paraCompleto, "\n\n"))
cat("\n\n")

knitr::asis_output(paste0("### ", subtitulo_biomarcador_rango, "\n\n"))
cat("\n\n")
knitr::include_graphics(png_barra)
cat("\n\n")
knitr::asis_output(paste0(mensaje, "\n\n"))
cat("\n\n")

knitr::asis_output(paste0("### ", titulotipopiel, "\n\n"))
cat("\n\n")
knitr::include_graphics(png_tipopiel)
cat("\n\n")
knitr::asis_output(paste0(mensaje_piel, "\n\n"))
cat("\n\n")


knitr::asis_output(paste0("### ", titulo_riesgo, "\n\n"))
cat("\n\n")
knitr::include_graphics(png_virus)
cat("\n\n")
knitr::asis_output(paste0(mensaje_virus, "\n\n"))
cat("\n\n")
knitr::asis_output(paste0(virus_paraCompleto, "\n\n"))
cat("\n\n")


knitr::asis_output(paste0("### ", titulo_fagos, "\n\n"))
cat("\n\n")
knitr::include_graphics(png_fagos)
cat("\n\n")
knitr::asis_output(paste0(mensaje_fagos, "\n\n"))
cat("\n\n")

knitr::asis_output(paste0("### ", titulo_devolucion, "\n\n"))
cat("\n\n")
knitr::asis_output(paste0(dev_scores, "\n\n"))
cat("\n\n")
knitr::asis_output(paste0(dev_edad, dev_piel, dev_virus, "\n\n"))
cat("\n\n")
cat("\n\n")
cat("\n\n")
knitr::asis_output(paste0(texto_pre_recomendacion, "\n\n"))
cat("\n\n")
knitr::asis_output(paste0("### ", titulo_recomendacion, "\n\n"))
cat("\n\n")
knitr::asis_output(paste0(texto_recomendacion, "\n\n"))
cat("\n\n")
knitr::asis_output(paste0(recomendaciones, "\n\n"))
cat("\n\n")
knitr::asis_output(paste0("### ", titulo_final, "\n\n"))
cat("\n\n")
knitr::asis_output(paste0("#### ",texto_final, "\n\n"))
cat("\n\n")

knitr::asis_output(paste0("### ", titulo_contacto, "\n\n"))
cat("\n\n")
knitr::asis_output(paste0("### ",texto_contacto, "\n\n"))


```
', patient_dir, patient_dir, patient_dir, patient_dir, patient_dir, patient_dir, patient_dir, patient_dir,patient_dir, patient_dir, patient_dir, patient_dir, patient_dir, patient_dir, patient_dir)

  rmd_file <- sprintf("%s/CompleteContenido_ID%s.Rmd", patient_dir, id_mr)
  writeLines(rmd_content, con = rmd_file)
  output_file <- sprintf("%s/CompleteContenido_ID%s.pdf", patient_dir, id_mr)
  rmarkdown::render(rmd_file, output_file = output_file)

  return(output_file)
}
