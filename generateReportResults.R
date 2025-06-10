#' @title Generate Results' Report
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
generateReportResults <- function(patient_dir) {

  id <- basename(patient_dir)
  patients_dir <- dirname(patient_dir)
  MetadataB <- as.data.frame(read_excel("~/Daniela/Biota/PipelineBiota/data/Metadata-soloColumnasUsables.xlsx"))
  BIOTALIFE_SKIN_Respuestas_ <- read_excel("~/Daniela/Biota/BIOTALIFE SKIN  (Respuestas).xlsx", sheet = 4)

  #Obtener info gral
  BIOTALIFE_SKIN_Respuestas_$ID  <- gsub("\\.0$", "", as.character(BIOTALIFE_SKIN_Respuestas_$ID))
  BIOTALIFE_SKIN_Respuestas_$ID[which(BIOTALIFE_SKIN_Respuestas_$ID == "118")] <- "118-1"
  BIOTALIFE_SKIN_Respuestas_$ID[which(BIOTALIFE_SKIN_Respuestas_$ID == "184")] <- "184-1"

  info <- BIOTALIFE_SKIN_Respuestas_[which(BIOTALIFE_SKIN_Respuestas_$ID == id),]
  colnames(info)
  info <- info[, which(colnames(info) %in% c("ID", "1. Nombre y Apellido","3. Sexo", "2. Fecha de Nacimiento", "Edad", "Fecha Cita" ))]
  info <- as.data.frame(t(info))
  info <- data.frame("Caracteristica" = rownames(info),"Informacion" = info$V1)
  info$Caracteristica <- c("ID", "Nombre y Apellido", "Fecha de nacimiento", "Sexo", "Fecha Cita", "Edad")

  # Obtener datos de las tablas
  counts_tax <- counts_Tax(patients_dir = patient_dir, source = "KRAKEN")
  reportCountsTax(id=id, MetadataB = MetadataB, patients_dir)

  tabla_subesp_mas_abundantes <- reportMasAbundantes(id = id, MetadataB = MetadataB)

  indicadorRangoEtario(id = id, MetadataB = MetadataB)

  titulo1 <- "1. ¿Qué es el microbioma cutáneo y por qué es clave para la salud de tu piel?"
  texto1 <- "El microbioma cutáneo es una comunidad compleja de microorganismos, como bacterias, hongos y virus, que habitan la superficie de la piel. Este ecosistema desempeña un papel esencial en la salud de la piel al actuar como una barrera protectora contra microorganismos dañinos, ayudar a mantener el equilibrio del pH y fortalecer la función inmunológica. La composición del microbioma puede variar según factores como la edad, el sexo biológico, la dieta y el estilo de vida. Un microbioma equilibrado es fundamental para mantener la piel sana y prevenir problemas cutáneos como el acné, la dermatitis atópica y la rosácea. Estudios recientes muestran que la diversidad y estabilidad del microbioma cutáneo están relacionadas con la integridad de la barrera cutánea y la respuesta inmunitaria local. Un desequilibrio en este ecosistema puede contribuir a la aparición de enfermedades de la piel y reducir su capacidad para combatir patógenos."

  titulo2 <- " 2. ¿Cómo interpretar tus resultados?"
  texto2 <- "Este informe analiza la abundancia relativa de los microorganismos que componen tu microbioma cutáneo. Los microorganismos se clasifican en diferentes niveles taxonómicos: filo, clase, orden, familia, género, especie y subespecie. Aquí, cada microorganismo se evalúa para determinar si su presencia es baja, equilibrada o alta en comparación con nuestra base de datos de referencia de piel saludable en Biotalife Skin. Es importante recordar que un valor fuera del rango deseado no implica necesariamente una enfermedad. La variabilidad individual en el microbioma es alta, y algunos desequilibrios pueden ser temporales o estar relacionados con factores externos como cambios en el entorno, uso de cosméticos o medicamentos. "

  titulo3 <- "3. ¿En qué consiste nuestro estudio?"
  texto3 <- "Este estudio utiliza tecnología de secuenciación de nueva generación (NGS) de Illumina para analizar el ADN de los microorganismos presentes en una muestra de la piel obtenida mediante un hisopado estéril. La secuenciación permite identificar con precisión los diferentes microorganismos, brindando una visión detallada de la composición del microbioma cutáneo. Esta tecnología nos permite comparar tu microbioma con una base de datos de referencia y evaluar su equilibrio."

  titulo4 <- "4. Descargo de responsabilidad"
  texto4 <- "Los resultados se basan en la comparación con la base de datos de referencia de Biotalife Skin, que incluye muestras de personas con piel saludable. Cada parámetro se clasifica como Bajo, Equilibrado o Alto utilizando un algoritmo desarrollado por nuestro equipo. Además, hemos identificado ciertas especies que actúan como biomarcadores de la edad de la piel, lo que nos permite predecir si la presencia y cantidad de estos biomarcadores corresponden a tu rango de edad biológica. Es importante entender que este informe es una herramienta complementaria y no un diagnóstico médico. Un valor fuera del rango deseado podría estar relacionado con factores temporales o individuales, y no necesariamente con una condición patológica."

  titulo5 <- "5. Descripción de la composición de tu microbioma cutáneo."
  texto5 <- "Tu microbioma cutáneo es único, como una huella digital bacteriana. En nuestra base de datos de Biotalife Skin, hemos identificado una amplia variedad de filos, familias, géneros y especies presentes en pieles saludables. En esta sección, detallaremos la diversidad y composición de tu microbioma, destacando los grupos más representativos."

  subtitulo_filos <- "5.1. Filos más abundantes"
  texto_filos <- "Los filos bacterianos más comunes en la piel son Actinobacteria y Firmicutes. Actinobacteria incluye géneros esenciales para la salud de la piel, como Cutibacterium y Corynebacterium. Firmicutes incluye Staphylococcus y Streptococcus, que pueden tener efectos beneficiosos o perjudiciales dependiendo de su equilibrio. Un desequilibrio en estos filos puede estar asociado con afecciones cutáneas como el acné o la dermatitis."

  subtitulo_generos <- "5.2. Géneros más abundantes"
  Genus_masAbundantes <- read_excel(sprintf("%s/Tabla_Genus_masAbundantes.xlsx", patient_dir))
  generos_top3 <- paste0(Genus_masAbundantes$Géneros[1:3], collapse = ", ")
  if(any(grepl("virus", Genus_masAbundantes$Géneros))) {
    warning_virus_generos <- "La presencia de virus dentro de los géneros más abundantes de la piel es un factor de riesgo para la salud cutánea. Se recomienda solicitar una consulta con un dermatólogo a la brevedad."
  } else {
    warning_virus_generos <- ""
  }
  texto_generos <- sprintf("La piel está colonizada por una variedad de géneros microbianos. En tu caso, %s son los más representativos. Estos géneros desempeñan roles importantes en la salud de la piel, desde la producción de ácidos grasos hasta la protección contra patógenos. Por ejemplo, Cutibacterium acnes produce ácidos grasos que ayudan a mantener el pH ácido de la piel, protegiéndola de microorganismos dañinos. %s", generos_top3, warning_virus_generos)

  subtitulo_especies <- "5.3. Especies más abundantes"
  esp_masAbundantes <- read_excel(sprintf("%s/Tabla_Species_masAbundantes.xlsx", patient_dir))
  if(any(grepl("virus", esp_masAbundantes$Especies))) {
    warning_virus_especies <- "La presencia de virus dentro de las especies más abundantes de la piel es un factor de riesgo para la salud cutánea. Se recomienda solicitar una consulta con un dermatólogo a la brevedad."
  } else {
    warning_virus_especies <- ""
  }
  texto_especies <- sprintf("Cada persona tiene un sello personal de especies bacterianas. Las especies más abundantes en tu microbioma brindan información crucial sobre la salud y diversidad de tu piel. Esta sección analiza las especies más destacadas y su posible impacto en la salud cutánea. Algunas especies pueden ser beneficiosas, mientras que otras, si están en desequilibrio, podrían estar asociadas con afecciones cutáneas. %s", warning_virus_especies)

  subtitulo_subespecies <- "5.4. Subespecies más abundantes"
  subesp_masAbundantes <- read_excel(sprintf("%s/Tabla_SubSpecies_masAbundantes.xlsx", patient_dir))
  if(any(grepl("virus", subesp_masAbundantes$SubEspecies))) {
    warning_virus <- "La presencia de virus dentro de las subespecies más abundantes de la piel es un factor de riesgo para la salud cutánea. Se recomienda solicitar una consulta con un dermatólogo a la brevedad."
  } else {
    warning_virus <- ""
  }
  texto_subespecies <- sprintf("En algunos casos, se pueden identificar subespecies, lo que proporciona una visión aún más detallada de tu microbioma. Las subespecies pueden variar entre personas y ofrecer información específica sobre la composición y variaciones en tu microbioma cutáneo. %s", warning_virus)

  subtitulo_biomarcador_rango <- "6. Biomarcadores - Envejecimiento de la piel."
  titulotipopiel <- "7. Biomarcadores - Tipo de piel."
  titulo_riesgo <- "8. Biomarcadores - Factores de riesgo."
  titulo_habitos <- "9. Biomarcadores - Hábitos y factores externos."
  titulo_fagos <- "10. Biomarcadores - Factores beneficiosos."

  titulo_final <- "11. Comentarios finales"
  texto_final <- "¡Gracias por participar en el estudio de Biotalife Skin! Tu colaboración es fundamental para avanzar en la comprensión del microbioma cutáneo. Estamos comprometidos con tu bienestar y con desarrollar soluciones personalizadas para el cuidado de la piel. Recuerda que cada pequeño cambio puede marcar una gran diferencia en la salud de tu piel."


  # Verificar si el informe ya existe
  #if (file.exists(sprintf("%s/ReportResults_Biota_ID%s.pdf", patient_dir, id))) {
  #  return(message("¡El informe de resultados de este paciente ya ha sido generado!"))
  #}

  # Contenido del archivo Rmd
  rmd_content <- sprintf('
---
title: "INFORME COMPLETO DE LA COMPOSICIÓN DEL MICROBIOMA CUTÁNEO - Muestra %s"
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

CantidadesTax <- read_excel(sprintf("%s/Tabla_CantidadesTax.xlsx", patient_dir))

SubEspecies_masAbundantes <- read_excel(sprintf("%s/Tabla_SubSpecies_masAbundantes.xlsx", patient_dir))
Especies_masAbundantes <- read_excel(sprintf("%s/Tabla_Species_masAbundantes.xlsx", patient_dir))
Genus_masAbundantes <- read_excel(sprintf("%s/Tabla_Genus_masAbundantes.xlsx", patient_dir))
Phylum_masAbundantes <- read_excel(sprintf("%s/Tabla_Phylum_masAbundantes.xlsx", patient_dir))

SubEspecies_RangoEtario <- read_excel(sprintf("%s/Tabla_SubEspecies_RangoEtario.xlsx", patient_dir))

c <- reportCountsTax(id=id, MetadataB = MetadataB)
cants_plot <- c[[2]]

out <- reportMasAbundantes(id = id, MetadataB = MetadataB)
subspecies_plot <- out[[1]]
species_plot <- out[[2]]
genus_plot <- out[[3]]
phylum_plot <- out[[4]]

png_barra <- sprintf("%s/rango_etario_barra.png", patient_dir)
png_biomarcadores <- (sprintf("%s/PieChart_BiomarcadoresRangoEtario.png", patient_dir))
out <- indicadorRangoEtario(id = id, MetadataB = MetadataB)
rango_predicho <- out[[2]]
plot_rango_real <- out[[3]]
combined_plot <- out[[4]]
mensaje <- out[[5]]
barra_plot <- out[[6]]

out <- indicadorTipoPiel(id = id, MetadataB = MetadataB)
mensaje_piel <- out[[2]]
png_tipopiel <- sprintf("%s/tipodepiel_barra.png", patient_dir)

out <- indicadorHabitos(id = id, MetadataB = MetadataB)
mensaje_alcohol <- out[[1]]
mensaje_tabaco <- out[[2]]
mensaje_maquillaje <- out[[3]]
mensaje_actividadfisica <- out[[4]]
mensajes_habitos <- paste0(mensaje_alcohol, mensaje_tabaco, mensaje_maquillaje, mensaje_actividadfisica, collapse=".")
png_habitos <- sprintf("%s/habitos_barra.png", patient_dir)

out <- indicadorVirus(id = id, MetadataB = MetadataB)
mensaje_virus <- out[[1]]
png_virus <- sprintf("%s/virus_barra.png", patient_dir)

out <- indicadorFagos(id = id, MetadataB = MetadataB)
mensaje_fagos <- out[[1]]
png_fagos <- sprintf("%s/fagos_barra.png", patient_dir)


titulo1 <- titulo1
texto1 <- texto1
titulo2 <- titulo2
texto2 <- texto2
titulo3 <- titulo3
texto3 <- texto3
titulo4 <- titulo4
texto4 <- texto4
titulo5 <- titulo5
texto5 <- texto5


subtitulo_filos <- subtitulo_filos
texto_filos <- texto_filos

subtitulo_generos <- subtitulo_generos
texto_generos <- texto_generos

subtitulo_especies <- subtitulo_especies
texto_especies <- texto_especies

subtitulo_subespecies <- subtitulo_subespecies
texto_subespecies <- texto_subespecies

titulotipopiel <- titulotipopiel
titulo_riesgo <- titulo_riesgo
titulo_habitos <- titulo_habitos
titulo_fagos <- titulo_fagos

titulo_final <- titulo_final
texto_final <- texto_final
titulo_contacto <- "12. Contactános"
texto_contacto <- "Si tenés alguna consulta relacionada con tu informe, no dudes en escribirnos a: contacto@biotalifeskin.com. Además, podés ingresar a nuestra página web para mayor información: www.biotalifeskin.com. Encontranos en instagram como @biotalife."

```

```{r, echo=FALSE}

knitr::asis_output(paste0("# ", "Índice", "\n\n"))
cat("\n\n")
knitr::asis_output(paste0(titulo1, "\n\n"))
knitr::asis_output(paste0(titulo2, "\n\n"))
knitr::asis_output(paste0(titulo3, "\n\n"))
knitr::asis_output(paste0(titulo4, "\n\n"))
knitr::asis_output(paste0(titulo5, "\n\n"))
knitr::asis_output(paste0(subtitulo_biomarcador_rango, "\n\n"))
knitr::asis_output(paste0(titulotipopiel, "\n\n"))
knitr::asis_output(paste0(titulo_habitos, "\n\n"))
knitr::asis_output(paste0(titulo_riesgo, "\n\n"))
knitr::asis_output(paste0(titulo_riesgo, "\n\n"))
knitr::asis_output(paste0(titulo_fagos, "\n\n"))
knitr::asis_output(paste0(titulo_final, "\n\n"))
knitr::asis_output(paste0(titulo_contacto, "\n\n"))

cat("\n\n")
knitr::asis_output(paste0(" ", "\n\n"))
cat("\n\n")
knitr::asis_output(paste0(" ", "\n\n"))
cat("\n\n")
knitr::asis_output(paste0(" ", "\n\n"))
cat("\n\n")
knitr::asis_output(paste0(" ", "\n\n"))
cat("\n\n")
knitr::asis_output(paste0(" ", "\n\n"))
cat("\n\n")
knitr::asis_output(paste0(" ", "\n\n"))
cat("\n\n")
knitr::asis_output(paste0(" ", "\n\n"))
cat("\n\n")
knitr::asis_output(paste0(" ", "\n\n"))
cat("\n\n")
knitr::asis_output(paste0(" ", "\n\n"))
cat("\n\n")
cat("\n\n")
knitr::asis_output(paste0(" ", "\n\n"))
cat("\n\n")
knitr::asis_output(paste0(" ", "\n\n"))
cat("\n\n")
knitr::asis_output(paste0(" ", "\n\n"))
cat("\n\n")
cat("\n\n")
knitr::asis_output(paste0(" ", "\n\n"))
cat("\n\n")
knitr::asis_output(paste0(" ", "\n\n"))
cat("\n\n")
knitr::asis_output(paste0(" ", "\n\n"))
cat("\n\n")
cat("\n\n")
knitr::asis_output(paste0(" ", "\n\n"))
cat("\n\n")
knitr::asis_output(paste0(" ", "\n\n"))
cat("\n\n")
knitr::asis_output(paste0(" ", "\n\n"))
cat("\n\n")

knitr::kable(InfoGral, caption = "Información general de la muestra")

knitr::asis_output(paste0("## ", titulo1, "\n\n"))
knitr::asis_output(paste0(texto1, "\n\n"))

knitr::asis_output(paste0("## ", titulo2, "\n\n"))
knitr::asis_output(paste0(texto2, "\n\n"))

knitr::asis_output(paste0("## ", titulo3, "\n\n"))
knitr::asis_output(paste0(texto3, "\n\n"))

knitr::asis_output(paste0("## ", titulo4, "\n\n"))
knitr::asis_output(paste0(texto4, "\n\n"))

knitr::asis_output(paste0("## ", titulo5, "\n\n"))
knitr::asis_output(paste0(texto5, "\n\n"))

cat("\n\n")

knitr::asis_output(paste0("## ", subtitulo_filos, "\n\n"))
knitr::asis_output(paste0(texto_filos, "\n\n"))
knitr::kable(Phylum_masAbundantes, caption = "Filos más abundantes")
plot(phylum_plot)
cat("\n\n")


knitr::asis_output(paste0("## ", subtitulo_generos, "\n\n"))
knitr::asis_output(paste0(texto_generos, "\n\n"))
knitr::kable(Genus_masAbundantes, caption = "Géneros más abundantes")
plot(genus_plot)
cat("\n\n")

knitr::asis_output(paste0("## ", subtitulo_especies, "\n\n"))
knitr::asis_output(paste0(texto_especies, "\n\n"))
knitr::kable(Especies_masAbundantes, caption = "Especies más abundantes")
plot(species_plot)
cat("\n\n")


knitr::asis_output(paste0("## ", subtitulo_subespecies, "\n\n"))
knitr::asis_output(paste0(texto_subespecies, "\n\n"))
knitr::kable(SubEspecies_masAbundantes, caption = "Subespecies más abundantes")
plot(subspecies_plot)
cat("\n\n")


knitr::asis_output(paste0("## ", subtitulo_biomarcador_rango, "\n\n"))
#knitr::kable(SubEspecies_RangoEtario, caption = "Biomarcadores Rango Etario")
#plot(combined_plot)
cat("\n\n")
knitr::include_graphics(png_barra)
cat("\n\n")
knitr::asis_output(paste0(mensaje, "\n\n"))
cat("\n\n")

knitr::asis_output(paste0("## ", titulotipopiel, "\n\n"))
cat("\n\n")
knitr::include_graphics(png_tipopiel)
cat("\n\n")
knitr::asis_output(paste0(mensaje_piel, "\n\n"))
cat("\n\n")

knitr::asis_output(paste0("## ", titulo_habitos, "\n\n"))
cat("\n\n")
knitr::include_graphics(png_habitos)
cat("\n\n")
knitr::asis_output(paste0(mensajes_habitos, "\n\n"))
cat("\n\n")

knitr::asis_output(paste0("## ", titulo_riesgo, "\n\n"))
cat("\n\n")
knitr::include_graphics(png_virus)
cat("\n\n")
knitr::asis_output(paste0(mensaje_virus, "\n\n"))
cat("\n\n")

knitr::asis_output(paste0("## ", titulo_fagos, "\n\n"))
cat("\n\n")
knitr::include_graphics(png_fagos)
cat("\n\n")
knitr::asis_output(paste0(mensaje_fagos, "\n\n"))
cat("\n\n")

knitr::asis_output(paste0("## ", titulo_final, "\n\n"))
cat("\n\n")
knitr::asis_output(paste0("### ",texto_final, "\n\n"))
cat("\n\n")

knitr::asis_output(paste0("## ", titulo_contacto, "\n\n"))
cat("\n\n")
knitr::asis_output(paste0("### ",texto_contacto, "\n\n"))


```
', id, patient_dir, patient_dir, patient_dir, patient_dir, patient_dir, patient_dir,patient_dir, patient_dir, patient_dir, patient_dir, patient_dir, patient_dir, patient_dir)

  rmd_file <- sprintf("%s/ReportResults_Biota_ID%s.Rmd", patient_dir, id)
  writeLines(rmd_content, con = rmd_file)
  output_file <-sprintf("%s/ReportResults_Biota_ID%s.pdf", patient_dir, id)
  rmarkdown::render(rmd_file, output_file = output_file)
  #output_file <- sprintf("%s/ReportResults_Biota.doc", patient_dir)
  #rmarkdown::render(rmd_file, output_file = output_file)

}






