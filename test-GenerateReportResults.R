#library(extrafont)
library(readxl)
library(data.table)
library(openxlsx)
library(tidyr)
library(PipelineBiota)
library(stringr)

font_import()  # Importar todas las fuentes del sistema
fonts()


#CONTROLAR INDICADOR DE RANGO ETARIO DE: 13, 197, 29, 33, 37, 40, 41-90-113-115-138-163-177 (muy jovenes), 118-130-133-140-141-159-162-176-182-207-246-248 (muy viejo)
#cambiar >70 que vaya de 70 a 80
#Todos los que tienne 60 los tira para arriba

generateReportResults(patient_dir = "/home/daniela/Daniela/Biota/Muestras/73m/4")
generateReportResults(patient_dir = "/home/daniela/Daniela/Biota/Muestras/73m/29")
generateReportResults(patient_dir = "/home/daniela/Daniela/Biota/Muestras/73m/140")
generateReportResults(patient_dir = "/home/daniela/Daniela/Biota/Muestras/73m/160")
generateReportResults(patient_dir = "/home/daniela/Daniela/Biota/Muestras/73m/269")

generateReportResults(patient_dir = "/home/daniela/Daniela/Biota/Muestras/73m/1")

#piel grasa
generateSimpleReport(patient_dir = "/home/daniela/Daniela/Biota/Muestras/73m/79")
generateCompleteReport(patient_dir = "/home/daniela/Daniela/Biota/Muestras/73m/79")
generateSimpleReport(patient_dir = "/home/daniela/Daniela/Biota/Muestras/73m/194")
generateSimpleReport(patient_dir = "/home/daniela/Daniela/Biota/Muestras/73m/192")
generateSimpleReport(patient_dir = "/home/daniela/Daniela/Biota/Muestras/73m/37")

#acne:
generateSimpleReport(patient_dir = "/home/daniela/Daniela/Biota/Muestras/acne/217")
generateCompleteReport(patient_dir = "/home/daniela/Daniela/Biota/Muestras/acne/217")
generateSimpleReport(patient_dir = "/home/daniela/Daniela/Biota/Muestras/acne/25")
generateCompleteReport(patient_dir = "/home/daniela/Daniela/Biota/Muestras/acne/25")

#CONOCIDOS:
generateSimpleReport(patient_dir = "/home/daniela/Daniela/Biota/Muestras/73m/273")
generateCompleteReport(patient_dir = "/home/daniela/Daniela/Biota/Muestras/73m/273")

generateSimpleReport(patient_dir = "/home/daniela/Daniela/Biota/Muestras/73m/1")
generateSimpleReport(patient_dir = "/home/daniela/Daniela/Biota/Muestras/73m/152")
generateSimpleReport(patient_dir = "/home/daniela/Daniela/Biota/Muestras/73m/156")
generateSimpleReport(patient_dir = "/home/daniela/Daniela/Biota/Muestras/73m/9")
generateSimpleReport(patient_dir = "/home/daniela/Daniela/Biota/Muestras/73m/79")

generateSimpleReport(patient_dir = "/home/daniela/Daniela/Biota/Muestras/73m/197")
generateSimpleReport(patient_dir = "/home/daniela/Daniela/Biota/Muestras/73m/45")


generateSimpleReport(patient_dir = "/home/daniela/Daniela/Biota/Muestras/73m/37")
generateCompleteReport(patient_dir = "/home/daniela/Daniela/Biota/Muestras/73m/37")


generateSimpleReport(patient_dir = "/home/daniela/Daniela/Biota/Muestras/3")
generateCompleteReport(patient_dir = "/home/daniela/Daniela/Biota/Muestras/3")

#piel seca
generateReportResults(patient_dir = "/home/daniela/Daniela/Biota/Muestras/73m/11")
generateReportResults(patient_dir = "/home/daniela/Daniela/Biota/Muestras/73m/13")
generateReportResults(patient_dir = "/home/daniela/Daniela/Biota/Muestras/73m/207")
generateReportResults(patient_dir = "/home/daniela/Daniela/Biota/Muestras/73m/248")
generateReportResults(patient_dir = "/home/daniela/Daniela/Biota/Muestras/73m/118-1")

#tabaco si:
generateReportResults(patient_dir = "/home/daniela/Daniela/Biota/Muestras/73m/108")
generateReportResults(patient_dir = "/home/daniela/Daniela/Biota/Muestras/73m/133")
generateReportResults(patient_dir = "/home/daniela/Daniela/Biota/Muestras/73m/164")
generateReportResults(patient_dir = "/home/daniela/Daniela/Biota/Muestras/73m/269")
generateReportResults(patient_dir = "/home/daniela/Daniela/Biota/Muestras/73m/27")


list_patients <- list.dirs("/home/daniela/Daniela/Biota/Muestras/73m", recursive = FALSE, full.names= TRUE)
list_patients <- list.dirs("~/Daniela/Biota/Muestras/acne", recursive = FALSE, full.names= TRUE)
list_patients <- list_patients[-c(1:55)]
list_patients <- list_patients[-c(1:15)]
list_patients <- list_patients[-c(1:4)]

list_patients <- list_patients[-c(74, 70, 54)]

#Resolver error con 197, 4, 45!

patient_dir <- list_patients[54]

source = "KRAKEN"
for(patient_dir in list_patients) {
  print(patient_dir)

  #tryCatch({
  #id <- basename(patient_dir)
  #id_mr <- str_pad(id, width = 4, side = "left", pad = "0")
  #report <- sprintf("%s/MR_%s.pdf", patient_dir, id_mr)

   #if(!(file.exists(report))) {
    generateSimpleReport(patient_dir = patient_dir)
    message(sprintf("Reporte Simple generado exitosamente para ID: %s", patient_dir))
  #}
      #}, error = function(e) {
  #  message(sprintf("Error al generar el reporte simple para ID: %s. Error: %s", patient_dir, e$message))
    # Continuar con el siguiente ID
  #})

  #tryCatch({
  #  generateCompleteReport(patient_dir = patient_dir)
  #  message(sprintf("Reporte Simple generado exitosamente para ID: %s", patient_dir))
  #}, error = function(e) {
  #  message(sprintf("Error al generar el reporte simple para ID: %s. Error: %s", patient_dir, e$message))
    # Continuar con el siguiente ID
  #})
}


#Subir los reportes a drive
library(googledrive)
patient_dir <- list_patients[1]
#carpeta de drive:
"https://drive.google.com/drive/folders/1Ublc3XIIRSkd39XKv8cepLR_Xls3kLNj?usp=drive_link"

for(patient_dir in list_patients) {
  print(patient_dir)
  id <- basename(patient_dir)
  id_mr <- str_pad(id, width = 4, side = "left", pad = "0")
  report <- sprintf("%s/MR_%s.pdf", patient_dir, id_mr)

  if(file.exists(report)) {
    folder_drive <- drive_get(as_id("1Ublc3XIIRSkd39XKv8cepLR_Xls3kLNj"))
    drive_upload(media = report, path = folder_drive)
  }

}



