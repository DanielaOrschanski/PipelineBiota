#' @import googledrive
#' @title Uploades fastq R1 and de-hosted fastq R1 from server to Drive
#' @description Uploads samples from a folder in googledrive and saves them respecting the structure needed.
#' @param folder_url is a url.
#' @param local_folder is the path where all the folders will be taken from.
#' @examples drive_folder_url <- "https://drive.google.com/drive/folders/1vcYZAqK74kxEZcZEe74bOQ_f9Hfu9nVW?usp=drive_link"
#' uploadFolderDrive(local_folder =  "~/Daniela/Biota/Muestras/73m", drive_folder_url = drive_folder_url)
#' @export

library(googledrive)
uploadFolderDrive <- function(local_folder, drive_folder_url) {
  # Obtén el ID de la carpeta de destino en Drive
  folder_id <- sub(".*/folders/(.*)\\?usp=drive_link", "\\1", drive_folder_url)

  # Obtén una lista de los archivos en el folder local
  folders <- list.dirs(local_folder, full.names = TRUE, recursive = FALSE)

  # Iterar sobre los archivos y subirlos
  for (folder in folders) {
    #folder <- folders[[1]]
    print(folder)
    id <- basename(folder)
    #fileR1 <- sprintf("%s/%s_S04_L001_R1_001.fastq.gz", folder, id)
    #file_dhR1 <- sprintf("%s/trimmed/%sDHBo_S04_L001_R1_001.fastq.gz", folder, id)

    txtR1 <- sprintf("%s/%s_S04_L001_R1_001_fastqc/fastqc_data.txt", folder, id)
    txt_dhR1 <- sprintf("%s/trimmed/%sDHBo_S04_L001_R1_001_fastqc/fastqc_data.txt", folder, id)

    if(!file.exists(txt_dhR1)) {
      zip_dhR1 <- sprintf("%s/trimmed/%sDHBo_S04_L001_R1_001_fastqc.zip", folder, id)
      #file.exists(zip_dhR1)
      unzip(zip_dhR1, exdir = sprintf("%s/trimmed", folder))

    }

    existing_files <- drive_ls(as_id(folder_id))$name

    #if(file.exists(fileR1) & !(basename(fileR1) %in% existing_files)) {
    #  drive_upload(fileR1, path = as_id(folder_id), overwrite = FALSE)
    #}
    #if(file.exists(file_dhR1) & !(basename(file_dhR1) %in% existing_files)) {
    #  drive_upload(file_dhR1, path = as_id(folder_id), overwrite = FALSE)
    #}

    nombre_en_drive <- sprintf("%s_fastqc_data.txt", id)
    if(file.exists(txtR1) & !(nombre_en_drive %in% existing_files)) {
      #drive_upload(txtR1, path = as_id(folder_id), overwrite = FALSE)
      drive_upload(txtR1, path = as_id(folder_id), name = nombre_en_drive, overwrite = FALSE)

    }

    nombre_en_drive <- sprintf("%sDHBo_fastqc_data.txt", id)
    if(file.exists(txt_dhR1) & !(nombre_en_drive %in% existing_files)) {
      drive_upload(txt_dhR1, path = as_id(folder_id), name = nombre_en_drive, overwrite = FALSE)
    }

  }
}
