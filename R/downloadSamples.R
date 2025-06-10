#' @import googledrive
#' @title Download a folder from Drive
#' @description Downloads samples from a folder in googledrive and saves them respecting the structure needed.
#' @param folder_url is a url.
#' @param local_folder is the path where all the folders will be saved.
#' @examples downloadFolderDrive(folder_url = "https://drive.google.com/drive/folders/1TayWJGJK9KJHExisiTNyk_cM-fxDCR4n?usp=drive_link", local_folder = "/media/4tb2/Daniela/Biota/Muestras")
#' @export

downloadFolderDrive <- function(folder_url, local_folder) {

  folder_id <- sub(".*/folders/(.*)\\?usp=drive_link", "\\1", folder_url)
  folder_files <- as.data.frame(drive_ls(as_id(folder_id)))

  # Iterar sobre la lista de archivos y descargarlos

  for (i in 1:nrow(folder_files)) {
    file_id <- folder_files[i,"id"]
    file_name <- folder_files[i,"name"]

    #Crea una carpeta donde se van a guardar las muestras
    id <- strsplit(file_name, split="_")[[1]][1]
    dir.create(paste(local_folder,"/", id, sep=""))

    gz <- ifelse(grepl("gz", file_name), ".gz", "")
    R <- ifelse(grepl("R1", file_name), "R1", "R2")
    local_path <- paste(local_folder,"/", id, "/", id, "_S04_L001_", R, "_001.fastq", gz, sep="")
    #setwd(local_path)
    print(file_name)
    if(!file.exists(local_path)) {
      drive_download(file_id, path = local_path, overwrite = TRUE)
    }
  }
}
