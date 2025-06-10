#' @title Format names
#' @description It changes the format in which samples are named in order to match with dragen requirments
#' @param patient_dir path to the folder that stores all the samples
#' @export
formatNames <- function(patients_dir) {

  #m <- muestras[[1]]
  muestras <-  list.dirs(path = patients_dir, full.names = TRUE, recursive = FALSE)
  for (m in muestras) {
    print(m)
    file_list <- list.files(path = m, full.names = TRUE, recursive = FALSE)
    gzip <- ifelse(length(nchar(file_list[endsWith(file_list, "R1.fastq.gz")])) == 0, "", ".gz")
    fileR1 <- file_list[endsWith(file_list, sprintf("R1.fastq%s", gzip))]
    fileR2 <- file_list[endsWith(file_list, sprintf("R2.fastq%s", gzip))]

    for (f in c(fileR1, fileR2)) {
      print(f)
      R <- ifelse(grepl("R1", basename(f)), "R1", "R2" )
      id <- basename(dirname(f))
      if (!(file.exists(paste(dirname(f), "/", id, "_S", id, "_L001_", R, "_001.fastq.gz", sep = "")))) {
        print("cambiado")
        file.rename(paste(dirname(f), "/", id,"_L001_", R, "_001.fastq.gz", sep = ""), paste(dirname(f), "/", id,"_S", id ,"_L001_", R, "_001.fastq.gz", sep = ""))
        file.rename(paste(dirname(f), "/", id,"_L002_", R, "_001.fastq.gz", sep = ""), paste(dirname(f), "/", id,"_S", id ,"_L002_", R, "_001.fastq.gz", sep = ""))
        file.rename(paste(dirname(f), "/", id,"_L003_", R, "_001.fastq.gz", sep = ""), paste(dirname(f), "/", id,"_S", id ,"_L003_", R, "_001.fastq.gz", sep = ""))
        file.rename(paste(dirname(f), "/", id,"_L004_", R, "_001.fastq.gz", sep = ""), paste(dirname(f), "/", id,"_S", id ,"_L004_", R, "_001.fastq.gz", sep = ""))

      }
    }
  }
}

