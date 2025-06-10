#' @title Split Lanes
#' @description kvjdfnkvjdf
#' @param patients_dir fskbfjsf
#' @examples splitLanes(patients_dir = "/media/16TBDisk/Daniela/Biota/Muestras2daTanda/5/concat5_04/yoconmuestraheritas")
#' @export

splitLanes <- function(patients_dir) {
  muestras <-  list.dirs(path = patients_dir, full.names = TRUE, recursive = FALSE)
  for (m in muestras) {
    #m <- "/media/4tb1/Daniela/Biota/5"
    file_list <- list.files(path = m, full.names = TRUE, recursive = FALSE)
    gzip <- ifelse(length(nchar(file_list[endsWith(file_list, "R1.fastq.gz")])) == 0, "", ".gz")
    fileR1 <- file_list[endsWith(file_list, sprintf("R1.fastq%s", gzip))]
    fileR2 <- file_list[endsWith(file_list, sprintf("R2.fastq%s", gzip))]

    #f <- fileR1
    for (f in c(fileR1, fileR2)) {
      print(f)
      R <- ifelse(grepl("R1", basename(f)), "R1", "R2" )
      id <- basename(dirname(f))

      if (!(file.exists(paste(dirname(f), "/", id, "_S", id, "_L001_", R, "_001.fastq.gz", sep = "")))) {
        fastq_data <- readFastq(f)

        #Dividir en los 4 lanes:
        indices_1 <- grep(":1:", id(fastq_data))
        indices_2 <- grep(":2:", id(fastq_data))
        indices_3 <- grep(":3:", id(fastq_data))
        indices_4 <- grep(":4:", id(fastq_data))

        # Crear dos nuevos objetos FASTQ separados basados en los ?ndices obtenidos
        lane1 <- fastq_data[indices_1]
        lane2 <- fastq_data[indices_2]
        lane3 <- fastq_data[indices_3]
        lane4 <- fastq_data[indices_4]

        head(id(lane1))
        head(id(lane2))
        head(id(lane3))
        head(id(lane4))
        #total_reads <- length(fastq_data)
        #reads_per_lane <- ceiling(total_reads / 4)
        #head(id(fastq_data[total_reads]))

        writeFastq(lane1, paste(dirname(f), "/", id,"_S", id ,"_L001_", R, "_001.fastq.gz", sep = ""))
        writeFastq(lane2, paste(dirname(f), "/", id,"_S", id, "_L002_", R, "_001.fastq.gz", sep = ""))
        writeFastq(lane3, paste(dirname(f), "/", id, "_S", id,"_L003_", R, "_001.fastq.gz", sep = ""))
        writeFastq(lane4, paste(dirname(f), "/", id, "_S", id, "_L004_", R, "_001.fastq.gz", sep = ""))
      }
    }

  }
}
