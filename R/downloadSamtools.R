#' @title downloadSamtools
#' @description Downloads and decompresses the Samtools software
#' @return The path where the .exe file is located
downloadSamtools <- function() {
  soft_directory <- sprintf("%s/PipelineBiota-Softwares", dirname(system.file(package = "PipelineBiota")))

  tryCatch(
    {
      system2(sprintf("%s/Samtools/samtools-1.16.1/samtools", soft_directory))
    },
    error = function(e) {
      message("The installation of Samtools will begin now. Check if all required packages are downloaded.")
      print(e)

      dir.create(sprintf("%s/Samtools", soft_directory))
      dir <- sprintf("%s/Samtools", soft_directory)
      URL <- "https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2"
      system2("wget", args = c(URL, "-P", dir), wait = TRUE, stdout = NULL, stderr = NULL)

      system2("bzip2", sprintf("-d %s/Samtools/samtools-1.16.1.tar.bz2", soft_directory))
      system2("tar", c("-xvf", sprintf("%s/Samtools/%s -C %s", soft_directory, list.files(sprintf("%s/Samtools", soft_directory)), dir)))

      file.remove(sprintf("%s/samtools-1.16.1.tar", dir))

      samtools_dir <- sprintf("%s/samtools-1.16.1", dir)
      system(paste("cd", shQuote(samtools_dir), "&& ./configure"))
      system(paste("cd", shQuote(samtools_dir), "&& make"))

      #Escribo el path en el txt
      softwares <- readLines(sprintf("%s/PipelineBiota-Softwares/path_to_soft.txt", dirname(system.file(package = "PipelineBiota"))))
      Samtools <- sprintf("%s/Samtools/samtools-1.16.1/samtools", soft_directory)
      softwares_actualizado <- c(softwares, sprintf("Samtools %s", Samtools))
      write(softwares_actualizado, file = sprintf("%s/PipelineBiota-Softwares/path_to_soft.txt", dirname(system.file(package = "PipelineBiota"))))

    },
    warning = function(w) {
      message("The installation of Samtools will begin now. Check if all required packages are downloaded.")

      dir.create(sprintf("%s/Samtools", soft_directory))
      dir <- sprintf("%s/Samtools", soft_directory)
      URL <- "https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2"
      system2("wget", args = c(URL, "-P", dir), wait = TRUE, stdout = NULL, stderr = NULL)

      system2("bzip2", sprintf("-d %s/Samtools/samtools-1.16.1.tar.bz2", soft_directory))
      system2("tar", c("-xvf", sprintf("%s/Samtools/%s -C %s", soft_directory, list.files(sprintf("%s/Samtools", omicsdo_sof)), dir)))

      file.remove(sprintf("%s/samtools-1.16.1.tar", dir))

      samtools_dir <- sprintf("%s/samtools-1.16.1", dir)
      system(paste("cd", shQuote(samtools_dir), "&& ./configure"))
      system(paste("cd", shQuote(samtools_dir), "&& make"))

      #Escribo el path en el txt
      softwares <- readLines(sprintf("%s/PipelineBiota-Softwares/path_to_soft.txt", dirname(system.file(package = "PipelineBiota"))))
      Samtools <- sprintf("%s/Samtools/samtools-1.16.1/samtools", soft_directory)
      softwares_actualizado <- c(softwares, sprintf("Samtools %s", Samtools))
      write(softwares_actualizado, file = sprintf("%s/PipelineBiota-Softwares/path_to_soft.txt", dirname(system.file(package = "PipelineBiota"))))

    },

    finally = {
      message("-.Message from Samtools")
      Samtools <<- sprintf("%s/Samtools/samtools-1.16.1/samtools", soft_directory)
    }
  )

  return(sprintf("%s/Samtools/samtools-1.16.1/samtools", soft_directory))
}
