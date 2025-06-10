# BWA
#' @title downloadBWA
#' @description Downloads and decompresses the BWA software
#' @return The path where the .exe file is
#' @export
downloadBWA <- function() {
  #omicsdo_sof <- sprintf("%s/OMICsdoSof", dirname(system.file(package = "OMICsdo")))
  soft_directory <- sprintf("%s/PipelineBiota-Softwares", dirname(system.file(package = "PipelineBiota")))

  omicsdo_sof <- soft_directory
  print(omicsdo_sof)

  tryCatch(
    expr = {
      system(sprintf('%s/BWA/usr/bin/bwa', soft_directory))
    },
    error = function(e) {
      message("Installation of BWA will now begin. Check if all the required packages are downloaded.")
      print(e)

      dir.create(sprintf("%s/BWA", soft_directory))

      bwa_url1 <- "https://download.opensuse.org/repositories/home:/vojtaeus/15.4/x86_64/bwa-0.7.17-lp154.6.1.x86_64.rpm"
      bwa_dir1 <- file.path(soft_directory, "BWA")
      system2("wget", args = c(bwa_url1, "-P", bwa_dir1), wait = TRUE, stdout = NULL, stderr = NULL)

      system2("rpm2cpio", sprintf("%s/bwa-0.7.17-lp154.6.1.x86_64.rpm | cpio -D %s -idmv", bwa_dir1, bwa_dir1), wait = TRUE)

      dir.create(sprintf("%s/BWA/bwa-0.7.17-lp154.6.1.src", omicsdo_sof))
      bwa_url2 <-"https://download.opensuse.org/repositories/home:/vojtaeus/15.4/src/bwa-0.7.17-lp154.6.1.src.rpm"
      bwa_dir2 <- file.path(omicsdo_sof, "BWA/bwa-0.7.17-lp154.6.1.src")
      system2("wget", args = c(bwa_url2, "-P", bwa_dir2), wait = TRUE, stdout = NULL, stderr = NULL)
      system2("rpm2cpio", sprintf("%s/bwa-0.7.17-lp154.6.1.src.rpm | cpio -D %s -idmv", bwa_dir2, bwa_dir2), wait = TRUE)

      bwa_url3 <-"https://download.opensuse.org/repositories/home:/vojtaeus/15.4/i586/bwa-0.7.17-lp154.6.1.i586.rpm"
      system2("wget", args = c(bwa_url3, "-P", bwa_dir1), wait = TRUE, stdout = NULL, stderr = NULL)
      system2("rpm2cpio", sprintf("%s/bwa-0.7.17-lp154.6.1.i586.rpm | cpio -D %s -idmv", bwa_dir1, bwa_dir1), wait = TRUE)

      #Escribo el path en el txt
      softwares <- readLines(sprintf("%s/OMICsdoSof/path_to_soft.txt", dirname(system.file(package = "OMICsdo"))))
      BWA <- sprintf('%s/BWA/usr/bin/bwa', omicsdo_sof)
      softwares_actualizado <- c(softwares, sprintf("BWA %s", BWA))
      write(softwares_actualizado, file = sprintf("%s/OMICsdoSof/path_to_soft.txt", dirname(system.file(package = "OMICsdo"))))

    },
    warning = function(w) {
      message("Installation of BWA will now begin. Check if all the required packages are downloaded.")

      dir.create(sprintf("%s/BWA", omicsdo_sof))

      bwa_url1 <- "https://download.opensuse.org/repositories/home:/vojtaeus/15.4/x86_64/bwa-0.7.17-lp154.6.1.x86_64.rpm"
      bwa_dir1 <- file.path(omicsdo_sof, "BWA")
      system2("wget", args = c(bwa_url1, "-P", bwa_dir1), wait = TRUE, stdout = NULL, stderr = NULL)
      system2("rpm2cpio", sprintf("%s/bwa-0.7.17-lp154.6.1.x86_64.rpm | cpio -D %s -idmv", bwa_dir1, bwa_dir1), wait = TRUE)

      dir.create(sprintf("%s/BWA/bwa-0.7.17-lp154.6.1.src", omicsdo_sof))
      bwa_url2 <-"https://download.opensuse.org/repositories/home:/vojtaeus/15.4/src/bwa-0.7.17-lp154.6.1.src.rpm"
      bwa_dir2 <- file.path(omicsdo_sof, "BWA/bwa-0.7.17-lp154.6.1.src")
      system2("wget", args = c(bwa_url2, "-P", bwa_dir2), wait = TRUE, stdout = NULL, stderr = NULL)
      system2("rpm2cpio", sprintf("%s/bwa-0.7.17-lp154.6.1.src.rpm | cpio -D %s -idmv", bwa_dir2, bwa_dir2), wait = TRUE)

      bwa_url3 <-"https://download.opensuse.org/repositories/home:/vojtaeus/15.4/i586/bwa-0.7.17-lp154.6.1.i586.rpm"
      system2("wget", args = c(bwa_url3, "-P", bwa_dir1), wait = TRUE, stdout = NULL, stderr = NULL)
      system2("rpm2cpio", sprintf("%s/bwa-0.7.17-lp154.6.1.i586.rpm | cpio -D %s -idmv", bwa_dir1, bwa_dir1), wait = TRUE)

      #Escribo el path en el txt
      softwares <- readLines(sprintf("%s/PipelineBiota-Softwares/path_to_soft.txt", dirname(system.file(package = "PipelineBiota"))))
      BWA <- sprintf('%s/BWA/usr/bin/bwa', omicsdo_sof)
      softwares_actualizado <- c(softwares, sprintf("BWA %s", BWA))
      write(softwares_actualizado, file = sprintf("%s/PipelineBiota-Softwares/path_to_soft.txt", dirname(system.file(package = "PipelineBiota"))))

    },

    finally = {
      #BWA <<- sprintf('%s/BWA/usr/bin/bwa', omicsdo_sof)
      message("-.Message from BWA")
    }
  )
  
  #Indexar con bwa la referencia:
  indexBWA <- sprintf("%s", soft_directory)
  if(dir.exists(sprintf("%s/HG38/indexBWA"))){
    message("The reference will be indexed by BWA")
    system(sprintf("bwa index %s", FastaHG38))
  }
  

  return(sprintf('%s/BWA/usr/bin/bwa', omicsdo_sof))
}
