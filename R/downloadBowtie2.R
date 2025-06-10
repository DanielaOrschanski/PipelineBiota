#' @title downloadBowtie2
#' @description Downloads and decompresses the Bowtie2 software
#' @return The path where the ejecution file is
#' @export
downloadBowtie2 <- function(){
  soft_directory <- sprintf("%s/PipelineBiota-Softwares", dirname(system.file(package = "PipelineBiota")))
  
  omicsdo_sof <- soft_directory
  print(omicsdo_sof)
  Bowtie2 <- sprintf('%s/Bowtie2/bowtie2-2.4.2-sra-linux-x86_64/bowtie2', soft_directory)
  
  #setwd(omicsdo_sof)
  if(!file.exists(Bowtie2)) {
  
      message("Installation of Bowtie2 will now begin. Check if all the required packages are downloaded.")
      
      dir.create(sprintf("%s/Bowtie2", soft_directory))
      
      download.file(
        #url = "https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.5.1/bowtie2-2.5.1-linux-x86_64.zip/download",
        url = "https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.5/bowtie2-2.4.5-linux-x86_64.zip/download",
        destfile = sprintf("%s/Bowtie2/bowtie2.zip", soft_directory), 
        mode = "wb"
      )
      unzip(sprintf("%s/Bowtie2/bowtie2.zip", soft_directory), exdir = sprintf("%s/Bowtie2", soft_directory))
      system(sprintf("export PATH=%s/Bowtie2/bowtie2-2.4.2-sra-linux-x86_64:$PATH", soft_directory))
      system(sprintf("chmod +x %s/Bowtie2/bowtie2-2.4.2-sra-linux-x86_64", soft_directory))
      
      
      #Escribo el path en el txt
      softwares <- readLines(sprintf("%s/PipelineBiota-Softwares/path_to_soft.txt", dirname(system.file(package = "OMICsdo"))))
      softwares_actualizado <- c(softwares, sprintf("Bowtie2 %s", Bowtie2))
      write(softwares_actualizado, file = sprintf("%s/PipelineBiota-Softwares/path_to_soft.txt", dirname(system.file(package = "OMICsdo"))))
      
      #system(sprintf("%s --version", Bowtie2))
      
  } else {
    message( "Bowtie2 was already downloaded")
  }
  
  #INDEX de ref con bowtie2:
  
  hg38_indexBowtie2 <- sprintf("%s/PipelineBiota-Softwares/HG38/index/HG38_index", dirname(system.file(package = "PipelineBiota")))
  dir_hg38_indexBowtie2 <- sprintf("%s/PipelineBiota-Softwares/HG38/index", dirname(system.file(package = "PipelineBiota")))
  if(!dir.exists(dir_hg38_indexBowtie2)) {
    #genero indice:
    out <- downloadHG38()
    hg38Fasta <- out[[2]]
    hg38_indexBowtie2 <- sprintf("%s/PipelineBiota-Softwares/HG38/index/HG38_index", dirname(system.file(package = "PipelineBiota")))
    system(sprintf("/home/juan/R/x86_64-pc-linux-gnu-library/4.1/PipelineBiota-Softwares/Bowtie2/bowtie2-2.4.2-sra-linux-x86_64/bowtie2-build %s HG38_index", hg38Fasta, hg38_indexBowtie2))
  } else {
    message("Bowtie2 Index was already generated")
  }
  
  return(list(Bowtie2, hg38_indexBowtie2))
}