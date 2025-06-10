#' @title downloadHG38
#' @description Downloads the FASTA and the annotation of the genome reference version HG38.
#' @return Paths of FASTA and GTF(annotation) from the genome reference.
#' @export
downloadHG38 <- function() {
  message("ESTOY EN DOWNLOAD HG38")
  #soft_directory <- sprintf("%s/OMICsdoSof", Sys.getenv('R_LIBS_USER'))
  soft_directory <- sprintf("%s/PipelineBiota-Softwares", dirname(system.file(package = "PipelineBiota")))


  #Checks if Annotation is already downloaded --------------------
  #softwares <- readLines(sprintf("%s/OMICsdoSof/path_to_soft.txt", Sys.getenv('R_LIBS_USER')))
  softwares <- readLines(sprintf("%s/PipelineBiota-Softwares/path_to_soft.txt", dirname(system.file(package = "PipelineBiota"))))

  linea_software <- grep("HG38Annotation", softwares, ignore.case = TRUE, value = TRUE)

  if(length(nchar(linea_software)) == 0) {
    #download Annotation
    message("HG38 annotation will be downloaded")
    dir.create(sprintf("%s/HG38", soft_directory))
    URL <- "https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz"

    dir <- sprintf("%s/HG38", soft_directory)
    setwd(dir)
    system2("wget", args = c(URL, "-P", dir), wait = TRUE, stdout = NULL, stderr = NULL)
    AnnotationHG38 <- sprintf("%s/HG38/Homo_sapiens.GRCh38.110.gtf.gz", soft_directory)

    #gunzip(AnnotationHG38, destname = gsub("[.]gz$", "", AnnotationHG38), overwrite = FALSE, remove = TRUE)
    output_file <- gsub("[.]gz$", "", AnnotationHG38)
    system2("gunzip", args = c("-k", "-c", AnnotationHG38), stdout = output_file)

    AnnotationHG38 <<- sprintf("%s/HG38/Homo_sapiens.GRCh38.110.gtf", soft_directory)

    #Writes down the paths in the txt
    softwares <- readLines(sprintf("%s/PipelineBiota-Softwares/path_to_soft.txt", dirname(system.file(package = "PipelineBiota"))))

    softwares_actualizado <- c(softwares, sprintf("HG38Annotation %s", AnnotationHG38))
    write(softwares_actualizado, file = sprintf("%s/PipelineBiota-Softwares/path_to_soft.txt", dirname(system.file(package = "PipelineBiota"))))

  } else {
    message("The annotation for HG38 was already downloaded")
    AnnotationHG38 <<- strsplit(linea_software, " ")[[1]][[2]]
  }

  #Checks if FASTA is already downloaded --------------------------------------------------------
  linea_software <- grep("HG38FASTA", softwares, ignore.case = TRUE, value = TRUE)

  if(length(nchar(linea_software)) == 0) {
    message("HG38 FASTA will be downloaded")
    #download Fasta HG38
   
    URL <- "https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz"
    system2("wget", args = c(URL, "-P", dir), wait = TRUE, stdout = NULL, stderr = NULL)
    
    FastaHG38 <- sprintf("%s/HG38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz", soft_directory)
    output_file <- gsub("[.]gz$", "", FastaHG38)
    system2("gunzip", args = c("-k", "-c", FastaHG38), stdout = output_file)

    FastaHG38 <<- sprintf("%s/HG38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa", soft_directory)
  
    #Indexar con bwa la referencia:
    
    #message("The reference will be indexed by BWA")
    #system(sprintf("bowtie2-build index %s", FastaHG38))
    #system(sprintf("bwa index %s", FastaHG38))
    
    #Writes down the paths in the txt
    softwares <- readLines(sprintf("%s/PipelineBiota-Softwares/path_to_soft.txt", dirname(system.file(package = "PipelineBiota"))))
    softwares_actualizado <- c(softwares, sprintf("HG38FASTA %s", FastaHG38))
    write(softwares_actualizado, file = sprintf("%s/PipelineBiota-Softwares/path_to_soft.txt", dirname(system.file(package = "PipelineBiota"))))

  } else {
    message("The fasta file for HG38 was already downloaded")
    FastaHG38 <<- strsplit(linea_software, " ")[[1]][[2]]
  }


 #index_dir_STAR <- indexRefSTAR(AnnotationHG38, FastaHG38)

  #return(c(AnnotationHG38, FastaHG38, index_dir_STAR))
  return(c(AnnotationHG38, FastaHG38))

}
