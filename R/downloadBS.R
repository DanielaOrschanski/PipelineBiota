#' @title Download BS (Base Space)
#' @description downloads and installs the executable of bs for DRAGEN metagenomics processes by CLI Base Space aplication
#' @param dir_bs path where to store the executable Base Space
#' @export
#' @return path of the bs executable
downloadBS <- function(dir_bs) {
  #dir_bs <- "~/Daniela/Biota"

  if(!file.exists(sprintf("%s/bs", dir_bs))) {
    setwd(dir_bs)
    system("wget 'https://launch.basespace.illumina.com/CLI/latest/amd64-linux/bs'")
    system(sprintf("chmod u+x %s/bs", dir_bs))
    system(sprintf("export PATH=$PATH:%s/bs", dir_bs))
    system(sprintf("%s/bs auth", dir_bs))
  }

  bs_path  <- sprintf("%s/bs", dir_bs)
  return(bs_path)
}
