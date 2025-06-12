setwd("/media/4tb2/Daniela/Biota/Muestras")

pipeline_biota_path <- "/media/4tb2/Daniela/Biota/PipelineBiota-master"
muestras_biota_path <- "/media/4tb2/Daniela/Biota/Muestras"

library(PipelineBiota)

RunPipelineIndividual(patient_dir = "/media/4tb2/Daniela/Biota/Muestras/269", runDRAGEN = FALSE)
