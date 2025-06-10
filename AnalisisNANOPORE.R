path_fastq <- "~/DatosDiscoExterno/BiotalifeRawSeqs/20221006/fastq5_pass/barcode12/AIE980_pass_barcode12_44a53ea4_64.fastq.gz"
lines <- readLines(con = gzfile(path_fastq), n = 4)
cat(lines, sep = "\n")

#@04ddcca3-3105-409a-a203-efed26798974
#runid=44a53ea41c39071fb46bf7d741b989391590906e
#read=38163 ch=68 start_time=2022-10-07T18:15:09.643088+00:00
#flow_cell_id=AIE980 protocol_group_id=2022-10-06-CiveServ-16S
#sample_id=2022-10-06-CiveServ-16S
#barcode=barcode12 barcode_alias=barcode12 parent_read_id=04ddcca3-3105-409a-a203-efed26798974
#basecall_model_version_id=2021-05-17_dna_r9.4.1_minion_96_29d8704b"

# 16S: 1600 bp por lectura deberian haber
# Deberian ser los mismos IDs que piel


