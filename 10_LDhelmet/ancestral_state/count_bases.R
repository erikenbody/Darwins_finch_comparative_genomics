require(Biostrings)


STF5_HiC <- readDNAStringSet("/Users/erikenbody/Documents/upload_genome_ENA/Camarhynchus_parvulus_V1.1.fasta")
alphabetFrequency(STF5_HiC, baseOnly=TRUE, collapse=TRUE)
