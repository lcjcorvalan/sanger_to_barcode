# Leonardo Corvalan
# automate DNA barcoder - pre
# 09/04/2026 -last update

# editar nomes
direc <- "Ab1_file"
arquivos_ab1 <- list.files(path = direc , pattern = "\\.ab1$", full.names = TRUE)

arquivos_fwd <- arquivos_ab1[grepl("rbcLB-F", arquivos_ab1)]
arquivos_rev <- arquivos_ab1[grepl("rbcLB-R", arquivos_ab1)]


# Os arquivos tem tamanho diferentes

arquivos_fwd <- arquivos_fwd[1:length(arquivos_rev)]
arquivos_rev

# codigo das especies

nomes_antes <- sub("-.*", "", arquivos_fwd)

# criando a tabela

df_ab1_file <-data.frame(Fwd_file=arquivos_fwd,
           Rev_file=arquivos_rev,
           Prefix= nomes_antes,
           row.names = NULL)

write.csv(file = "Ab1.csv", df_ab1_file, row.names = F)
?write.csv
