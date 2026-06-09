# Leonardo Corvalan
# automate DNA barcoder - pre
# 18/05/2026 -last update

# rcbL A
# editar nomes
direc_F <- "Ab1_file/rbcL_A_F"
direc_R <- "Ab1_file/rbcL_A_R"

arquivos_ab1_F <- list.files(path = direc_F , pattern = "\\.ab1$", full.names = TRUE)
arquivos_ab1_R <- list.files(path = direc_R , pattern = "\\.ab1$", full.names = TRUE)

arquivos_fwd <- arquivos_ab1_F[grepl("A_F", arquivos_ab1_F)]
arquivos_rev <- arquivos_ab1_R[grepl("A_R", arquivos_ab1_R)]


# Os arquivos tem tamanho diferentes

length(arquivos_fwd) == length(arquivos_rev)

# remover tudo antes do 
nomes_antes <- sub("Ab1_file/rbcL_A_F/", "", arquivos_fwd)

nomes_depois <- sub("_.*", "", nomes_antes)

# cria df vazio
df_ab1_file <- data.frame(
  Fwd_file = character(),
  Rev_file = character(),
  Prefix   = character()
)

# loop
nomes_depois <- paste0(nomes_depois,"_")

for (id in nomes_depois) {
  
  tmp_fwd <- arquivos_fwd[grepl(id, arquivos_fwd)]
  tmp_rev <- arquivos_rev[grepl(id, arquivos_rev)]
  
  # pula IDs com números diferentes
  if(length(tmp_fwd) != length(tmp_rev)) {
    next
  }
  
  tmp_prefix <- rep(paste0("rbcL_A/", id), length(tmp_fwd))
  
  df_tmp <- data.frame(
    Fwd_file = tmp_fwd,
    Rev_file = tmp_rev,
    Prefix   = tmp_prefix
  )
  
  df_ab1_file <- rbind(df_ab1_file, df_tmp)
}
# codigo das especies

View(df_ab1_file)
write.csv(file = "rbcL_A_Ab1.csv", df_ab1_file, row.names = F)


# rcbL B
# editar nomes
direc_F <- "Ab1_file/rbcL_B_F"
direc_R <- "Ab1_file/rbcL_B_R"

arquivos_ab1_F <- list.files(path = direc_F , pattern = "\\.ab1$", full.names = TRUE)
arquivos_ab1_R <- list.files(path = direc_R , pattern = "\\.ab1$", full.names = TRUE)

arquivos_fwd <- arquivos_ab1_F[grepl("B_F", arquivos_ab1_F)]
arquivos_rev <- arquivos_ab1_R[grepl("B_R", arquivos_ab1_R)]


# Os arquivos tem tamanho diferentes

length(arquivos_fwd) == length(arquivos_rev)

# remover tudo antes do 
nomes_antes <- sub("Ab1_file/rbcL_B_F/", "", arquivos_fwd)

nomes_depois <- sub("_.*", "", nomes_antes)

nomes_depois <- paste0(nomes_depois,"_")

# cria df vazio
df_ab1_file <- data.frame(
  Fwd_file = character(),
  Rev_file = character(),
  Prefix   = character()
)

# loop

for (id in nomes_depois) {
  
  tmp_fwd <- arquivos_fwd[grepl(id, arquivos_fwd)]
  tmp_rev <- arquivos_rev[grepl(id, arquivos_rev)]
  
  # pula IDs com números diferentes
  if(length(tmp_fwd) != length(tmp_rev)) {
    next
  }
  
  tmp_prefix <- rep(paste0("rbcL_B/", id), length(tmp_fwd))
  
  df_tmp <- data.frame(
    Fwd_file = tmp_fwd,
    Rev_file = tmp_rev,
    Prefix   = tmp_prefix
  )
  
  df_ab1_file <- rbind(df_ab1_file, df_tmp)
}
# codigo das especies


write.csv(file = "rbcL_B_Ab1.csv", df_ab1_file, row.names = F)
?write.csv
