# Leonardo Corvalan
# automate DNA barcoder
# 09/04/2026 -last update

options(show.error.locations = TRUE)
options(error = traceback)

# Pacotes
library(Biostrings)

# receber input

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  cat("Uso: Rscript script.r <arquivo.tsv>  <output.tsv>\n")
  quit(status = 2)
}

input <- args[1]
output <- args[2]

if (!file.exists(input)) {
  stop(paste("Arquivo não encontrado:", arquivo))
}

# parametros 

cutoff <- 30
window <- 4 
step <- 2
AA_lib <- "models.aa.fasta"
Bold_db <- "/media/lgbio-nas1/lcorvalan/Sanger_barcode/BOLD_db"

#input <- "Ab1.csv"
# criando AA db para o blast

# cria o db # mover para antes do loop
#makedb <- paste0("makeblastdb -in ", AA_lib, " -dbtype prot -out tmp_lib_aa")

#system(makedb)

# entrada 
# aqui eu crio um loop para roda para toda lista

ab1_tbl <- read.csv(input)

tbl_final <- data.frame()

for (i in 1:length(ab1_tbl[,1])) {
  Fwd_input <- ab1_tbl$Fwd_file[i]
  Rev_input <- ab1_tbl$Rev_file[i]
  output_prefix<- ab1_tbl$Prefix[i]

  ################################################################################
  ### step1: Convertendo dados do ab1 e gerendo consenso #########################
  ################################################################################
  
  # Ler arquivo 
  # os nomes nao podem conter ()
  # editando nome
  
  #rename <- "for f in *.ab1; do mv \"$f\" \"$(echo $f | tr '()' '--')\"; done"
  #system(rename)
  
  # dados ab1
  ### Fwd 
  # Fwd
  # Qualidade do sequenciamento
  # ab1 to fq
  
  Fwd_ab1_out <- paste0(output_prefix,"_F_raw.fq")
  ab1_fq <- paste0("docker run --rm -v $(pwd):/data geargenomics/tracy tracy basecall -f fastq -o /data/",Fwd_ab1_out, " /data/", Fwd_input)
  
  #ab1_fq <- "docker run --rm -v $(pwd):/data geargenomics/tracy tracy basecall -f fastq -o /data/PAV296_rbcLB-F_A01_BCPlan-06-S.fq /data/PAV296_rbcLB-F_A01_BCPlan-06-S.ab1"
  
  system(ab1_fq)
  
  # plota
  #fastqc -t --noextract PAV296_rbcLB-F_A01_BCPlan-06-S.fq
  
  # Informacoes fwd
  seq_info_q <- paste0("/usr/bin/python3 seq_q.py ", Fwd_ab1_out," ", cutoff, " " , window," ", step, " > Fwd_inf.tsv")
  #seq_info_q <-"/usr/bin/python3 seq_q.py PAV296_rbcLB-F_A01_BCPlan-06-S.fq 30 4 2 > Fwd_inf.tsv"
  print(seq_info_q)
  
  system(seq_info_q)
  
  Fwd_raw <- read.table("Fwd_inf.tsv", sep = "\t", header = T, stringsAsFactors = FALSE)
  
  Fwd_ab1_trim_out <- paste0(output_prefix,"_F_trim.fq")
  print(paste("raw len:", Fwd_raw$read_len))
  
  if (Fwd_raw$read_len>200) { # testando tamanho
    Fwd_raw$last_before[1] <- ifelse((Fwd_raw$last_before[1]=="None" | Fwd_raw$last_before[1]<= 10), yes = 10, no = Fwd_raw$last_before[1])
    Fwd_raw$first_after[1] <- ifelse((Fwd_raw$first_after[1]=="None" | Fwd_raw$first_after[1]<= 0), yes = Fwd_raw$read_len[1], no = Fwd_raw$first_after[1])
    
    fwd_trim <- paste0("/usr/local/bin/Conda_apps/seqkit subseq -r ", Fwd_raw$last_before,":",Fwd_raw$first_after, " ",Fwd_ab1_out, " > ", Fwd_ab1_trim_out)
    
    print(fwd_trim)
    system(fwd_trim)
    print("ok")
    
    # qualidade pos filtro
    
    seq_info_q <- paste0("/usr/bin/python3 seq_q.py ", Fwd_ab1_trim_out," ", cutoff, " " , window," ", step, " > Fwd_inf_trim.tsv")
    print(seq_info_q)
    #seq_info_q <- "/usr/bin/python3 seq_q.py PAV296_rbcLB-F.fq 30 4 2 > Fwd_inf_trim.tsv"
    system(seq_info_q)
    
    Fwd_trim <- read.table("Fwd_inf_trim.tsv", sep = "\t", header = T, stringsAsFactors = FALSE)
    
    print(Fwd_trim)
    
    if (Fwd_trim$read_len>200) {
      Fwd_status <- "ok"
      
    } else {
      Fwd_status <- " trim short sequence"
    }
  } else {
    Fwd_trim <- data.frame(
      id = "primary",
      read_len = "short sequence",
      mean_q = "short sequence",
      last_before = "short sequence",
      first_after = "short sequence",
      row.names = NULL
    )
    Fwd_status <- "short sequence"
  }
  
  print(paste(output_prefix,"Fwd status",Fwd_status ))
  
  # cortando .fq
  
  
  # se nao tiver que corta no inicio

  # dados ab1
  ### Rev
  # Qualidade do sequenciamento
  # ab1 to fq
  Rev_ab1_out <- paste0(output_prefix,"_R_raw.fq")
  ab1_fq <- paste0("docker run --rm -v $(pwd):/data geargenomics/tracy tracy basecall -f fastq -o /data/",Rev_ab1_out, " /data/", Rev_input)
  
  #ab1_fq <- "docker run --rm -v $(pwd):/data  geargenomics/tracy tracy basecall -f fastq -o /data/PAV296_rbcLB-R_A01_BCPlan-07-S.fq  /data/PAV296_rbcLB-R_A01_BCPlan-07-S.ab1"
  
  system(ab1_fq)
  
  # info rev
  Rev_ab1_trim_out <- paste0(output_prefix,"_R_trim.fq")
  print(Rev_ab1_trim_out)
  seq_info_q <- paste0("/usr/bin/python3 seq_q.py ", Rev_ab1_out," ", cutoff, " " , window," ", step, " > Rev_inf.tsv")
  
  system(seq_info_q)
    
  Rev_raw <- read.table("Rev_inf.tsv", sep = "\t", header = T, stringsAsFactors = FALSE)
  
  if (Rev_raw$read_len>200) { # testando tamanho bruuto
    
    # cortando .fq
    Rev_ab1_trim_out <- paste0(output_prefix,"_R_trim.fq")
    # se nao tiver que corta no inicio
    Rev_raw$last_before[1] <- ifelse((Rev_raw$last_before[1]=="None" | Rev_raw$last_before[1]<= 10), yes = 10, no = Rev_raw$last_before[1])
    Rev_raw$first_after[1] <- ifelse((Rev_raw$first_after[1]=="None" | Rev_raw$first_after[1] <= 0), yes = Rev_raw$read_len[1], no = Rev_raw$first_after[1])
    
    
    Rev_trim <- paste0("/usr/local/bin/Conda_apps/seqkit subseq -r ", Rev_raw$last_before,":",Rev_raw$first_after, " ",Rev_ab1_out, " > ", Rev_ab1_trim_out)
    system(Rev_trim)
    
    # qualidade pos filtro
    seq_info_q <- paste0("/usr/bin/python3 seq_q.py ", Rev_ab1_trim_out," ", cutoff, " " , window," ", step, " > Rev_inf_trim.tsv")
    
    system(seq_info_q)
    
    Rev_trim <- read.table("Rev_inf_trim.tsv", sep = "\t", header = T, stringsAsFactors = FALSE)
    
    if (Rev_trim$read_len>200) { #tamanho apos filtro
      Rev_status <- "ok"
      
    } else {
      Rev_status <- " trim short sequence"
      }
    } else { # quando o tamanho da sequencia e menor antes
      Rev_trim <- data.frame(
        id = "primary",
        read_len = "short sequence",
        mean_q = "short sequence",
        last_before = "short sequence",
        first_after = "short sequence",
        row.names = NULL
      )
      Rev_status <- "short sequence"
    }
    
    print(paste(output_prefix,"Rev status",Rev_status ))
  ################################################################################
  ### step2: Merge fwd and Rev  ##################################################
  ################################################################################
  
  if(Rev_status== "ok" & Fwd_status== "ok") { # teste para ver se pode juntar
    print(paste(output_prefix,"Tamanho ok merging"))
    cut <- as.numeric(Fwd_raw$read_len)-as.numeric(Fwd_raw$first_after)
    
    # pegando tamanho 
    resultado <- ifelse(
      is.na(suppressWarnings(as.numeric(Rev_raw$read_len))) |
        is.na(suppressWarnings(as.numeric(Rev_raw$first_after))),
      0,
      as.numeric(Rev_raw$read_len) - as.numeric(Rev_raw$first_after)
    )
    
    merge <- paste0("docker run --rm -v $(pwd):/data geargenomics/tracy tracy consensus -t 0 -q ",
                    Fwd_raw$last_before," -u ", cut, 
                    " -r ", Rev_raw$last_before, " -s ", resultado,
                    " -o /data/", output_prefix, " /data/", Fwd_input, " /data/", Rev_input)
    
    system(merge)
    merge_status<- "ok" 
  } else {
    merge_status<- "short sequence"
  }
  
  print(paste0("Merge status: ", merge_status))
  merged_file <- paste0(output_prefix,".fa")
  
  ################################################################################
  ### step3: Frame and stop codom  ###############################################
  ################################################################################
  
  # identificando frame de leitura
  # blastX com sequencia de aminoacido de rbcL de Arabidopsis 
  # https://www.ncbi.nlm.nih.gov/protein/NP_051067.1/
  
  # aqui vc pode alterar o seus base dade dados de AA
  if (merge_status == "ok") {
    # blast
    
    blastx <- paste0('blastx -query ',merged_file, ' -db tmp_lib_aa -out Blast_test.txt -max_target_seqs 1 -evalue 1e-25 -outfmt "6 qseqid sseqid pident length evalue bitscore qstart qend sstart send qframe"')
    
    system(blastx)
    
    # lendo blast
    
    blast_aa_results <- read.table("Blast_test.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    
    # add headers
    
    colnames(blast_aa_results) <- c("qseqid", "sseqid", "pident" ,"length" ,"evalue" ,"bitscore" ,"qstart" ,"qend", "sstart", "send", "qframe")
    blast_aa_results <- blast_aa_results[order(-blast_aa_results$bitscore), ]
    blast_aa_results <- blast_aa_results[1, ]
    
    frame <-  blast_aa_results$qframe
    
    # ler sequencia 
    seq_tmp <- readDNAStringSet(merged_file)
    
    # Pegando sentido correto
    if (frame > 0) {
      seq_tmp <- seq_tmp
    } else {
      seq_tmp <- reverseComplement(seq_tmp)
    }
    # corrigindo fram
    seq_tmp <- subseq(seq_tmp, start= abs(frame))
    
    # transformar em AA para ver stopcondons
    
    # genetic code
    GENETIC_CODE_Chloroplast <- getGeneticCode("Bacterial, Archaeal and Plant Plastid",full.search=T)
    
    # Traduzindo
    seq_tmp_aa <- translate(seq_tmp, genetic.code = GENETIC_CODE_Chloroplast, if.fuzzy.codon= "X")
    
    
    # Procurar stopcodon
    
    stop_codon_site <- function(x) {
      stop_positions <- which(unlist(strsplit(x, "")) == "*") # pegas as posicoes do stopcodon
      if (length(stop_positions) == 0) { # se nao tiver stop
        return("No_stop_codon")
      } else if (length(stop_positions) == 1) { # testar se stop codon esta apenas no final
        if (stop_positions == nchar(x)) {
          return("Final_stop_codon")
        } else {
          return("Middle_stop_codon")
        }
      } else {
        return("Mult_stop_condons")
      }
    }
    
    
    stop_codon <- stop_codon_site(as.character(seq_tmp_aa))
    
    # numero de N ou IUPAC
    
    Nucleotide_ATCG <-sapply(gregexpr("[AaTtCcGg]",seq_tmp), function(x) sum(x > 0))
    
    final_len <- width(seq_tmp)
    
    # IUPAC
    N_iupac <- final_len-Nucleotide_ATCG
    P_iupac <- N_iupac/final_len
    
  } else {
    blast_aa_results <- data.frame("qseqid"="NA", "sseqid" = "NA", "pident"= "NA",
                                   "length"= "NA", "evalue" = "NA" ,"bitscore"= "NA",
                                   "qstart"= "NA", "qend"= "NA", "sstart"= "NA",
                                   "send" = "NA", "qframe"= "NA",row.names = NULL)
    frame <-  blast_aa_results$qframe
    seq_tmp <- data.frame( Consensus = "NA")
    #seq_tmp$Consensus <- "NA"
    seq_tmp_aa <- data.frame( Consensus = "NA")
    stop_codon <- "NA"
    Nucleotide_ATCG <- "NA"
    final_len <- "NA"
    N_iupac <- "NA"
    P_iupac <- "NA"
  }
  
  
  
  ################################################################################
  ### step4: Blast in BOLD DB  ###################################################
  ################################################################################
  
  # criar db do bold
  # vou deixar esse db construido
  #make_db <- "makeblastdb -in bold.fa -dbtype nucl -out BOLD_db"
  if (merge_status == "ok") {
    
    #system(make_db)
    out_final_seq <- paste0(output_prefix,"_final.fasta")
    writeXStringSet(seq_tmp, out_final_seq)
    
    # blastando com bold
    blast_tmp <-  paste0("blastn -query ",out_final_seq,
           " -db " , Bold_db, " -out seq_bold.txt -max_target_seqs 1 -evalue 1e-25 -outfmt 6")
    
    system(blast_tmp)
    
    # carregando
    
    blast_bold <- read.table("seq_bold.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    
    # add headers
    
    colnames(blast_bold) <- c("qseqid", "sseqid" ,"pident", "length", "mismatch", "gapopen", 
                              "qstart", "qend", "sstart", "send", "evalue", "bitscore")
  } else {
    blast_bold <- data.frame("qseqid"= "NA", "sseqid"= "NA" ,"pident"= "NA", 
                             "length" = "NA", "mismatch"= "NA", "gapopen" = "NA", 
                      "qstart" = "NA", "qend" = "NA", "sstart" = "NA", "send" = "NA",
                      "evalue"= "NA", "bitscore"= "NA",
                      row.names = NULL)
  }
  
  
  ################################################################################
  ### step5: Result compile   ####################################################
  ################################################################################
  
  # ajustando cabecalho
  
  
  colnames(Fwd_trim) <- paste0("Fwd_trim_", colnames(Fwd_trim))
  
  colnames(Fwd_raw) <- paste0("Fwd_raw_", colnames(Fwd_raw))
  
  colnames(Rev_trim) <- paste0("Rev_trim_", colnames(Rev_trim))
  
  colnames(Rev_raw) <- paste0("Rev_raw_", colnames(Rev_raw))
  
  # data frame merged data
  df_mergseq<- data.frame(
    prefix = output_prefix,
    Merged_seq = as.character(seq_tmp$Consensus),
    Merged_length = final_len,
    Merged_n_iupac = N_iupac,
    Merged_stop_codon = stop_codon,
    Merge_aa_seq = as.character(seq_tmp_aa$Consensus),
    row.names = NULL
  )
  
  
  # add blas
  print(paste0("linha: ", i, " ok"))
  df_tmp <- cbind(Fwd_raw, Fwd_trim[,2:3], Rev_raw, Rev_trim[,2:3], df_mergseq, blast_bold)
  
  tbl_final <- rbind(tbl_final, df_tmp)
}

write.csv(tbl_final, output)
