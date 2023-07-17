#!/usr/local/bin/Rscript
library(org.Hs.eg.db)

translate_ensemble <- function(x){
  xx <- as.data.frame(org.Hs.egENSEMBL2EG)
  x1 <- gsub("\\.[0-9]+$", "", rownames(x), perl = TRUE)
  
  x2 <- left_join(tibble(ensembl_id = x1, id = 1), 
                  xx %>% group_by(ensembl_id) %>% mutate(id = row_number()), 
                  by = c("ensembl_id", "id")) %>% pull(gene_id)
  x2[duplicated(x2)] <- NA
  x <- x[!is.na(x2),]
  rownames(x) <- x2[!is.na(x2)]
  return(x)
}


# summarize differential analyses
summarise_table <- function(res, test, group, pval_col, stat_col){
  require(tidyverse)
  my_summary <- res %>% group_by_at(group) %>% 
    summarise(n_diff = sum(!!rlang::sym(pval_col) < 0.05), 
              total = length(!!rlang::sym(pval_col)), 
              perc = sum(!!rlang::sym(pval_col) < 0.05)/total,
              up = sum(!!rlang::sym(pval_col) < 0.05 & !!rlang::sym(stat_col) > 0), 
              perc_up = up/total, 
              down = sum(!!rlang::sym(pval_col) < 0.05 & !!rlang::sym(stat_col) < 0), 
              perc_down = down/total) %>%
    add_column(.after = 2, diff_test = test)
  return(my_summary)
}


#Transcription Factor Target Gene Analysis (TFTEA) FUNCTION
uvGsa <- function (rankstat, annotation, p.adjust.method = "BH") {
  
  t0 <- proc.time ()
  
  ## variable names (colnames)
  statsnames <- colnames (rankstat)[1]
  if (is.null (statsnames)) statsnames <- "X"
  
  ## transform vector, matrix or data.frame keeping rownames
  rankstat <- as.matrix (rankstat)[,1] 
  
  ## LOOP over GOs
  annotation <- as.data.frame(annotation)
  GOs <- unique (annotation[,2])
  ncolRES <- 7 #3 + 4
  RES <- matrix (nrow = 0, ncol = ncolRES)
  separador <- "." 
  
  colnames (RES) <- c ("size", "conv", "error",                     #COLS: 1:3
                       paste ("LOR", statsnames, sep = separador),  #COLS: 4
                       paste ("sd",  statsnames, sep = separador),  #COLS: 5
                       paste ("z",   statsnames, sep = separador),  #COLS: 6
                       paste ("p",   statsnames, sep = separador))  #COLS: 7
  
  for (go in GOs) {
    
    geneswithgo <- annotation[annotation[,2] == go, 1]
    GObinario   <- as.numeric (names (rankstat) %in% geneswithgo)
    blocksize   <- sum (GObinario)
    
    G  <- try (glm (GObinario ~ rankstat, family = binomial))
    SG <- try (summary (G))
    
    try.error <- "try-error" %in% class (G) | "try-error" %in% class (SG)
    
    if (try.error) {
      RES <- rbind (RES, c(blocksize, NA, try.error, rep (NA, times = ncolRES - 3)))
    } else {
      RES <- rbind (RES, c(blocksize, G$converged, try.error, as.vector (SG$coefficients[-1,])))
    }
  }
  
  rownames (RES) <- GOs
  
  ## adjust p-values
  if (!is.null (p.adjust.method)) {
    adj.p.matrix <- as.matrix (p.adjust (RES[,7], method = p.adjust.method))
    colnames (adj.p.matrix) <- paste ("adj", statsnames, sep = separador)  #COLS: 8
    RES <- cbind (RES, adj.p.matrix)
  }
  
  RES <- as.data.frame (RES, stringsAsFactors = FALSE)
  
  t1 <- proc.time ()
  print ("timing in seconds")
  print (t1-t0)
  
  ## results
  return (RES)
}

## plot by FC
plot_tftea_by_fc <- function(rankstat, GObinario){
  breaks <- seq(min(rankstat), max(rankstat), 
                (max(rankstat) - min(rankstat))/50)
  freqs <- mat.or.vec(nr = length(breaks) - 1, nc = 3)
  for (i in seq(1, length(breaks) - 1)) {
    freqs[i,1] <- breaks[i]
    gb <- as.logical(GObinario)
    freqs[i,2] <- length(rankstat[gb][rankstat[gb] >= breaks[i] & rankstat[gb] < breaks[i+1]])
    freqs[i,3] <- length(rankstat[rankstat >= breaks[i] & rankstat < breaks[i+1]])
    
  }
  colnames(freqs) <- c("t", "tf", "bg")
  freqs[,"tf"] <- freqs[,"tf"]/max(freqs[,"tf"])
  freqs[,"bg"] <- freqs[,"bg"]/max(freqs[,"bg"])
  long_freq <- as.data.frame(freqs) %>% pivot_longer(c("tf", "bg"), names_to = "source")
  p <- ggplot(long_freq) + 
    geom_line(aes(x = t, y = value, color = source), size = 2)
  return(p)
}


## plot by rank
plot_tftea_by_rank <- function(rankstat, GObinario){
  breaks <- seq(1, length(rankstat), length(rankstat)/50)
  freqs <- mat.or.vec(nr = length(breaks) - 1, nc = 3)
  for (i in seq(1, length(breaks) - 1)) {
    freqs[i,1] <- breaks[i]
    gb <- as.logical(GObinario)[breaks[i]:breaks[i+1]]
    freqs[i,2] <- length(rankstat[breaks[i]:breaks[i+1]][gb])
    freqs[i,3] <- length(rankstat[breaks[i]:breaks[i+1]])
    
  }
  colnames(freqs) <- c("rank", "tf", "bg")
  freqs[,"tf"] <- freqs[,"tf"]/max(freqs[,"tf"])
  freqs[,"bg"] <- freqs[,"bg"]/max(freqs[,"bg"])
  long_freq <- as.data.frame(freqs) %>% pivot_longer(c("tf", "bg"), names_to = "source")
  p <- ggplot(long_freq) + 
    geom_line(aes(x = rank, y = value, color = source), size = 2)
  return(p)
}

