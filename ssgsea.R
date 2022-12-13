library(IOBR)
ssgsea.array <- function(x){
  ciber<-deconvo_tme(eset = x, method = "cibersort", arrays = TRUE, perm = 200)

  sigHallmark<-calculate_sig_score(pdata           = NULL,
                                   eset            = x,
                                   signature       = hallmark,
                                   method          = "ssgsea",
                                   mini_gene_count = 5,
                                   parallel.size = 10)
  
  sigKegg <-calculate_sig_score(pdata           = NULL,
                                eset            = x,
                                signature       = kegg,
                                method          = "ssgsea",
                                mini_gene_count = 5,
                                parallel.size = 10)
  ssgsea_result <- estimate %>% 
    inner_join(ciber,by="ID") %>%
    inner_join(sigHallmark, by="ID") %>%
    inner_join(sigKegg,by="ID")
  return(ssgsea_result)
}
ssgsea_combat <- ssgsea.array(combined.array.combat)
