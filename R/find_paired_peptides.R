
find_paired_peptides = function(){
  # Finding paired peptide ids
  paired_peptides = mspms::processed_qf %>%  mspms::mspms_tidy("peptides_norm") %>% 
    dplyr::select(quantCols,library_id,peptide,cleavage_seq,peptides_norm) %>% 
    tidyr::pivot_wider(names_from = peptide,
                       values_from = peptides_norm) %>% 
    dplyr::mutate(non_na_count = rowSums(!is.na(.))-3) %>% 
    dplyr::select(quantCols,cleavage_seq,library_id,non_na_count) %>%
    dplyr::filter(non_na_count > 1) %>% 
    dplyr::mutate(library_id_cleavage_seq = paste0(library_id,".",cleavage_seq)) %>% 
    dplyr::pull(library_id_cleavage_seq)
  
  # Counting how many paired peptides were detected
  n_paired = length(unique(paired_peptides))
  
  
  ## Looking a the correlation of Processed values
 pn =  mspms::processed_qf %>%  mspms::mspms_tidy("peptides_norm") %>% 
   dplyr::mutate(library_id_cleavage_seq = paste0(library_id,".",cleavage_seq)) %>% 
   filter(library_id_cleavage_seq %in% paired_peptides) %>% 
   mutate(cleavage_side = case_when(grepl("^._.*",peptide) ~ "cterm",
                                    TRUE ~ "nterm")) %>% 
   dplyr::select(-peptide) %>% 
   tidyr::pivot_wider(names_from = cleavage_side,
                      values_from = peptides_norm)
 
 library(ggplot2)
 
 lm = lm(data = pn, nterm ~cterm)
 summary(lm)
 p1 = ggplot(aes(x = nterm, y = cterm), data = pn) + 
   geom_point() + 
   theme_bw() + 
   xlab("N-terminal peptide") + 
   ylab("C-terminal peptide") + 
   ggtitle(paste0("Paired peptides detected: ",n_paired))+
   geom_smooth(method = "lm")+
   geom_abline(slope = 1, intercept = 0,linetype = "dashed")
    
p1

# Looking at how well the signicance matches
## Looking a the correlation of Processed values
sig=  mspms::processed_qf %>% mspms::limma_stats() %>% 
  dplyr::mutate(library_id_cleavage_seq = paste0(library_id,".",cleavage_seq)) %>% 
  dplyr::filter(library_id_cleavage_seq %in% paired_peptides) %>% 
  dplyr::mutate(cleavage_side = case_when(grepl("^._.*",peptide) ~ "cterm",
                                   TRUE ~ "nterm")) %>% 
  dplyr::select(-peptide) %>% 
  dplyr::select(library_id,comparison,cleavage_seq,p.adj,cleavage_side) %>% 
  tidyr::pivot_wider(names_from = cleavage_side,
                     values_from = p.adj)

p2 = ggplot(aes(x = -log10(nterm), y = -log10(cterm)), data = sig) + 
  geom_point() + 
  theme_bw() + 
  xlab("N-terminal peptide P.adj") + 
  ylab("C-terminal peptide P.adj") + 
  ggtitle(paste0("Paired peptides detected: ",n_paired))+
  geom_smooth(method = "lm")+
  geom_vline(xintercept = -log10(0.05),linetype = "dashed")+
  geom_hline(yintercept = -log10(0.05),linetype = "dashed")

p2

## Looking a the correlation of Processed values
log2fc=  mspms::processed_qf %>% mspms::limma_stats() %>% 
  dplyr::mutate(library_id_cleavage_seq = paste0(library_id,".",cleavage_seq)) %>% 
  dplyr::filter(library_id_cleavage_seq %in% paired_peptides) %>% 
  dplyr::mutate(cleavage_side = case_when(grepl("^._.*",peptide) ~ "cterm",
                                          TRUE ~ "nterm")) %>% 
  filter(time == 90) %>% 
  mutate(length = nchar(peptide)-2) %>% 
  dplyr::select(-peptide) %>% 
  dplyr::select(library_id,comparison,cleavage_seq,p.adj,log2fc,length,cleavage_side) %>% 
  tidyr::pivot_wider(names_from = cleavage_side,
                     values_from = c(p.adj,log2fc,length)) %>% 
  mutate(length_dif = abs(length_nterm - length_cterm))


p3 = ggplot(aes(x = -log10(p.adj_nterm), y = -log10(p.adj_cterm),color = length_dif), data = log2fc) + 
  geom_point() + 
  theme_bw() + 
  ggtitle(paste0("Paired peptides detected: ",n_paired))+
  geom_smooth(method = "lm")+
  geom_vline(xintercept = 3,linetype = "dashed")+
  geom_hline(yintercept = 3,linetype = "dashed")

p3


}

