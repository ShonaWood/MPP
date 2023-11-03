# libraries

library(tidyverse)
library(fuzzyjoin)

### This is how the tfbs_per_gene_promoter.tsv is made - this does not need to be run again if you have the tsv####

all.fimo.folders = dir("fimo_out/", full.names = T)

# keep only vertebrate core
vertebrate.core <- readLines("jaspar_core_tfbs_vertebrates.txt") %>% grep(">", x = ., value = T) %>% 
    gsub(">", "", .) %>% gsub("\t.*", "", .)


all.fimo.folders <- all.fimo.folders[gsub("fimo_out/", "", all.fimo.folders) %in% vertebrate.core]

## run supercode
all.fimo.files = sapply(all.fimo.folders, function(i) {
                        grep("\\.tsv", dir(i, full.names = T),  value = T)
                        })
names(all.fimo.files) <- gsub("fimo_out\\/\\/", "", names(all.fimo.files))

genes <- prom %>% distinct(gene_id)


fimo_promoters <- lapply(all.fimo.files, function(i){
                      print("\n\n")
                      print(i)
                      tfbs <-  read_tsv(i[[1]])
                      if(nrow(tfbs)==2) {return(NULL)}
                      prom_tfbs <- prom %>% genome_inner_join(tfbs,
                                                 by=c('sequence_name','start_promoter'='start','end_promoter'='stop')) %>%
                        mutate(tfbsID = paste(sequence_name.x, start, stop)) %>% # NOW sort by pvalue and remove identical matches on different strands.  
                        filter(start_promoter<=start & end_promoter>=stop) %>% 
                        arrange(`q-value`) %>% select(-sequence_name.y) %>% rename(sequence_name = sequence_name.x) %>% 
                        distinct(tfbsID, .keep_all = T) %>% 
                        filter(`p-value`< 0.05) %>% 
                        group_by(gene_id) %>% # NOTE: you need to make a choice about how to handle multiple promoter regions from multiple transcrips
                        summarize(motif_count = n()) 
                      
                      
                        id <- gsub("fimo_out\\/\\/|\\/fimo.tsv", "", i)
                        colnames(prom_tfbs)[2] <- gsub("fimo_out\\/\\/|\\/fimo.tsv", "", i)
                        
                        names <- tfbs %>% select(1:2) %>% slice(1) # getting the human readable jaspar motif name
                
                        # joining promoter count info for all genes with jaspar names
                        prom_tfbs_final <- genes %>% left_join(prom_tfbs, by="gene_id")
                        colnames(prom_tfbs_final)[2] <- paste0(colnames(prom_tfbs_final)[2], "_", names$motif_alt_id)
                        prom_tfbs_final  
                        }
                      )


names(fimo_promoters) <- sapply(fimo_promoters, function(i) colnames(i)[2])

length(fimo_promoters) # this many total fimo motifs  
length(compact(fimo_promoters)) # this many fimo motifs have been found ... 


## make final table 
fimo_promoters <- compact(fimo_promoters) # removing NULL values from list
fimo_promoters_final <- cbind(fimo_promoters[[1]][-2], do.call(cbind, lapply(fimo_promoters, `[[`, 2)))

#### Replace NAs with zeros
fimo_promoters_final <- as_tibble(fimo_promoters_final)  
fimo_promoters_final_part2 <- fimo_promoters_final %>% select(-1) %>% replace(is.na(.), 0)
fimo_promoters_final_part1 <- fimo_promoters_final %>% select(1) 
fimo_promoters_final <- as_tibble(cbind(fimo_promoters_final_part1, fimo_promoters_final_part2))

write_tsv(fimo_promoters_final, "tfbs_per_gene_promoter.tsv")


##############START HERE if you have the BED file for promoters and TSV file for tfbs per promoter###############################
# read in bed with promoters

library(tidyverse)

prom <- read_tsv("sib-ham-polished-partials-included-fix.gtf_promoters_u2000_d200.bed")
prom <- prom %>% rename(sequence_name = seqname)

tfbs <- read_tsv("tfbs_per_gene_promoter.tsv")

plot(colSums(tfbs[,-1])) # number of TFBS per gene


interestingGenes <- read_delim("Gene-list.txt", delim = " ") #%>% select(2) %>% rename(gene_id = logFC)
# weird lacking "E" in front of GTF gene names..
interestingGenes$gene_id <- substr(interestingGenes$gene_id, 2, 100)

otherGenes <- genes %>% filter(!gene_id %in% interestingGenes$gene_id) #tfbs$transcript_id[!tfbs$transcript_id %in% interestingGenes]
tfbs$FocusGenes <- tfbs$gene_id %in% interestingGenes$gene_id

tfbs.testlist <- lapply(grep("MA", colnames(tfbs), value = T), function(i) {
    #i = "MA0001.2"
    print(i)
    test <- tfbs %>% select(i, FocusGenes) %>% rename(motif = 1) %>% data.frame()
    
    focus_Hasmotif <- sum(test$motif[which(test$FocusGenes == T)] > 0)
    focus_Notmotif <- sum(test$motif[which(test$FocusGenes == T)] == 0)
    other_Hasmotif <-  sum(test$motif[which(test$FocusGenes == F)] > 0)
    other_Notmotif <- sum(test$motif[which(test$FocusGenes == F)] == 0)
    
    test.matr <- cbind(c(focus_Hasmotif, other_Hasmotif), 
                        c(focus_Notmotif , other_Notmotif))
    
    test.res <- fisher.test(test.matr, alternative = "greater")
    list(test.matr, test.res)
  }
)

names(tfbs.testlist) <- grep("MA", colnames(tfbs), value = T)
length(tfbs.testlist)



pvalues <- sapply(tfbs.testlist, function(i){
              i[[2]][[1]]
})

OR <- sapply(tfbs.testlist, function(i){
  i[[2]][[3]]
})

summary(p.adjust(pvalues, method = "fdr"))
tfbs.testlist[which.min(pvalues)] # most significant motif

table(pvalues<0.01)
names(pvalues<0.01)


#### make final table:

final.pvalues <- data.frame(motif = gsub("_.*", "", names(pvalues)),
           motif.name = gsub(".*_", "", names(pvalues)),
           pvalue = pvalues,
           fdr.pvalue=p.adjust(pvalues, method = "fdr"),
           odds.ratio = OR) %>% as_tibble()


save(final.pvalues, file="TFBS_pvalues_all.RData")

all <- final.pvalues %>% arrange(fdr.pvalue, desc(odds.ratio))
write.csv(all, "output-motifs.csv")

