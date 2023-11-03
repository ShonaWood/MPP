library(tidyverse)

upstream = 2000
downstream = 200
gft_file <- "hamster_gtf_gff/sib-ham-polished-partials-included-fix.gtf"

gff_file <- "hamster_gtf_gff/sib-ham-partials-included-polished.GFF3"


gffTbl <- read_tsv(gff_file, col_names = c("seqname","source","feature","start","end","dot1","strand","dot2","attribute"))  %>%
  # get only feature that has transcript  id
  filter(grepl("^ID=transcript",attribute)) %>%
  # add temporary index
  mutate(idx = 1:n()) %>% 
  # split the attributes
  mutate(attribute = strsplit(attribute,split = ";")) %>% 
  unnest(attribute) %>% 
  separate(attribute,c("key","value"), sep="=") %>% 
  spread(key = key,value = value) %>% 
  mutate(gene_id=gsub("gene:","",Parent)) %>%
  select(transcript_id, gene_id, Name, biotype)

gtfTbl <- read_tsv(gft_file, col_names = c("seqname","source","feature","start","end","dot1","strand","dot2","attributes"))  %>%
  filter(feature == "transcript") %>%
  mutate("transcript_id" = gsub("transcript_id \"transcript:([^\"]+).+","\\1",attributes)) %>%
  mutate("start_promoter" = ifelse(strand == "+", start - upstream, end - downstream)) %>%
  mutate("end_promoter" = ifelse(strand == "+", start + downstream, end + upstream)) %>%
  rename(start_transcript=start,end_transcript=end) %>%
  select(seqname,start_promoter,end_promoter,strand,transcript_id,start_transcript,end_transcript) 

res <- left_join(gtfTbl, gffTbl)

write_tsv(res, paste0("sib-ham-polished-partials-included-fix.gtf_promoters_u",upstream,"_d",downstream,".bed"))

