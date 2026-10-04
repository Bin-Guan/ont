args <- commandArgs(trailingOnly=TRUE)

# args <- c("Z:/exome/blueprint/AutoMap/G4V9_1_BP93806/G4V9_1_BP93806.HomRegions.tsv",
#           "Z:/exome/blueprint/AutoMap/test.annotSV.tsv", "Z:/resources/OGLpanelGeneDxORcandidate.xlsx", "Z:/exome/blueprint/AutoMap/test.annotSV.output.tsv")
sample <- args[1]
annovar_file <- args[2]
output_file <- args[3]
#geneCategory_file <- args[3]
#annotated_file <- args[4]

library(tidyverse)
library(readxl)

annovar <- read_tsv(annovar_file, col_names = TRUE, na = c("NA", "", "None", "."), col_types = cols(.default = col_character())) %>%
  mutate(INFO = case_when( INFO == "P" ~ "pileup",
                           INFO == "F" ~ "full-alignment",
                           TRUE ~ INFO)) %>%
  type.convert() %>% 
  mutate(Sample=sample) %>% 
  mutate(Note = ifelse(POS %in% c(38285414, 38298269, 38299739), "Homopolymer", "")) %>% 
  unite("refgenewithver", GeneDetail.refGeneWithVer, AAChange.refGeneWithVer, sep = ",", remove = TRUE, na.rm = TRUE) %>% 
  select(Sample, CHROM:GT_FIELDS, Note, `Gene.refGeneWithVer`, Func.refGeneWithVer, ExonicFunc.refGeneWithVer, refgenewithver )

plof <- filter(annovar, grepl("splicing", Func.refGeneWithVer) | grepl("^frameshift|stop|start", ExonicFunc.refGeneWithVer) )

openxlsx::write.xlsx(list("pLoF" = plof, "orf15" = annovar), file = output_file, firstRow = TRUE, firstCol = FALSE)
