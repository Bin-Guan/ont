args <- commandArgs(trailingOnly=TRUE)
#args <- c( "S45", "Z:/projects/bcmORF15/mosdepth/S45.md.regions.bed.gz", "Z:/genome/nisc23-1/bcmlocus/mosdepth/D1695_02.md.regions.bed.gz", "Z:/genome/nisc23-1/bcmlocus/D1695_02.bcm.cn.tsv", "Z:/genome/nisc23-1/bcmlocus/D1695_02.bcm.cn.wide.tsv")

sampleName <- args[1]
coverage_file <- args[2]
coverage_wide_file <- args[3]
# CN calculation; rearrangedGemini_file <- args[3]

library(tidyverse)
library(readxl)
library(vroom)

coverage <- vroom(coverage_file, col_names = FALSE) %>% select(X4, X5) %>% 
  pivot_wider(names_from = X4, values_from = X5) %>% 
  mutate(Sample = sampleName) %>% 
  select(Sample, everything()) %>% 
  mutate(OPN1=pmax(OPN1LW, OPN1MW))

write_tsv(coverage, file = coverage_wide_file)



