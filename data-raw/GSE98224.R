## GSE98224 ---------------------------------------------------------------------------
library(GEOquery)
library(dplyr)
library(readr)
library(readxl)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(janitor)

# processed data collected by Amy Inkster
GSE98224_dnam <- readRDS('Z:/Amy/Project_PE/Data/Cox 48 Sample Subset with DNAme/cox48_450K.RDS')
GSE98224_expr <- readRDS('Z:/Amy/Project_PE/Data/Cox 48 Sample Subset with DNAme/cox48_expression.RDS')
GSE98224_pDat <- readRDS('Z:/Amy/Project_PE/Data/Cox 48 Sample Subset with DNAme/cox48_pData.RDS')

# sample key to link sample IDs between datasets
samplekey <- read_xlsx('Z:/Victor/Projects/DNAm-Ethnicity-Predictor/data/Cox/2019-03-13 From Katie - Sample to GEO ID.xlsx') %>% clean_names


                       # clean pdata
GSE98224_pDat <- GSE98224_pDat %>% as_tibble %>% 
  select(geo_accession, diagnosis, tissue, maternal_age, maternal_bmi, maternal_ethnicity, 
         `ga_(week)`, `ga_(day)`) %>%
  dplyr::rename(ga_weeks = `ga_(week)`, ga_days = `ga_(day)`) %>%
  mutate_all(function(x){gsub('\\s', '', x)}) %>%
  mutate(maternal_age = as.numeric(maternal_age),
         maternal_bmi = as.numeric(maternal_bmi),
         ga_weeks = as.numeric(ga_weeks),
         ga_days = as.numeric(ga_days)) 

# add methylation geo accessions
GSE98224_pDat$geo_accession %in% samplekey$ge_geo_id %>% all() #T
GSE98224_pDat <- GSE98224_pDat %>% left_join(samplekey[,3:4], by = c('geo_accession' = 'ge_geo_id'))

# combine with some DNA methylation and expression data
# use PE-associated hits to subset dna methylation data
hits <- read_table2('Z:/Victor/Projects/DNAm-Ethnicity-Predictor/Results/08_Cox_Roblab_651overlappingPEhits.txt')

rownames(GSE98224_dnam) <- GSE98224_dnam[,'ID_REF']
GSE98224_dnam <- t(GSE98224_dnam[intersect(GSE98224_dnam[,'ID_REF'], hits$cpg),-1]) %>% as.data.frame() %>%
  mutate(geo_accession = rownames(.)) %>%
  as_tibble()
set.seed(1)
GSE98224_dnam <- GSE98224_dnam[,652] %>% bind_cols(GSE98224_dnam[,sample(1:651, 100)])

# convert to numeric
GSE98224_dnam <- bind_cols(GSE98224_dnam %>% select(geo_accession), 
          as_tibble(sapply(GSE98224_dnam[,-1], function(x) (as.numeric(as.character(x)))))) %>%
  as_tibble()

# sample 50 random transcripts, and transpose
set.seed(1)
GSE98224_expr <- GSE98224_expr %>% sample_n(50) 
rn <- GSE98224_expr$ID_REF
GSE98224_expr <- GSE98224_expr[,-1] %>% t() %>% as.data.frame() %>% mutate(geo_accession = rownames(.)) %>%
  as_tibble() %>% select(geo_accession, everything())
colnames(GSE98224_expr)[2:51] <- paste0('transcript_', rn)

# join dnam, expr, and pdat
GSE98224_pDat <- GSE98224_pDat %>% left_join(GSE98224_expr) 
GSE98224_pDat <- GSE98224_pDat %>% left_join(GSE98224_dnam, by = c('meth_geo_id' = 'geo_accession')) %>%
  dplyr::rename(expr_geo_id = geo_accession) %>%
  select(expr_geo_id, meth_geo_id, everything()) 

write.csv(GSE98224_pDat,'../data/GSE98224.csv')

#---------------------------------------------------------------------------------------------------
#get annotation, not used yet
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno_dnam <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno_dnam <- anno_dnam %>% # subset to data
  as_tibble() %>% 
  filter(Name %in% GSE98224_dnam[,'ID_REF'])

anno_expr <- readRDS('Z:/Amy/Project_PE/Data/Cox 48 Sample Subset with DNAme/Affymetrix10HumanSTArray_Anno.RDS') 
anno_expr <- anno_expr %>% as_tibble %>%
  filter(ID %in% GSE98224_expr$ID_REF)

