### R command to count the number of tumor and normal samples in DNA methylation
### Same way we can use for other data type :: data_category=='DNA Methylation'
pacman::p_load(DT, GenomicDataCommons, magrittr, dplyr) 

####### How many total samples in GDC 
res = cases() %>% facet("project.project_id") %>% aggregations()
res

##### How many total cases in different classes
files() %>% facet('type') %>% aggregations() ### Total in GDC
### Total in specific cancer with gene expression
files() %>% facet('type') %>% filter( ~ cases.project.project_id == 'TCGA-CHOL' & type == 'gene_expression')%>%aggregations()
#Specific cases with DNA Methylation
files() %>% facet('type') %>% filter( ~ cases.project.project_id == 'TCGA-CHOL' &  type=="DNA Methylation")%>%aggregations()

## For DNA methylation 
q = files() %>%
  GenomicDataCommons::select(available_fields('files')) %>%
  filter(~ cases.project.project_id=='TCGA-CHOL' &
           data_category=='DNA Methylation')
q %>% facet(c('platform', 'cases.samples.sample_type')) %>% aggregations() ## This will print total


q = files() %>%
  GenomicDataCommons::select(available_fields('files')) %>%
  filter(~ cases.project.project_id=='TCGA-CHOL' &
           data_category=='DNA Methylation')
q %>% facet(c('cases.submitter_id')) %>% aggregations()

########### For RNASeq HTSeq - Counts
q = files() %>%
  GenomicDataCommons::select(available_fields('files')) %>%
  filter(~ cases.project.project_id=='TCGA-CHOL' &
           data_type=='Gene Expression Quantification')
q %>% facet(c('analysis.workflow_type', 'cases.samples.sample_type')) %>% aggregations()


q = files() %>%
  GenomicDataCommons::select(available_fields('files')) %>%
  filter(~ cases.project.project_id=='TCGA-CHOL' &
           data_type=='Gene Expression Quantification')
q %>% facet(c('analysis.workflow_type', 'cases.submitter_id')) %>% aggregations() ## With TCGA id

########## for miRNA
q = files() %>%
  GenomicDataCommons::select(available_fields('files')) %>%
  filter(~ cases.project.project_id=='TCGA-CHOL' &
           data_type=='miRNA Expression Quantification')
q %>% facet(c('analysis.workflow_type', 'cases.samples.sample_type')) %>% aggregations()


q = files() %>%
  GenomicDataCommons::select(available_fields('files')) %>%
  filter(~ cases.project.project_id=='TCGA-CHOL' &
           data_type=='miRNA Expression Quantification')
q %>% facet(c('analysis.workflow_type', 'cases.submitter_id')) %>% aggregations() ## With TCGA ID