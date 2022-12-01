library(limma)
library(Mus.musculus)
library(tidyverse)
library(scales)


library(tidyverse)
# library(edgeR)
library(DESeq2)
library(pheatmap)
library(annotables)
grcm38 # mouse genes



rawcount_list <- list()
temp = list.files(path = "manuscript/brain",pattern="*_counts.csv")
temp0 <- paste0("manuscript/brain/",temp)
rawcount_list = lapply(temp0, read_csv)
regions <- gsub("_counts.csv","",temp)

# behavior data cleaned from rna_seq github repo 
names(rawcount_list) <- regions
lapply(rawcount_list, head)

rawcount_list<- rawcount_list %>% map(~mutate(., ensgene = X)) %>% map(~dplyr::select(., -X))

counts <- rawcount_list %>%
  reduce(full_join, by = 'ensgene')

# total gene counts per sample ========================================

counts %>% 
  summarize_if(is.numeric,sum,na.rm = T) %>% 
  t() %>% 
  as.data.frame %>% 
  rename(genecounts = V1) %>% 
  rownames_to_column(var = 'sampleID') -> genecounts


genecounts %>% 
  ggplot(aes(genecounts)) +
  geom_histogram(bins = 40,alpha =0.5,color = 'grey') +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_y_continuous(breaks = pretty_breaks()) +
  theme_minimal() -> p_genecounts

print(p_genecounts)
# ggsave(p_genecounts,filename = 'results_figures/p_genecounts.png',width = 8, height = 5)
genecounts %>% 
  filter(genecounts < 800000)                    


genecounts %>% 
  filter(genecounts > 3000000) 



# create behav dataframe
#Getting metadata ready 
coldata <- read_csv("manuscript/brain/sample70min_table.csv")
head(coldata)
str(coldata)

# Normalizing cort data
# df <- transform(df, N = (N - min(N)) / (max(N) - min(N))

coldata <-coldata %>% mutate(post_Ncort =  (mean_con_ng_ul - min(mean_con_ng_ul,na.rm=T))/(max(mean_con_ng_ul,na.rm=T)-min(mean_con_ng_ul,na.rm=T)))
coldata <- coldata %>%  dplyr::select(-group, -period)
coldata$post_Ncort[is.na(coldata$post_Ncort)] = 0

##getting condition1
# CDOM, RDOM to Descenders (DOM to SUB)(4->1)
# CSUB, SUB to Ascenders (Sub to DOM)  (1->4)
table(coldata$condition)

coldata$condition1 <- ifelse(coldata$condition == "same" & coldata$Prerank == 1, "DOM", coldata$condition)
coldata$condition1 <- ifelse(coldata$condition1 == "same" & coldata$Prerank == 4, "SUB", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "descenders" & coldata$Postrank == 4 & coldata$Prerank == 1, "DES", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "ascenders" & coldata$Prerank == 4 & coldata$Postrank == 1, "ASC", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "control" & coldata$Postrank == 4, "CSUB", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "control" & coldata$Postrank == 1, "CDOM", coldata$condition1)

coldata <- coldata %>% 
  filter(Postrank != 3) %>% 
  filter(condition1 != 'ascenders')  %>%
  dplyr::select(-mean_con_ng_ul, -condition)

write_csv(coldata, "manuscript/brain/manuscript70/results/tables/coldata_70min.csv")

table(coldata$condition1)

# First AMY data 
amy_data<- coldata %>% filter(region == "AMY")%>% arrange(condition1)

row.names <- amy_data$SampleName

amy_data1 <- amy_data%>% dplyr::select(-SampleName)
row.names(amy_data1) <- row.names #Assigning row names from as sample names  
head(amy_data1)


# Bring in count data for mAMY
a_countdata <- read_csv("manuscript/brain/AMY_counts.csv")
a_countdata$X -> nrows
a_count <- a_countdata[,-1]
row.names(a_count) <- nrows
a_count <-as.data.frame(a_count)

#filter
a_count <- a_count[rowSums(a_count >= 10) > round((length(amy_data1))*0.9), ]

#checks
all (row.names(amy_data1) %in% colnames(a_count)) #check 

a_count <- a_count[,rownames(amy_data1)]
all(rownames(amy_data1) == colnames(a_count)) #check



# Second mpfc data data 
pfc_data<- coldata %>% filter(region != "AMY")%>% arrange(condition1) 

row.names <- pfc_data$SampleName

pfc_data1 <- pfc_data%>% dplyr::select(-SampleName)

row.names(pfc_data1) <- row.names #Assigning row names from as sample names  
head(pfc_data1)

# Bring in count data for mpfc
p_countdata <- read_csv("manuscript/brain/PFC_counts.csv")
p_countdata$X -> nrows
p_count <- p_countdata[,-1]
row.names(p_count) <- nrows
p_count <-as.data.frame(p_count)

#filter
p_count <- p_count[rowSums(p_count >= 10) > round((length(pfc_data1))*0.9), ]
# pkeep <- rowSums(p_count >= 10)
# p_count <- p_count[pkeep,]

all (row.names(pfc_data1) %in% colnames(p_count)) #check 

p_count <- p_count[, rownames(pfc_data1)]
all(rownames(pfc_data1) == colnames(p_count)) #check



#Just for DOMS 

#MeA
dag <- DESeqDataSetFromMatrix(countData = a_count,
                              colData = amy_data1,
                              design = ~condition1)



dag_sub <- dag[, dag$condition1 %in% c("CDOM", "DOM", "DES")]
dag_sub$condition1 <- droplevels(dag_sub$condition1)
dag_sub$condition1
dag<- dag_sub
dag$condition1


nrow(dag) # 14811 lost about half

#mPFC
dpg <- DESeqDataSetFromMatrix(countData = p_count,
                              colData = pfc_data1,
                              design = ~condition1)


dpg_sub <- dpg[, dpg$condition1 %in% c("CDOM", "DOM", "DES")]
dpg_sub$condition1 <- droplevels(dpg_sub$condition1)
dpg_sub$condition1
dpg<- dpg_sub
dpg$condition1


nrow(dpg) # 14875

#dds
da_dom <- DESeq(dag)
dp_dom <- DESeq(dpg)

DESeq_domlist<- list(MeA=da_dom, mPFC=dp_dom)

# adjusted pval and fold change threshold to 'define' DEG
adj_pval_cutoff = 0.1
threshold = 0.9  # 90% of the subject
LFC_threshold = 0.2
pvalue_threshold = 0.05 # didn't use it for now



get_LFC_results <- function(dds)
{ 
  res1 <- results(dds, contrast=c("condition1","CDOM","DOM")) %>% 
    as.data.frame() %>% 
    rownames_to_column('ensgene') %>% 
    as_tibble() %>% 
    dplyr::select(ensgene, log2FoldChange, pvalue, padj) %>% 
    mutate(contrast = "CDOM - DOM") %>% 
    filter(abs(log2FoldChange)>=LFC_threshold) %>%
    filter(padj <=adj_pval_cutoff ) %>%
    arrange(desc(abs(log2FoldChange))) %>% 
    left_join(x=., y =grcm38[,c("ensgene", "symbol", "description")],by = "ensgene") %>% 
    unique()
  
  res2 <- results(dds, contrast=c("condition1","DOM","DES")) %>% 
    as.data.frame() %>% 
    rownames_to_column('ensgene') %>% 
    as_tibble() %>% 
    dplyr::select(ensgene, log2FoldChange, pvalue, padj) %>% 
    mutate(contrast = "DOM - DES") %>% 
    filter(abs(log2FoldChange)>=LFC_threshold) %>%
    filter(padj <=adj_pval_cutoff ) %>% filter(pvalue <=pvalue_threshold) %>%
    arrange(desc(abs(log2FoldChange))) %>% 
    left_join(., y =grcm38[,c("ensgene", "symbol", "description")],by = "ensgene") %>% 
    unique()
  
  
  
  res3 <- results(dds, contrast=c("condition1","CDOM","DOM")) %>% 
    as.data.frame() %>% 
    rownames_to_column('ensgene') %>% 
    as_tibble() %>% 
    dplyr::select(ensgene, log2FoldChange, pvalue, padj) %>% 
    mutate(contrast = "CDOM - DOM") %>% 
    filter(abs(log2FoldChange)>=LFC_threshold) %>%
    filter(padj <=adj_pval_cutoff ) %>% filter(pvalue <=pvalue_threshold) %>%
    arrange(desc(abs(log2FoldChange))) %>% 
    left_join(., y =grcm38[,c("ensgene", "symbol", "description")],by = "ensgene") %>% 
    unique()
  
  res4 <- results(dds, contrast=c("condition1","CDOM","DES")) %>% 
    as.data.frame() %>% 
    rownames_to_column('ensgene') %>% 
    as_tibble() %>% 
    dplyr::select(ensgene, log2FoldChange, pvalue, padj) %>% 
    mutate(contrast = "CDOM - DES") %>% 
    filter(abs(log2FoldChange)>=LFC_threshold) %>%
    filter(padj <=adj_pval_cutoff ) %>% filter(pvalue <=pvalue_threshold) %>%
    arrange(desc(abs(log2FoldChange))) %>% 
    left_join(., y =grcm38[,c("ensgene", "symbol", "description")],by = "ensgene") %>% 
    unique()
  
  res5 <- results(dds, contrast=c("condition1","CDOM","DES")) %>% 
    as.data.frame() %>% 
    rownames_to_column('ensgene') %>% 
    as_tibble() %>% 
    dplyr::select(ensgene, log2FoldChange, pvalue, padj) %>% 
    mutate(contrast = "CDOM - DES") %>% 
    filter(abs(log2FoldChange)>=LFC_threshold) %>%
    filter(padj <=adj_pval_cutoff ) %>% filter(pvalue <=pvalue_threshold) %>%
    arrange(desc(abs(log2FoldChange))) %>% 
    left_join(., y =grcm38[,c("ensgene", "symbol", "description")],by = "ensgene") %>% 
    unique()
  
  res6 <- results(dds, contrast=c("condition1","DOM","DES")) %>% 
    as.data.frame() %>% 
    rownames_to_column('ensgene') %>% 
    as_tibble() %>% 
    dplyr::select(ensgene, log2FoldChange, pvalue, padj) %>% 
    mutate(contrast = "DOM - DES") %>% 
    filter(abs(log2FoldChange)>=LFC_threshold) %>%
    filter(padj <=adj_pval_cutoff ) %>% filter(pvalue <=pvalue_threshold) %>%
    arrange(desc(abs(log2FoldChange))) %>% 
    left_join(., y =grcm38[,c("ensgene", "symbol", "description")],by = "ensgene") %>% 
    unique()
  
  
  rbind(res1, res2, res3, res4,res5, res6) %>% 
    arrange(desc(abs(log2FoldChange))) -> result
  
  return(result)
}

results_domlist <- DESeq_domlist %>% map(get_LFC_results)
results_dom <- results_domlist %>% map2_df(.,names(.), ~mutate(.x,region = .y))%>%  arrange(desc(log2FoldChange)) %>% unique(.)
head(results_dom)
tail(results_dom)

table(results_dom$contrast, results_dom$region)

# write.csv(results_dom, 'manuscript/brain/manuscript70/results/ControlDom_DEG_table.csv' ,row.names = F)



#Just for Subs 

#MeA
dag <- DESeqDataSetFromMatrix(countData = a_count,
                              colData = amy_data1,
                              design = ~condition1)



dag_sub <- dag[, dag$condition1 %in% c("CSUB", "SUB", "ASC")]
dag_sub$condition1 <- droplevels(dag_sub$condition1)
dag_sub$condition1
dag<- dag_sub
dag$condition1


nrow(dag) # 14811 lost about half

#mPFC
dpg <- DESeqDataSetFromMatrix(countData = p_count,
                              colData = pfc_data1,
                              design = ~condition1)


dpg_sub <- dpg[, dpg$condition1 %in% c("CSUB", "SUB", "ASC")]
dpg_sub$condition1 <- droplevels(dpg_sub$condition1)
dpg_sub$condition1
dpg<- dpg_sub
dpg$condition1


nrow(dpg) # 14875

#dds
da_sub <- DESeq(dag)
dp_sub <- DESeq(dpg)

DESeq_sublist<- list(MeA=da_sub, mPFC=dp_sub)

# adjusted pval and fold change threshold to 'define' DEG
adj_pval_cutoff = 0.1
threshold = 0.9  # 90% of the subject
LFC_threshold = 0.2
pvalue_threshold = 0.05 # didn't use it for now



get_LFC_results <- function(dds)
{ 
  res1 <- results(dds, contrast=c("condition1","CSUB","SUB")) %>% 
    as.data.frame() %>% 
    rownames_to_column('ensgene') %>% 
    as_tibble() %>% 
    dplyr::select(ensgene, log2FoldChange, pvalue, padj) %>% 
    mutate(contrast = "CSUB- SUB") %>% 
    filter(abs(log2FoldChange)>=LFC_threshold) %>%
    filter(padj <=adj_pval_cutoff ) %>%
    arrange(desc(abs(log2FoldChange))) %>% 
    left_join(x=., y =grcm38[,c("ensgene", "symbol", "description")],by = "ensgene") %>% 
    unique()
  
  res2 <- results(dds, contrast=c("condition1","CSUB","ASC")) %>% 
    as.data.frame() %>% 
    rownames_to_column('ensgene') %>% 
    as_tibble() %>% 
    dplyr::select(ensgene, log2FoldChange, pvalue, padj) %>% 
    mutate(contrast = "CSUB - ASC") %>% 
    filter(abs(log2FoldChange)>=LFC_threshold) %>%
    filter(padj <=adj_pval_cutoff ) %>% filter(pvalue <=pvalue_threshold) %>%
    arrange(desc(abs(log2FoldChange))) %>% 
    left_join(., y =grcm38[,c("ensgene", "symbol", "description")],by = "ensgene") %>% 
    unique()
  
  
  
  res3 <- results(dds, contrast=c("condition1","SUB","ASC")) %>% 
    as.data.frame() %>% 
    rownames_to_column('ensgene') %>% 
    as_tibble() %>% 
    dplyr::select(ensgene, log2FoldChange, pvalue, padj) %>% 
    mutate(contrast = "SUB - ASC") %>% 
    filter(abs(log2FoldChange)>=LFC_threshold) %>%
    filter(padj <=adj_pval_cutoff ) %>% filter(pvalue <=pvalue_threshold) %>%
    arrange(desc(abs(log2FoldChange))) %>% 
    left_join(., y =grcm38[,c("ensgene", "symbol", "description")],by = "ensgene") %>% 
    unique()
  
  res4 <- results(dds, contrast=c("condition1","CSUB","SUB")) %>% 
    as.data.frame() %>% 
    rownames_to_column('ensgene') %>% 
    as_tibble() %>% 
    dplyr::select(ensgene, log2FoldChange, pvalue, padj) %>% 
    mutate(contrast = "CSUB - SUB") %>% 
    filter(abs(log2FoldChange)>=LFC_threshold) %>%
    filter(padj <=adj_pval_cutoff ) %>% filter(pvalue <=pvalue_threshold) %>%
    arrange(desc(abs(log2FoldChange))) %>% 
    left_join(., y =grcm38[,c("ensgene", "symbol", "description")],by = "ensgene") %>% 
    unique()
  
  res5 <- results(dds, contrast=c("condition1","CSUB","ASC")) %>% 
    as.data.frame() %>% 
    rownames_to_column('ensgene') %>% 
    as_tibble() %>% 
    dplyr::select(ensgene, log2FoldChange, pvalue, padj) %>% 
    mutate(contrast = "CSUB - ASC") %>% 
    filter(abs(log2FoldChange)>=LFC_threshold) %>%
    filter(padj <=adj_pval_cutoff ) %>% filter(pvalue <=pvalue_threshold) %>%
    arrange(desc(abs(log2FoldChange))) %>% 
    left_join(., y =grcm38[,c("ensgene", "symbol", "description")],by = "ensgene") %>% 
    unique()
  
  res6 <- results(dds, contrast=c("condition1","SUB","ASC")) %>% 
    as.data.frame() %>% 
    rownames_to_column('ensgene') %>% 
    as_tibble() %>% 
    dplyr::select(ensgene, log2FoldChange, pvalue, padj) %>% 
    mutate(contrast = "SUB - ASC") %>% 
    filter(abs(log2FoldChange)>=LFC_threshold) %>%
    filter(padj <=adj_pval_cutoff ) %>% filter(pvalue <=pvalue_threshold) %>%
    arrange(desc(abs(log2FoldChange))) %>% 
    left_join(., y =grcm38[,c("ensgene", "symbol", "description")],by = "ensgene") %>% 
    unique()
  
  
  rbind(res1, res2, res3, res4,res5, res6) %>% 
    arrange(desc(abs(log2FoldChange))) -> result
  
  return(result)
}

results_sublist <- DESeq_sublist %>% map(get_LFC_results)
results_sub <- results_sublist %>% map2_df(.,names(.), ~mutate(.x,region = .y))%>%  arrange(desc(log2FoldChange)) %>% unique(.)
head(results_sub)
tail(results_sub)

table(results_sub$contrast, results_sub$region)

 # write.csv(results_sub, 'manuscript/brain/manuscript70/results/ControlSUB_DEG_table.csv' ,row.names = F)
