library(tidyverse)


datExpr <- readRDS("manuscript/brain/manuscript70/results/wgcna/datExpr_mPFC_CDOM.RDS") 
net <- readRDS("manuscript/brain/manuscript70/results/wgcna/net_mPFC_CDOM_Power5.RDS")
MEs <- readRDS("manuscript/brain/manuscript70/results/wgcna/mPFC_CDOM_MEs.RDS")


ME_df <-MEs%>% data.frame() %>% 
  tibble::rownames_to_column(var = "SampleName") %>%
  dplyr::select(SampleName, MEyellowgreen, MEcyan, MEfloralwhite, MEsienna3, MEviolet, MEsteelblue) %>% 
  pivot_longer(cols = 2:7, names_to = "Module") %>% 
  full_join(coldata)

ME_df$Module <- gsub("ME", "", ME_df$Module)


ME_df <- ME_df %>%
  mutate(status = condition1) %>% 
  filter(condition1 != "ASC") %>%
  filter(condition1 != "SUB") %>%
  filter(condition1 != "CSUB")


ME_df$status <- factor(ME_df$status, levels = c("CDOM", "DOM", "DES"))

ME_df$Module <- factor(ME_df$Module, levels =c("yellowgreen", "cyan", "violet", "sienna3", "floralwhite", "steelblue"))
ME_df <- na.omit(ME_df)
table(ME_df$Module)
source("functions/geom_boxjitter.R")

ME_df %>%
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = viridis::viridis(3))+
  scale_fill_manual(values = viridis::viridis(3))+
  facet_wrap(~Module, scales = "free_y", ncol = 7)+
  labs(x = "Social status",
       y = "Module eigengene",
       title = "Module eigengene across social status in mPFC")+ 
  theme_classic(base_size = 15) +
  theme(legend.position = "none") -> p2

p2


ggsave(filename = "manuscript/brain/manuscript70/results/results_figures/mPFC_CDOM_eigengene_significant.png",
       p2,
       height = 8, width = 18, dpi = 300)





### DEG overlap with kIN 



kin <- readRDS("manuscript/brain/manuscript70/results/kIN/mPFC_CDOM_kIN_dataframe_.RDS")

kin_df <- kin %>%  
  arrange(desc(kWithin)) %>% 
  mutate(kINrank = rank(kWithin))


round(nrow(kin_df)*0.20) -> cutoff

kin_df %>% 
  head(cutoff) %>% 
  .$symbol -> y



my_logFC_threshold = 0.2

#mPFC data

limma_list<- readRDS("manuscript/brain/manuscript70/results/RDS/limma_mPFC_ControlDD.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= 0.2)) %>%
  map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol)))


x_cdom <- limma_list$controldom$symbol
x_cdes <- limma_list$controldes$symbol
x_dd  <- limma_list$domdes$symbol


x_cdom[x_cdom %in% y] %>% as.tibble()-> olap_cdom
olap_cdom
  
x_cdes[x_cdes %in% y] %>% as.tibble()-> olap_cdes
olap_cdes

x_dd[x_dd %in% y]%>% as.tibble() -> olap_dd
olap_dd


s.genes <- rbind(olap_cdom, olap_cdes, olap_dd) %>% unique(.)

s.genes <- t(s.genes)


kin %>% 
    filter(symbol %in% s.genes) -> key_reg_list
  


 
write.csv(key_reg_list,"manuscript/brain/manuscript70/results/tables/mPFC_CDOM_kIN_DEG_ovelap.csv",row.names = F)

