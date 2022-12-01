library(tidyverse)


datExpr <- readRDS("manuscript/brain/manuscript70/results/wgcna/datExpr_MEA_CDOM.RDS") 
net <- readRDS("manuscript/brain/manuscript70/results/wgcna/net_MEA_CDOM_Power9.RDS")
MEs <- readRDS("manuscript/brain/manuscript70/results/wgcna/MEA_CDOM_MEs.RDS")


ME_df <-MEs%>% data.frame() %>% 
  tibble::rownames_to_column(var = "SampleName") %>%
  pivot_longer(cols = 2:23, names_to = "Module") %>% 
  full_join(coldata)

ME_df$Module <- gsub("ME", "", ME_df$Module)

ME_df <- ME_df %>%
  mutate(status = condition1) %>% 
  filter(condition1 != "ASC") %>%
  filter(condition1 != "SUB") %>%
  filter(condition1 != "CSUB") %>% filter(Module != "grey") %>% 
  filter(Module != "grey60") %>% filter(Module != "	magenta") %>% 
  filter(Module != "turquoise") %>% filter(Module != "lightgreen") %>% 
  filter(Module != "lightyellow") %>% filter(Module != "blue") %>% 
  filter(Module != "cyan")


ME_df$status <- factor(ME_df$status, levels = c("CDOM", "DOM", "DES"))

ME_df$Module <- factor(ME_df$Module, levels =c("darkred", 'black', "brown", "green","tan", "purple",
                                               "greenyellow","midnightblue", "lightcyan","red", 
                                                "pink","royalblue",
                                               'yellow', "salmon"))
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
       title = "Module eigengene across social status in MEA")+ 
  theme_classic(base_size = 15) +
  theme(legend.position = "none") -> p2

p2


ggsave(filename = "manuscript/brain/manuscript70/results/results_figures/MEA_CDOM_eigengene_significant.png",
       p2,
       height = 8, width = 18, dpi = 300)





### DEG overlap with kIN 



kin <- readRDS("manuscript/brain/manuscript70/results/kIN/MEA_CDOM_kIN_dataframe_.RDS")

kin_df <- kin %>%  
  arrange(desc(kWithin)) %>% 
  mutate(kINrank = rank(kWithin))


round(nrow(kin_df)*0.20) -> cutoff

kin_df %>% 
  head(cutoff) %>% 
  .$symbol -> y



my_logFC_threshold = 0.2

#MEA data

limma_list<- readRDS("manuscript/brain/manuscript70/results/RDS/limma_MEA_ControlDD.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= 0.2)) %>%
  map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol))) %>% 
  map(~dplyr::select(.,engene)) 


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
  


 
write.csv(key_reg_list,"manuscript/brain/manuscript70/kIN_DEG_ovelap.csv",row.names = F)

