library(tidyverse)


datExpr <- readRDS("manuscript/brain/manuscript70/results/wgcna/datExpr_mPFC_CSUB.RDS") 
net <- readRDS("manuscript/brain/manuscript70/results/wgcna/net_mPFC_CSUB_Power6.RDS")
MEs <- readRDS("manuscript/brain/manuscript70/results/wgcna/mPFC_CSUB_MEs.RDS")


ME_df <-MEs%>% data.frame() %>% 
  tibble::rownames_to_column(var = "SampleName") %>%
  dplyr::select(SampleName, MEblack, MEcyan, MEmagenta, MEsienna3, MEblue, MEsteelblue,MEskyblue, MElightgreen, MEdarkgreen, MEpink) %>% 
  pivot_longer(cols = 2:11, names_to = "Module") %>% 
  full_join(coldata)

ME_df$Module <- gsub("ME", "", ME_df$Module)


ME_df <- ME_df %>%
  mutate(status = condition1) %>% 
  filter(condition1 != "DES") %>%
  filter(condition1 != "DOM") %>%
  filter(condition1 != "CDOM")  %>% 
  filter(condition1 != "ascenders") %>%
  filter(condition1 != "descenders") %>%
  filter(condition1 != "control") %>%
  filter(condition1 != "same") 

ME_df$status <- factor(ME_df$status, levels = c("CSUB", "SUB", "ASC"))


ME_df$Module <- factor(ME_df$Module, levels =c("black", "skyblue", "magenta", "cyan", "lightgreen", "steelblue", "sienna3", "darkgreen","pink", "blue"))
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
  facet_wrap(~Module, scales = "free_y", ncol = 5)+
  labs(x = "Social status",
       y = "Module eigengene",
       title = "Module eigengene across social status in mPFC")+ 
  theme_classic(base_size = 15) +
  theme(legend.position = "none") -> p2

p2


ggsave(filename = "manuscript/brain/manuscript70/results/results_figures/mPFC_CSUB_eigengene_significant.png",
       p2,
       height = 8, width = 18, dpi = 300)





### DEG overlap with kIN 



kin <- readRDS("manuscript/brain/manuscript70/results/kIN/mPFC_CSUB_kIN_dataframe_.RDS")


kin_df <- kin %>%  
  arrange(desc(kWithin)) %>% 
  mutate(kINrank = rank(kWithin))


round(nrow(kin_df)*0.20) -> cutoff

kin_df %>% 
  head(cutoff) %>% 
  .$symbol -> y



my_logFC_threshold = 0.2

#mPFC data

limma_list<- readRDS("manuscript/brain/manuscript70/results/RDS/limma_mPFC_ControlSUB.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= 0.2)) %>%
  map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol)))


x_csub <- limma_list$controlsub$symbol
x_casc <- limma_list$controlasc$symbol
x_sa  <- limma_list$subasc$symbol


x_csub[x_csub %in% y] %>% tibble()-> olap_csub
olap_csub

x_casc[x_casc %in% y] %>% tibble()-> olap_casc
olap_casc

x_sa[x_sa %in% y]%>% tibble() -> olap_sa
olap_sa


s.genes <- rbind(olap_csub, olap_casc, olap_sa) %>% unique(.)

s.genes <- t(s.genes)


kin %>% 
  filter(symbol %in% s.genes) -> key_reg_list

write.csv(key_reg_list,"manuscript/brain/manuscript70/mPFC_CSUB_kIN_DEG_ovelap.csv",row.names = F)
