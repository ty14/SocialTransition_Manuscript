library(tidyverse)

#transition rank genes 

des <- read_csv('manuscript/brain/manuscript70/results/tables/descenders_tranisiton_MEA_genes.csv')

asc <- read_csv("manuscript/brain/manuscript70/results/tables/ascenders_tranisiton_MEA_genes.csv")

asc$a_reg <- asc$reg

asc <- asc %>% select(-reg)


trans <- des %>% full_join(asc)

na.omit(trans)

trans_genes <- trans %>% filter(reg == a_reg)

#####reorganized genes

csub <- read_csv('manuscript/brain/manuscript70/results/tables/csub_reorganization_MEA_genes.csv')

cdom <- read_csv('manuscript/brain/manuscript70/results/tables/cdom_reorganization_MEA_genes.csv')

csub$s_reg <- csub$reg
csub <- csub %>%  select(-reg)

reorg <- cdom %>% full_join(csub)

na.omit(reorg)

reorg_genes <- reorg %>% filter(reg == s_reg)


## stayed same rank genes
cdom <- cdom %>% select(symbol, reg)
csub$s_reg <- csub$reg
csub <- csub %>%  select(symbol, s_reg)

stay <- cdom %>% full_join(csub)
stay
stay_genes <- stay %>% filter(reg == s_reg) %>% arrange(reg) %>% unique(.)


## social defeat genes
## looking for genes different from SUB: 
cdd <- sa %>% filter(reg == "up")

cds <- csub %>% filter(reg == "down")

up_csub <- cdd$symbol[(cdd$symbol %in% cds$symbol)] %>% as.data.frame(.)
colnames(up_csub)<-"symbol"

up_csub <- up_csub %>% mutate(reg = "Up")

cadd <- sa %>% filter(reg == "down")

cdsub<- csub %>% filter(reg == "up")

down_csub <- cadd$symbol[(cadd$symbol %in% cdsub$symbol)] %>% as.data.frame(.)
colnames(down_csub)<-"symbol"

down_csub <- down_csub %>% mutate(reg = "Down")

#up csub is 166, down csub is 124
sub_gene <- up_csub %>% full_join(down_csub)

write.csv(sub_gene, 'manuscript/brain/manuscript70/results/tables/SUB_defeatstress_MEA_genes.csv' ,row.names = F)

sub_gene$s_reg <- sub_gene$reg 
sub <- sub_gene %>% select (-reg)


defeat <- des %>% full_join(sub)
defeat
defeat_genes <- defeat %>% filter(reg == s_reg) %>% arrange(reg) %>% unique(.)

#########

cdesx <- cdes %>% select(symbol, reg)
csubx  <- csub %>%  select (symbol, s_reg = reg)

defeatx <- cdesx %>% full_join(csubx)
defeatx
defeat_genes <- defeatx %>% filter(reg == s_reg) %>% arrange(reg) %>% unique(.)

table(defeat_genes$reg)


######
#aggressive genes?
## looking for genes different from DOM: 
cdd <- domdes %>% filter(reg == "up")

cds <- cdom %>% filter(reg == "down")

up_csub <- cdd$symbol[(cdd$symbol %in% cds$symbol)] %>% as.data.frame(.)
colnames(up_csub)<-"symbol"

up_csub <- up_csub %>% mutate(reg = "Up")

cadd <- domdes %>% filter(reg == "down")

cdsub<- cdom %>% filter(reg == "up")

down_csub <- cadd$symbol[(cadd$symbol %in% cdsub$symbol)] %>% as.data.frame(.)
colnames(down_csub)<-"symbol"

down_csub <- down_csub %>% mutate(reg = "Down")

#up csub is 166, down csub is 124
dom_gene <- up_csub %>% full_join(down_csub)

write.csv(dom_gene, 'manuscript/brain/manuscript70/results/tables/DOM_aggression_MEA_genes.csv' ,row.names = F)

dom_gene$d_reg <- dom_gene$reg 
dom <- dom_gene %>% select (-reg)

agg <- asc %>% full_join(dom)
dom
agg_genes <- agg %>% filter(reg == d_reg) %>% arrange(reg) %>% unique(.)


cdomx <- cdom %>% select(symbol, reg)
cascx <- casc %>% select(symbol, s_reg = reg)

aggx <- cdomx %>% full_join(cascx)
agg_genesx <- aggx %>% filter(reg == s_reg) %>% arrange(reg) %>% unique(.)



