library(dplyr)
library(pairwiseAdonis)
library(ggplot2)
library(yarrr)

piska_df <- read.table("clipboard",header=TRUE,sep="\t",na.strings="NA",dec=".",strip.white=TRUE) #piska database full

#split datasets into env and species matrices separate for each year

piska_df2 <- piska_df
colnr <- ncol(piska_df2)
piska_df_2007 <- piska_df2 %>% filter(year==2007)
rownames(piska_df_2007) <- piska_df_2007[,1]
piska_df_2011 <- piska_df2 %>% filter(year==2011)
rownames(piska_df_2011) <- piska_df_2011[,1]
piska_df_2013 <- piska_df2 %>% filter(year==2013)
rownames(piska_df_2013) <- piska_df_2013[,1]
piska_df_2015 <- piska_df2 %>% filter(year==2015)
rownames(piska_df_2015) <- piska_df_2015[,1]
piska_df_2019 <- piska_df2 %>% filter(year==2019)
rownames(piska_df_2019) <- piska_df_2019[,1]
piska_df_2020 <- piska_df2 %>% filter(year==2020)
rownames(piska_df_2020) <- piska_df_2020[,1]
piska_df_2021 <- piska_df2 %>% filter(year==2021)
rownames(piska_df_2021) <- piska_df_2021[,1]
piska_df_env_2007 <- piska_df_2007[,1:5]
piska_df_species_2007 <- piska_df_2007[,8:colnr]
piska_df_env_2011 <- piska_df_2011[,1:5]
piska_df_species_2011 <- piska_df_2011[,8:colnr]
piska_df_env_2013 <- piska_df_2013[,1:5]
piska_df_species_2013 <- piska_df_2013[,8:colnr]
piska_df_env_2015 <- piska_df_2015[,1:5]
piska_df_species_2015 <- piska_df_2015[,8:colnr]
piska_df_env_2019 <- piska_df_2019[,1:5]
piska_df_species_2019 <- piska_df_2019[,8:colnr]
piska_df_env_2020 <- piska_df_2020[,1:5]
piska_df_species_2020 <- piska_df_2020[,8:colnr]
piska_df_env_2021 <- piska_df_2021[,1:5]
piska_df_species_2021 <- piska_df_2021[,8:colnr]


#PERMANOVA - done on the full dataset

pairwise_adonizing <- function(year) {
  var1 <- paste("piska_df_species_", year, sep="")
  var1 <- get(var1)
  var2 <- paste("piska_df_env_", year, sep="")
  var2 <- get(var2)
  print(year)
  adonis <- pairwise.adonis(var1, var2$habitat,
                            sim.function = "vegdist",
                            sim.method = "bray",
                            p.adjust.m = "bonferroni",
                            perm = 9999)
  print(adonis)
}
fullpic <- data.frame(pairs=character(),
                      stringsAsFactors=FALSE)
for (i in c(2007,2011,2013,2015,2019,2020,2021)) {
  set.seed(1)
  adonis <- pairwise_adonizing(i)
  pairs <- adonis$pairs
  F.val <- adonis$F
  P <- adonis$p.adjusted
  year <- i
  ready_adonis <- cbind(pairs, F.val,P, year)
  fullpic <- rbind(fullpic, ready_adonis)
}

fullpic$year_since_dist = as.numeric(fullpic$year) - 2002

pdf("permanova_plot.pdf",width=8,height=6)
fullpic$F.val <- as.numeric(fullpic$F.val)
plot <- ggplot(fullpic, aes(x= as.numeric(year_since_dist), y=F.val)) + 
  geom_point(size=3) +
  geom_smooth(method= "lm") +
  theme(panel.background = element_blank(),
        axis.line = element_line(linewidth = 0.5)) +
  scale_x_continuous(breaks=c(8,13,18)) +
  scale_y_continuous(breaks=c(0,5,10,15,20))
plot + facet_wrap(~ pairs)
dev.off()


#permanova analysis on the balanced subsets of data (n=20 subsets, each consisted of n=29 random points per habitat)


fullpic <- data.frame(pairs=character(),
                      stringsAsFactors=FALSE)
for (i in c(2007,2011,2013,2015,2019,2020,2021)) {
    for (j in c(1:20)) {
      
      
  set.seed(j)
  year <- i
  var1 <- paste("piska_df_", year, sep="")
  var12 <- get(var1)
      balanced_design <- var12 %>% 
        group_by(habitat, year) %>%
        sample_n(29, replace = FALSE) %>%
        ungroup() %>%
        as.data.frame()

      env <- balanced_design[,1:11]
      species <- balanced_design[,14:colnr]
      
      adonis <- pairwise.adonis(species, env$habitat,
                                sim.function = "vegdist",
                                sim.method = "bray",
                                p.adjust.m = "bonferroni",
                                perm = 9999)
      
  pairs <- adonis$pairs
  F <- adonis$F.Model
  R2 <- adonis$R2
  P <- adonis$p.adjusted
  ready_adonis <- cbind(pairs, F, R2,P, year)
  fullpic <- rbind(fullpic, ready_adonis)
}
}

permanova_bal <- fullpic
permanova_bal$year_since_dist <- as.numeric(permanova_bal$year) - 2002
permanova_bal$year_since_dist <- as.factor(permanova_bal$year_since_dist)
permanova_bal$P <- as.numeric(permanova_bal$P)
permanova_bal$F <- as.numeric(permanova_bal$F)
permanova_bal$pairs <- as.factor(permanova_bal$pairs)
pdf(file = "permanova-bal_pirateplots.pdf", 
    width = 12,
    height = 8)

cbpal <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")


pirateplot(formula = F ~ year_since_dist + pairs,
           data = permanova_bal,
           theme = 3,
           pal = cbpal,
           main = "effect size (F-stat) distribution (n=20)")

pirateplot(formula = P ~ year_since_dist + pairs,
           data = permanova_bal,
           theme = 3,
           pal = cbpal,
           main = "p-values distribution (n=20)") +
  abline(h=0.05) + abline(h=0.1) + abline(h=0.15) + abline(h=0.2) + abline(h=0.25) + abline(h=0.3)



dev.off()

