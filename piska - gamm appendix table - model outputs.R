piska_df <- read.table("clipboard",header=TRUE,sep="\t",na.strings="NA",dec=".",strip.white=TRUE)
View(piska_df)
lesne <- as.vector(read.table("clipboard",header=FALSE,sep="\t",na.strings="NA",dec=".",strip.white=TRUE))
lesne <- lesne$V1
farmlandowe <- as.vector(read.table("clipboard",header=FALSE,sep="\t",na.strings="NA",dec=".",strip.white=TRUE))
farmlandowe <- farmlandowe$V1

gatunki_klucz <- read.table("clipboard",header=TRUE,sep="\t",na.strings="NA",dec=".",strip.white=TRUE)


library(ecotraj)
library(vegan)
library(grid)
library(gridExtra)
library(lme4)
library(dplyr)
library(mgcv)
library(ggeffects)
library(ggplot2)
library(logistf)
library(ncf)
library(cowplot)
library(effects)
library(DAMisc)
library(itsadug)
library(emmeans)

piska_df <- piska_df %>% mutate(abu_forest = rowSums(piska_df[,lesne]))
piska_df <- piska_df %>% mutate(rich_forest = rowSums(piska_df[,lesne] > 0))
piska_df <- piska_df %>% mutate(abu_farm = rowSums(piska_df[,farmlandowe]))
piska_df <- piska_df %>% mutate(rich_farm = rowSums(piska_df[,farmlandowe] > 0))
View(piska_df)


#filtering out species with total sum <50
modelledspecies <- c()
for(i in 12:length(colnames(piska_df))) { if(is.numeric(piska_df[,i]) == TRUE) { if(sum((piska_df[,i])) > 50){
modelledspecies <- c(modelledspecies, colnames(piska_df)[i]) }}}
piska_df_gamms <- cbind(piska_df[,1:11], filter(piska_df[,modelledspecies]))


piska_df_gamms$rok_factor <- as.factor(piska_df_gamms$rok)
piska_df_gamms$srodowisko <- as.factor(piska_df_gamms$srodowisko)
piska_df_gamms$srodowisko <- relevel(piska_df_gamms$srodowisko, c("control"))
piska_df_gamms$plot <- as.factor(piska_df_gamms$ID)

#ENGLISH COLUMN NAMES

kolumny <- colnames(piska_df_gamms)
gatunki_klucz <- as.data.frame(gatunki_klucz)
colnames(gatunki_klucz) <- c("name", "name_latin")
kolumny <- as.data.frame(kolumny)
colnames(kolumny) <- "name"
kolumny_ang <- left_join(kolumny, gatunki_klucz, by="name")
kolumny_ang <- kolumny_ang[!duplicated(kolumny_ang$name), ]
kolumny_ang[42,2] <- "Phylloscopus sibilatrix"
kolumny_ang$total <- coalesce(kolumny_ang$name_latin, kolumny_ang$name)
kolumny_ang$total[3] <- "habitat" 
kolumny_ang$total[6] <- "y_after_dist"
kolumny_ang$total[49] <- "abundance_for"
kolumny_ang$total[50] <- "richness_for"
kolumny_ang$total[51] <- "abundance_farm"
kolumny_ang$total[52] <- "richness_farm"
kolumny_ang$total[53] <- "year_factor"
kolumny_ang$total[54] <- "pointcount_id"

colnames(piska_df_gamms) <- kolumny_ang$total
piska_df_gamms <- piska_df_gamms %>% relocate(abundance_for,richness_for,abundance_farm,richness_farm, .after = richness)
sp_alphabetical <- sort(colnames(piska_df_gamms[18:52]))
piska_df_gamms <- piska_df_gamms[ , c(colnames(piska_df_gamms[1:17]), sp_alphabetical, colnames(piska_df_gamms[53:54]))]

#THE TABLE CREATION


pdf("Appendix table 2 - model outputs.pdf", height=8, width=8)
for(i in 12:(length(colnames(piska_df_gamms)) - 2)) { 
  piska_df_gamms$habitat <- relevel(piska_df_gamms$habitat, c("control"))
  g1<-gam(piska_df_gamms[,i] ~ s(pointcount_id,bs="re") + 
            habitat*y_after_dist + s(year_factor,bs='re'), data=piska_df_gamms, family=poisson, method="ML")
  g2<- gam(piska_df_gamms[,i] ~ s(pointcount_id,bs="re") + 
             habitat + y_after_dist + s(year_factor,bs='re'), data=piska_df_gamms, family=poisson, method="ML")
  aicinter <- round(AIC(logLik(g1)), 2)
  aicnointer <- round(AIC(logLik(g2)), 2)
  if (AIC(logLik(g1)) > AIC(logLik(g2))) {
    g_good <- g2
  } else {
    g_good <- g1 }
  summary(g_good)
  output<-capture.output(summary(g_good), file=NULL,append=FALSE)
  output[1] <- colnames(piska_df_gamms)[i]
  output[4] <- paste("AIC interaction:", aicinter, "    ", "AIC w/o interaction:", aicnointer)
  plot.new()
  tt1 <- ttheme_default()
  grid.table(output)
  while(!par('page')) plot.new()
}
dev.off()

