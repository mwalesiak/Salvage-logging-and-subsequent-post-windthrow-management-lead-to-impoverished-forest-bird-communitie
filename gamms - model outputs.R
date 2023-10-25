library(dplyr)
library(ggplot2)
library(mgcv)
library(grid)
library(gridExtra)

#datasets upload

piska_df <- read.table("clipboard",header=TRUE,sep="\t",na.strings="NA",dec=".",strip.white=TRUE) #piska database full
species_habitats <- read.table("clipboard",header=TRUE,sep="\t",na.strings="NA",dec=".",strip.white=TRUE) #species_habitats_EU

#dataset preparation

forestbirds <- as.vector(species_habitats$forestbirds.EU)
farmlandbirds <- as.vector(species_habitats$farmland.EU)[1:9]

piska_df <- piska_df %>% mutate(abundance_for = rowSums(piska_df[,forestbirds]))
piska_df <- piska_df %>% mutate(richness_for = rowSums(piska_df[,forestbirds] > 0))
piska_df <- piska_df %>% mutate(abundance_farm = rowSums(piska_df[,farmlandbirds]))
piska_df <- piska_df %>% mutate(richness_farm = rowSums(piska_df[,farmlandbirds] > 0))

#filtering out species with total sum <50
modelledspecies <- c()
for(i in 6:length(colnames(piska_df))) { if(is.numeric(piska_df[,i]) == TRUE) { if(sum((piska_df[,i])) > 50){
  modelledspecies <- c(modelledspecies, colnames(piska_df)[i]) }}}
piska_df_gamms <- cbind(piska_df[,1:5], filter(piska_df[,modelledspecies]))

#preparing variables
piska_df_gamms$year_factor <- as.factor(piska_df_gamms$year)
piska_df_gamms$habitat <- as.factor(piska_df_gamms$habitat)
piska_df_gamms$habitat <- relevel(piska_df_gamms$habitat, c("control"))
piska_df_gamms$pointcount_id <- as.factor(piska_df_gamms$ID)

#sorting the columns
piska_df_gamms <- piska_df_gamms %>% relocate(abundance_for,richness_for,abundance_farm,richness_farm, .after = richness)
sp_alphabetical <- sort(colnames(piska_df_gamms[12:46]))
piska_df_gamms <- piska_df_gamms[ , c(colnames(piska_df_gamms[1:11]), sp_alphabetical, colnames(piska_df_gamms[47:48]))]


pdf("Supplementary table 2 - model outputs.pdf", height=8, width=8)
for(i in 6:(length(colnames(piska_df_gamms)) - 2)) { 
  piska_df_gamms$habitat <- relevel(piska_df_gamms$habitat, c("control"))
  g1<-gam(piska_df_gamms[,i] ~ s(pointcount_id,bs="re") + 
            habitat*year_since_dist + s(year_factor,bs='re'), data=piska_df_gamms, family=poisson, method="ML")
  g2<- gam(piska_df_gamms[,i] ~ s(pointcount_id,bs="re") + 
             habitat + year_since_dist + s(year_factor,bs='re'), data=piska_df_gamms, family=poisson, method="ML")
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
  grid.table(output)
  while(!par('page')) plot.new()
}
dev.off()