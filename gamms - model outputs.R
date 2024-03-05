library(dplyr)
library(ggplot2)
library(mgcv)
library(grid)
library(gridExtra)


#datasets upload

piska_df <- read.table("clipboard",header=TRUE,sep="\t",na.strings="NA",dec=".",strip.white=TRUE) #piska database full
species_habitats <- read.table("clipboard",header=TRUE,sep="\t",na.strings="NA",dec=".",strip.white=TRUE) #species_habitats_EU
View(piska_df)
#dataset preparation

forestbirds <- as.vector(species_habitats$forestbirds.EU)
farmlandbirds <- as.vector(species_habitats$farmland.EU)[1:9]

piska_df <- piska_df %>% mutate(abundance_for = rowSums(piska_df[,forestbirds]))
piska_df <- piska_df %>% mutate(richness_for = rowSums(piska_df[,forestbirds] > 0))
piska_df <- piska_df %>% mutate(abundance_farm = rowSums(piska_df[,farmlandbirds]))
piska_df <- piska_df %>% mutate(richness_farm = rowSums(piska_df[,farmlandbirds] > 0))

#filtering out species with total sum <50
modelledspecies <- c()
for(i in 8:length(colnames(piska_df))) { if(is.numeric(piska_df[,i]) == TRUE) { if(sum((piska_df[,i]), na.rm=TRUE) > 50){
  modelledspecies <- c(modelledspecies, colnames(piska_df)[i]) }}}
piska_df_gamms <- cbind(piska_df[,1:7], filter(piska_df[,modelledspecies]))
View(piska_df_gamms)

#preparing variables
piska_df_gamms$year_factor <- as.factor(piska_df_gamms$year)
piska_df_gamms$habitat <- as.factor(piska_df_gamms$habitat)
piska_df_gamms$habitat <- relevel(piska_df_gamms$habitat, c("control"))
piska_df_gamms$pointcount_id <- as.factor(piska_df_gamms$ID)

#sorting the columns
piska_df_gamms <- piska_df_gamms %>% relocate(abundance_for,richness_for,abundance_farm,richness_farm, 
                                              richness.rarefied, richness.for.rarefied, richness.farm.rarefied, .after = richness)
sp_alphabetical <- sort(colnames(piska_df_gamms[17:51]))
piska_df_gamms <- piska_df_gamms[ , c(colnames(piska_df_gamms[1:16]), sp_alphabetical, colnames(piska_df_gamms[52:53]))]


##scaling the N,E variables separately for each habitat 

piska_df_gamms  <- transform(piska_df_gamms, N_scaled=ave(N, habitat, FUN=scale))
piska_df_gamms  <- transform(piska_df_gamms, E_scaled=ave(E, habitat, FUN=scale))

View(piska_df_gamms)

table_list <- list()

#original models (no spatial autocorrelation)

for(i in c(8:13,17:51)) { 
  piska_df_gamms$habitat <- relevel(piska_df_gamms$habitat, c("unmanaged"))
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
  if (i==45) {
    output[30] <- "CAUTION: estimates for control unreliable: inflated standard errors (due to 0 observations in control)"
    t1 <- ttheme_default(core=list(
      fg_params=list(fontface=c(rep("plain", 29), "bold.italic")),
      bg_params = list(fill=c(rep(c("grey95", "grey90"),
                                  length.out=29), "firebrick1"),
                       alpha = rep(c(1,0.5), each=30))
    ))
    plot.new()
    grid.table(output,theme = t1)
  } else {plot.new()
    grid.table(output)}
  model_sum <- recordPlot()
  j = i - 7
  table_list[[j]] <- model_sum
}

#Rarefied plots, nb distribution (14,15,16)

for(i in c(14,15,16)) { 
  piska_df_gamms$habitat <- relevel(piska_df_gamms$habitat, c("unmanaged"))
  g1<-gam(piska_df_gamms[,i] ~ s(pointcount_id,bs="re") + 
            habitat*year_since_dist + s(year_factor,bs='re'), data=piska_df_gamms, family=nb, method="ML")
  g2<- gam(piska_df_gamms[,i] ~ s(pointcount_id,bs="re") + 
             habitat + year_since_dist + s(year_factor,bs='re'), data=piska_df_gamms, family=nb, method="ML")
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
  if (i==45) {
    output[30] <- "CAUTION: estimates for control unreliable: inflated standard errors (due to 0 observations in control)"
    t1 <- ttheme_default(core=list(
      fg_params=list(fontface=c(rep("plain", 29), "bold.italic")),
      bg_params = list(fill=c(rep(c("grey95", "grey90"),
                                  length.out=29), "firebrick1"),
                       alpha = rep(c(1,0.5), each=30))
    ))
    plot.new()
    grid.table(output,theme = t1)
  } else {plot.new()
    grid.table(output)}
  model_sum <- recordPlot()
  j = i - 7
  table_list[[j]] <- model_sum
}

#autocorrelation corrected - scaled NE, automatic k

for(i in c(8,9,17,33)) { 
  piska_df_gamms$habitat <- relevel(piska_df_gamms$habitat, c("unmanaged"))
  g1<-gam(piska_df_gamms[,i] ~ s(pointcount_id,bs="re") + habitat*year_since_dist + 
            s(year_factor,bs='re') + s(N_scaled, E_scaled), data=piska_df_gamms, family=poisson, method="ML")
  g2<- gam(piska_df_gamms[,i] ~ s(pointcount_id,bs="re") + habitat + year_since_dist + 
             s(year_factor,bs='re') + s(N_scaled, E_scaled), data=piska_df_gamms, family=poisson, method="ML")
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
  if (i==17) {
    output[30] <- "CAUTION: estimates for control unreliable: inflated standard errors (due to 0 observations in control)"
    t1 <- ttheme_default(core=list(
      fg_params=list(fontface=c(rep("plain", 29), "bold.italic")),
      bg_params = list(fill=c(rep(c("grey95", "grey90"),
                                  length.out=29), "firebrick1"),
                       alpha = rep(c(1,0.5), each=30))
    ))
    plot.new()
    grid.table(output,theme = t1)
  } else {plot.new()
    grid.table(output)}
  model_sum <- recordPlot()
  j = i - 7
  table_list[[j]] <- model_sum
}

#autocorrelation corrected - scaled NE, k increased

for(i in c(13,27)) { 
  piska_df_gamms$habitat <- relevel(piska_df_gamms$habitat, c("unmanaged"))
  g1<-gam(piska_df_gamms[,i] ~ s(pointcount_id,bs="re") + habitat*year_since_dist + s(year_factor,bs='re') + 
            s(N_scaled, E_scaled, k=60,fx=TRUE), data=piska_df_gamms, family=poisson, method="ML")
  g2<- gam(piska_df_gamms[,i] ~ s(pointcount_id,bs="re") + habitat + year_since_dist + s(year_factor,bs='re') + 
             s(N_scaled, E_scaled, k=60,fx=TRUE), data=piska_df_gamms, family=poisson, method="ML")
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
  model_sum <- recordPlot()
  j = i - 7
  table_list[[j]] <- model_sum
}

dev.off()
pdf("model_outputs_spatiallycorrected_rarefied.pdf")
pdf.options(width = 14, height = 10)
print("supplementary table XX")
for (i in 1:length(table_list)){
  print(table_list[[i]])
}
dev.off()
