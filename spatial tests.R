library(dplyr)
library(ncf)
library(mgcv)
library(cowplot)

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

#The initial test for(i in c(8:13,17:51)


plot_list <- list()

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
  set.seed(1)
  splines.auto <- ncf::spline.correlog(x=piska_df_gamms$N,y= piska_df_gamms$E,z=resid(g_good), xmax=0.1) #xmax is in 2/3 of all distances
  name <- paste(colnames(piska_df_gamms[i]))
  par(xpd = NA, # switch off clipping, necessary to always see axis labels
      bg = "transparent", # switch off background to avoid obscuring adjacent plots
      oma = c(2, 2, 0, 0)) # move plot to the right and up
  plot(splines.auto, main = name, ylim = c(-0.5, 0.5), yaxt = "n")
  axis(2, at=seq(-0.5, 0.5, by = 0.2))
  plotplease <- recordPlot()
  plotplease <- ggdraw(plotplease)
  j = i - 7
  plot_list[[j]] = plotplease
}

#Rarefied plots, nb distribution (14,15,16)

for(i in c(14,15,16)) { 
  piska_df_gamms$habitat <- relevel(piska_df_gamms$habitat, c("unmanaged"))
  g1<-gam(piska_df_gamms[,i] ~ s(pointcount_id,bs="re") + habitat*year_since_dist + s(year_factor,bs='re') + 
            s(N_scaled, E_scaled), data=piska_df_gamms, family=nb, method="ML")
  g2<- gam(piska_df_gamms[,i] ~ s(pointcount_id,bs="re") + habitat + year_since_dist + s(year_factor,bs='re') + 
             s(N_scaled, E_scaled), data=piska_df_gamms, family=nb, method="ML")
  aicinter <- round(AIC(logLik(g1)), 2)
  aicnointer <- round(AIC(logLik(g2)), 2)
  if (AIC(logLik(g1)) > AIC(logLik(g2))) {
    g_good <- g2
  } else {
    g_good <- g1 }
  print(summary(g_good))
  set.seed(1)
  splines.auto <- ncf::spline.correlog(x=piska_df_gamms$N,y= piska_df_gamms$E,z=resid(g_good), xmax=0.1) #xmax is in 2/3 of all distances
  name <- paste(colnames(piska_df_gamms[i]))
  par(xpd = NA, # switch off clipping, necessary to always see axis labels
      bg = "transparent", # switch off background to avoid obscuring adjacent plots
      oma = c(2, 2, 0, 0)) # move plot to the right and up
  plot(splines.auto, main = paste(name), ylim = c(-0.5, 0.5), yaxt = "n", xaxt = "n")
  axis(2, at=seq(-0.5, 0.5, by = 0.2), cex.axis=0.8)
  axis(1, cex.axis=0.8)
  plotplease2 <- recordPlot()
  plotplease2 <- ggdraw(plotplease2)
  j = i - 7
  plot_list[[j]] <- plotplease2
}

#autocorrelation corrected - scaled NE, automatic k (8,9,17,33)

for(i in c(8,9,17,33)) { 
  piska_df_gamms$habitat <- relevel(piska_df_gamms$habitat, c("unmanaged"))
  g1<-gam(piska_df_gamms[,i] ~ s(pointcount_id,bs="re") + habitat*year_since_dist + s(year_factor,bs='re') + 
            s(N_scaled, E_scaled), data=piska_df_gamms, family=poisson, method="ML")
  g2<- gam(piska_df_gamms[,i] ~ s(pointcount_id,bs="re") + habitat + year_since_dist + s(year_factor,bs='re') + 
             s(N_scaled, E_scaled), data=piska_df_gamms, family=poisson, method="ML")
  aicinter <- round(AIC(logLik(g1)), 2)
  aicnointer <- round(AIC(logLik(g2)), 2)
  if (AIC(logLik(g1)) > AIC(logLik(g2))) {
    g_good <- g2
  } else {
    g_good <- g1 }
  print(summary(g_good))
  set.seed(1)
  splines.auto <- ncf::spline.correlog(x=piska_df_gamms$N,y= piska_df_gamms$E,z=resid(g_good), xmax=0.1) #xmax is in 2/3 of all distances
  name <- paste(colnames(piska_df_gamms[i]))
  par(xpd = NA, # switch off clipping, necessary to always see axis labels
      bg = "transparent", # switch off background to avoid obscuring adjacent plots
      oma = c(2, 2, 0, 0)) # move plot to the right and up
  plot(splines.auto, main = paste(name, "corrected for autocorr."), ylim = c(-0.5, 0.5), yaxt = "n", xaxt = "n")
  axis(2, at=seq(-0.5, 0.5, by = 0.2), cex.axis=0.8)
  axis(1, cex.axis=0.8)
  plotplease2 <- recordPlot()
  plotplease2 <- ggdraw(plotplease2)
  j = i - 7
  plot_list[[j]] <- plot_grid(plot_list[[j]],plotplease2)
}

#autocorrelation corrected - scaled NE, k increased to 60, for(i in c(13,27))

for(i in c(13,27)) { 
  piska_df_gamms$habitat <- relevel(piska_df_gamms$habitat, c("unmanaged"))
  g1<-gam(piska_df_gamms[,i] ~ s(pointcount_id,bs="re") + habitat*year_since_dist + s(year_factor,bs='re') + 
            s(N_scaled, E_scaled,k=60, fx=TRUE), data=piska_df_gamms, family=poisson, method="ML")
  g2<- gam(piska_df_gamms[,i] ~ s(pointcount_id,bs="re") + habitat + year_since_dist + s(year_factor,bs='re') + 
             s(N_scaled, E_scaled,k=60, fx=TRUE), data=piska_df_gamms, family=poisson, method="ML")
  aicinter <- round(AIC(logLik(g1)), 2)
  aicnointer <- round(AIC(logLik(g2)), 2)
  if (AIC(logLik(g1)) > AIC(logLik(g2))) {
    g_good <- g2
  } else {
    g_good <- g1 }
  print(summary(g_good))
  set.seed(1)
  splines.auto <- ncf::spline.correlog(x=piska_df_gamms$N,y= piska_df_gamms$E,z=resid(g_good), xmax=0.1) #xmax is in 2/3 of all distances
  name <- paste(colnames(piska_df_gamms[i]))
  par(xpd = NA, # switch off clipping, necessary to always see axis labels
      bg = "transparent", # switch off background to avoid obscuring adjacent plots
      oma = c(2, 2, 0, 0)) # move plot to the right and up
  plot(splines.auto, main = paste(name, "corrected for autocorr."), ylim = c(-0.5, 0.5), yaxt = "n", xaxt = "n")
  axis(2, at=seq(-0.5, 0.5, by = 0.2), cex.axis=0.8)
  axis(1, cex.axis=0.8)
  plotplease2 <- recordPlot()
  plotplease2 <- ggdraw(plotplease2)
  j = i - 7
  plot_list[[j]] <- plot_grid(plot_list[[j]],plotplease2)
}

dev.off()
pdf("ncf-spatial.correlationplots_final.pdf")
pdf.options(width = 12, height = 9)
for (i in 1:length(plot_list)){
  print(plot_list[[i]])
}
dev.off()