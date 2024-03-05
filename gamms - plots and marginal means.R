library(dplyr)
library(emmeans)
library(ggplot2)
library(mgcv)
library(grid)
library(gridExtra)
library(ggeffects)

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


#all graphs for species & bird diversities
            #(later in the publication we presented only plots for selected species: 
              #those for which any difference between managed windthrow/unmanaged windthrow was confirmed by GAMMs)

plot_list <- list()

for(i in c(8:13,17:51)) {
  piska_df_gamms$habitat <- relevel(piska_df_gamms$habitat, c("control"))
  g1<-gam(piska_df_gamms[,i] ~ s(pointcount_id,bs="re") + 
            habitat*year_since_dist + s(year_factor,bs='re'), data=piska_df_gamms, family=poisson, method="ML")
  g2<- gam(piska_df_gamms[,i] ~ s(pointcount_id,bs="re") + 
             habitat + year_since_dist + s(year_factor,bs='re'), data=piska_df_gamms, family=poisson, method="ML")
  aicinter <- round(AIC(logLik(g1)), 2)
  aicnointer <- round(AIC(logLik(g2)), 2)
  if (AIC(logLik(g1)) > AIC(logLik(g2))) {
    g_good <- g2
    model_chosen <- "NO internaction env*years"
  } else {
    g_good <- g1
    model_chosen <- "WITH interaction env*years"
  }
  print(summary(g_good))
  print(colnames(piska_df_gamms)[i])
  emms <- emmeans(ref_grid(g_good, at = list(year_since_dist = seq(5, 19, length = 30)), regrid = "response"), "habitat", by = "year_since_dist")
  contrasts <- contrast(emms,
                        "trt.vs.ctrl", infer = c(TRUE, TRUE))
  contrasts <- as.data.frame(contrasts)
  levels(contrasts$contrast) <- c("managed", "unmanaged")
  control <- setNames(data.frame(matrix(ncol = length(colnames(contrasts)), nrow = 30)), c(colnames(contrasts)))
  control$contrast <- "control"
  control$year_since_dist <- seq(5, 19, length = 30)
  control[is.na(control)] <- 0
  contrasts <- rbind(contrasts,control)
  colorpal <- c("#F0E442", "#0072B1", "#00D198")
  plotplease <- ggplot(contrasts, aes(x = year_since_dist, y = estimate, colour = factor(contrast))) +
    theme_classic() +
    geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL, fill = factor(contrast)), alpha = 0.15, colour = NA) +
    scale_colour_manual(values=colorpal) +
    scale_fill_manual(values=colorpal) +
    geom_line(linewidth=1.5) +
    ggtitle(paste(colnames(piska_df_gamms)[i]), model_chosen) + xlab("year") +
    xlab("years after disturbance") +
    ylab("difference from control")
  j = i-7
  plot_list[[j]] = plotplease
} 

#Rarefied plots, nb distribution (14,15,16)

for(i in c(14,15,16)) {
  piska_df_gamms$habitat <- relevel(piska_df_gamms$habitat, c("control"))
  g1<-gam(piska_df_gamms[,i] ~ s(pointcount_id,bs="re") + 
            habitat*year_since_dist + s(year_factor,bs='re'), data=piska_df_gamms, family=nb, method="ML")
  g2<- gam(piska_df_gamms[,i] ~ s(pointcount_id,bs="re") + 
             habitat + year_since_dist + s(year_factor,bs='re'), data=piska_df_gamms, family=nb, method="ML")
  aicinter <- round(AIC(logLik(g1)), 2)
  aicnointer <- round(AIC(logLik(g2)), 2)
  if (AIC(logLik(g1)) > AIC(logLik(g2))) {
    g_good <- g2
    model_chosen <- "NO internaction env*years"
  } else {
    g_good <- g1
    model_chosen <- "WITH interaction env*years"
  }
  print(summary(g_good))
  print(colnames(piska_df_gamms)[i])
  emms <- emmeans(ref_grid(g_good, at = list(year_since_dist = seq(5, 19, length = 30)), regrid = "response"), "habitat", by = "year_since_dist")
  contrasts <- contrast(emms,
                        "trt.vs.ctrl", infer = c(TRUE, TRUE))
  contrasts <- as.data.frame(contrasts)
  levels(contrasts$contrast) <- c("managed", "unmanaged")
  control <- setNames(data.frame(matrix(ncol = length(colnames(contrasts)), nrow = 30)), c(colnames(contrasts)))
  control$contrast <- "control"
  control$year_since_dist <- seq(5, 19, length = 30)
  control[is.na(control)] <- 0
  contrasts <- rbind(contrasts,control)
  colorpal <- c("#F0E442", "#0072B1", "#00D198")
  plotplease <- ggplot(contrasts, aes(x = year_since_dist, y = estimate, colour = factor(contrast))) +
    theme_classic() +
    geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL, fill = factor(contrast)), alpha = 0.15, colour = NA) +
    scale_colour_manual(values=colorpal) +
    scale_fill_manual(values=colorpal) +
    geom_line(linewidth=1.5) +
    ggtitle(paste(colnames(piska_df_gamms)[i]), model_chosen) + xlab("year") +
    xlab("years after disturbance") +
    ylab("difference from control")
  j = i-7
  plot_list[[j]] = plotplease
} 

#autocorrelation corrected - scaled NE, automatic k

for(i in c(8,9,17,33)) {
  piska_df_gamms$habitat <- relevel(piska_df_gamms$habitat, c("control"))
  g2<- gam(piska_df_gamms[,i] ~ s(pointcount_id,bs="re") + habitat + year_since_dist + 
             s(year_factor,bs='re') + s(N_scaled, E_scaled), data=piska_df_gamms, family=poisson, method="ML")
    g_good <- g2 #we already tested previously that excluding interaction lowered AIC for this whole group
    model_chosen <- "NO internaction env*years"
  print(summary(g_good))
  print(colnames(piska_df_gamms)[i])
  contrasts <- contrast(emmeans(ref_grid(g_good, at = list(year_since_dist = seq(5, 19, length = 30)), regrid = "response"), "habitat", by = "year_since_dist"),
                        "trt.vs.ctrl", infer = c(TRUE, FALSE))
  contrasts <- as.data.frame(contrasts)
  levels(contrasts$contrast) <- c("managed", "unmanaged")
  control <- setNames(data.frame(matrix(ncol = length(colnames(contrasts)), nrow = 30)), c(colnames(contrasts)))
  control$contrast <- "control"
  control$year_since_dist <- seq(5, 19, length = 30)
  control[is.na(control)] <- 0
  contrasts <- rbind(contrasts,control)
  colorpal <- c("#F0E442", "#0072B1", "#00D198")
  plotplease <- ggplot(contrasts, aes(x = year_since_dist, y = estimate, colour = factor(contrast))) +
    theme_classic() +
    geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL, fill = factor(contrast)), alpha = 0.15, colour = NA) +
    scale_colour_manual(values=colorpal) +
    scale_fill_manual(values=colorpal) +
    geom_line(linewidth=1.5) +
    ggtitle(paste(colnames(piska_df_gamms)[i]), model_chosen) + xlab("year") +
    xlab("years after disturbance") +
    ylab("difference from control")
  j = i-7
  plot_list[[j]] = plotplease
} 

#autocorrelation corrected - scaled NE, k increased

for(i in c(13,27)) {
  piska_df_gamms$habitat <- relevel(piska_df_gamms$habitat, c("control"))
  g1<-gam(piska_df_gamms[,i] ~ s(pointcount_id,bs="re") + habitat*year_since_dist + s(year_factor,bs='re') + 
            s(N_scaled, E_scaled, k=60,fx=TRUE), data=piska_df_gamms, family=poisson, method="ML")
    g_good <- g1 #we already tested previously that including interaction lowered AIC for this whole group
    model_chosen <- "WITH interaction env*years"
  print(summary(g_good))
  print(colnames(piska_df_gamms)[i])
  contrasts <- contrast(emmeans(ref_grid(g_good, at = list(year_since_dist = seq(5, 19, length = 30)), regrid = "response"), "habitat", by = "year_since_dist", regrid = "response"),
                        method="trt.vs.ctrl", infer = c(TRUE, TRUE))
  contrasts <- as.data.frame(contrasts)
  levels(contrasts$contrast) <- c("managed", "unmanaged")
  control <- setNames(data.frame(matrix(ncol = length(colnames(contrasts)), nrow = 30)), c(colnames(contrasts)))
  control$contrast <- "control"
  control$year_since_dist <- seq(5, 19, length = 30)
  control[is.na(control)] <- 0
  contrasts <- rbind(contrasts,control)
  colorpal <- c("#F0E442", "#0072B1", "#00D198")
  plotplease <- ggplot(contrasts, aes(x = year_since_dist, y = estimate, colour = factor(contrast))) +
    theme_classic() +
    geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL, fill = factor(contrast)), alpha = 0.15, colour = NA) +
    scale_colour_manual(values=colorpal) +
    scale_fill_manual(values=colorpal) +
    geom_line(linewidth=1.5) +
    ggtitle(paste(colnames(piska_df_gamms)[i]), model_chosen) + xlab("year") +
    xlab("years after disturbance") +
    ylab("difference from control")
  j = i-7
  plot_list[[j]] = plotplease
} 


dev.off()
plot_list <- plot_list[!sapply(plot_list,is.null)]
plot_list
glist <- lapply(plot_list, ggplotGrob)
ggsave("gamm effectplots_contrasts.pdf", marrangeGrob(glist, nrow = 3, ncol = 3), height=12, width=12, dpi = 300) 








#Additional pairwise comparisons of marginal means conducted for those 4 richness/abundance models 
                              #in which adding year*habitat interaction was used in the model (increased their parsimony) 
                              #(i.e. forest abundance, forest richness, farmland abundance, farmland richness). 
                              #We added periods of significant difference to the graphs manually (CorelDRAW), using results from below

emms_total <- data.frame()

for (i in 10:12) {
  g1<-gam(piska_df_gamms[,i] ~ s(pointcount_id,bs="re") + habitat*year_since_dist + s(year_factor,bs='re'), data=piska_df_gamms, family=poisson, method="ML")
  g_good <- g1
  plotme<-ggemmeans(g_good,terms=c("year_since_dist","habitat"))
  plot(plotme)
  rg <- ref_grid(g_good, cov.keep = "year_since_dist", at = list("year_since_dist" = seq(min(piska_df$year_since_dist),max(piska_df$year_since_dist), length.out = 15))) 
  emms <- emmeans(rg, "habitat", "year_since_dist")
  pairs(emms)
  emms2 <- as.data.frame(pairs(emms))
  emms2$pair <- colnames(piska_df_gamms[i])
  emms_total <- rbind(emms_total, emms2)
}

for (i in 13) {
piska_df_gamms$habitat <- relevel(piska_df_gamms$habitat, c("control"))
g1<-gam(piska_df_gamms[,i] ~ s(pointcount_id,bs="re") + habitat*year_since_dist + s(year_factor,bs='re') + 
            s(N_scaled, E_scaled, k=60,fx=TRUE), data=piska_df_gamms, family=poisson, method="ML")
g_good <- g1
plotme<-ggemmeans(g_good,terms=c("year_since_dist","habitat"))
plot(plotme)
rg <- ref_grid(g_good, cov.keep = "year_since_dist", at = list("year_since_dist" = seq(min(piska_df$year_since_dist),max(piska_df$year_since_dist), length.out = 15))) 
emms <- emmeans(rg, "habitat", "year_since_dist")
pairs(emms)
emms2 <- as.data.frame(pairs(emms))
emms2$pair <- colnames(piska_df_gamms[i])
emms_total <- rbind(emms_total, emms2)
}

emms_total
write.csv(emms_total, "GAMMs - marginal means pairwise test results.csv")

#results from analysis above are presented in supplementary material. Plots below were used to easier identifies periods of differences visually

plot <- ggplot(emms_total, aes(x = year_since_dist, y = p.value, colour = factor(contrast))) +
  theme_classic() +
  geom_point() +
  geom_line() +
  geom_hline(yintercept=0.05)
plot + facet_wrap(~ pair)


###### additional plots for appendix, no contrast #####

plot_list <- list()

for(i in c(8:13,17:51)) {
  piska_df_gamms$habitat <- relevel(piska_df_gamms$habitat, c("control"))
  g1<-gam(piska_df_gamms[,i] ~ s(pointcount_id,bs="re") + 
            habitat*year_since_dist + s(year_factor,bs='re'), data=piska_df_gamms, family=poisson, method="ML")
  g2<- gam(piska_df_gamms[,i] ~ s(pointcount_id,bs="re") + 
             habitat + year_since_dist + s(year_factor,bs='re'), data=piska_df_gamms, family=poisson, method="ML")
  aicinter <- round(AIC(logLik(g1)), 2)
  aicnointer <- round(AIC(logLik(g2)), 2)
  if (AIC(logLik(g1)) > AIC(logLik(g2))) {
    g_good <- g2
    model_chosen <- "NO internaction env*years"
  } else {
    g_good <- g1
    model_chosen <- "WITH interaction env*years"
  }
  print(summary(g_good))
  print(colnames(piska_df_gamms)[i])
  emms <- emmeans(ref_grid(g_good, at = list(year_since_dist = seq(5, 19, length = 30))), "habitat", by = "year_since_dist", type="response")
  emms <- as.data.frame(emms)
  colorpal <- c("#00D198","#F0E442", "#0072B1")
  plotplease <- ggplot(emms, aes(x = year_since_dist, y = response, colour = factor(habitat))) +
    theme_classic() +
    geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL, fill = factor(habitat)), alpha = 0.15, colour = NA) +
    scale_colour_manual(values=colorpal) +
    scale_fill_manual(values=colorpal) +
    geom_line(linewidth=1.5) +
    ggtitle(paste(colnames(piska_df_gamms)[i]), model_chosen) + xlab("year") +
    xlab("years after disturbance") +
    ylab("response")
  j = i-7
  plot_list[[j]] = plotplease
} 

#Rarefied plots, nb distribution (14,15,16)

for(i in c(14,15,16)) {
  piska_df_gamms$habitat <- relevel(piska_df_gamms$habitat, c("control"))
  g1<-gam(piska_df_gamms[,i] ~ s(pointcount_id,bs="re") + 
            habitat*year_since_dist + s(year_factor,bs='re'), data=piska_df_gamms, family=nb, method="ML")
  g2<- gam(piska_df_gamms[,i] ~ s(pointcount_id,bs="re") + 
             habitat + year_since_dist + s(year_factor,bs='re'), data=piska_df_gamms, family=nb, method="ML")
  aicinter <- round(AIC(logLik(g1)), 2)
  aicnointer <- round(AIC(logLik(g2)), 2)
  if (AIC(logLik(g1)) > AIC(logLik(g2))) {
    g_good <- g2
    model_chosen <- "NO internaction env*years"
  } else {
    g_good <- g1
    model_chosen <- "WITH interaction env*years"
  }
  print(summary(g_good))
  print(colnames(piska_df_gamms)[i])
  emms <- emmeans(ref_grid(g_good, at = list(year_since_dist = seq(5, 19, length = 30))), "habitat", by = "year_since_dist", type="response")
  emms <- as.data.frame(emms)
  colorpal <- c("#00D198","#F0E442", "#0072B1")
  plotplease <- ggplot(emms, aes(x = year_since_dist, y = response, colour = factor(habitat))) +
    theme_classic() +
    geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL, fill = factor(habitat)), alpha = 0.15, colour = NA) +
    scale_colour_manual(values=colorpal) +
    scale_fill_manual(values=colorpal) +
    geom_line(linewidth=1.5) +
    ggtitle(paste(colnames(piska_df_gamms)[i]), model_chosen) + xlab("year") +
    xlab("years after disturbance") +
    ylab("response")
  j = i-7
  plot_list[[j]] = plotplease
} 


#autocorrelation corrected - scaled NE for(i in c(8,9,17,33))

for(i in c(8,9,17,33)) {
  piska_df_gamms$habitat <- relevel(piska_df_gamms$habitat, c("control"))
  g2<- gam(piska_df_gamms[,i] ~ s(pointcount_id,bs="re") + habitat + year_since_dist + 
             s(year_factor,bs='re') + s(N_scaled, E_scaled), data=piska_df_gamms, family=poisson, method="ML")
  g_good <- g2 #we already tested previously that excluding interaction lowered AIC for this whole group
  model_chosen <- "NO internaction env*years"
  print(summary(g_good))
  print(colnames(piska_df_gamms)[i])
  emms <- emmeans(ref_grid(g_good, at = list(year_since_dist = seq(5, 19, length = 30))), "habitat", by = "year_since_dist", type="response")
  emms <- as.data.frame(emms)
  colorpal <- c("#00D198","#F0E442", "#0072B1")
  plotplease <- ggplot(emms, aes(x = year_since_dist, y = response, colour = factor(habitat))) +
    theme_classic() +
    geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL, fill = factor(habitat)), alpha = 0.15, colour = NA) +
    scale_colour_manual(values=colorpal) +
    scale_fill_manual(values=colorpal) +
    geom_line(linewidth=1.5) +
    ggtitle(paste(colnames(piska_df_gamms)[i]), model_chosen) + xlab("year") +
    xlab("years after disturbance") +
    ylab("response")
  j = i-7
  plot_list[[j]] = plotplease
} 

#autocorrelation corrected - scaled NE, k increased for(i in c(13,27))

for(i in c(13,27)) {
  piska_df_gamms$habitat <- relevel(piska_df_gamms$habitat, c("control"))
  g1<-gam(piska_df_gamms[,i] ~ s(pointcount_id,bs="re") + habitat*year_since_dist + s(year_factor,bs='re') + 
            s(N_scaled, E_scaled, k=60,fx=TRUE), data=piska_df_gamms, family=poisson, method="ML")
  g_good <- g1 #we already tested previously that including interaction lowered AIC for this whole group
  model_chosen <- "WITH interaction env*years"
  print(summary(g_good))
  print(colnames(piska_df_gamms)[i])
  emms <- emmeans(ref_grid(g_good, at = list(year_since_dist = seq(5, 19, length = 30))), "habitat", by = "year_since_dist", type="response")
  emms <- as.data.frame(emms)
  colorpal <- c("#00D198","#F0E442", "#0072B1")
  plotplease <- ggplot(emms, aes(x = year_since_dist, y = response, colour = factor(habitat))) +
    theme_classic() +
    geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL, fill = factor(habitat)), alpha = 0.15, colour = NA) +
    scale_colour_manual(values=colorpal) +
    scale_fill_manual(values=colorpal) +
    geom_line(linewidth=1.5) +
    ggtitle(paste(colnames(piska_df_gamms)[i]), model_chosen) + xlab("year") +
    xlab("years after disturbance") +
    ylab("response")
  j = i-7
  plot_list[[j]] = plotplease
} 


dev.off()
plot_list <- plot_list[!sapply(plot_list,is.null)]
plot_list
glist <- lapply(plot_list, ggplotGrob)
ggsave("gamm effectplots_response.pdf", marrangeGrob(glist, nrow = 3, ncol = 3), height=12, width=12, dpi = 300) 


