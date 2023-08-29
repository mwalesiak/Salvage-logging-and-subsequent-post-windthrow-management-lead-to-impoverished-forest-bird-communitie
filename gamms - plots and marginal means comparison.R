library(dplyr)
library(emmeans)
library(ggplot2)
library(mgcv)
library(grid)
library(gridExtra)
library(ggeffects)

#dataset preparation

piska_df <- read.table("clipboard",header=TRUE,sep="\t",na.strings="NA",dec=".",strip.white=TRUE)
View(piska_df)
species_habitats <- read.table("clipboard",header=TRUE,sep="\t",na.strings="NA",dec=".",strip.white=TRUE)
forestbirds <- as.vector(species_habitats$forestbirds.EU.)
farmlandbirds <- as.vector(species_habitats$farmland.EU.)[1:9]
piska_df <- piska_df %>% mutate(abu_forest = rowSums(piska_df[,forestbirds]))
piska_df <- piska_df %>% mutate(rich_forest = rowSums(piska_df[,forestbirds] > 0))
piska_df <- piska_df %>% mutate(abu_farm = rowSums(piska_df[,farmlandbirds]))
piska_df <- piska_df %>% mutate(rich_farm = rowSums(piska_df[,farmlandbirds] > 0))

piska_df_gamms <- piska_df

piska_df_gamms$year_factor <- as.factor(piska_df_gamms$year)
piska_df_gamms$habitat <- as.factor(piska_df_gamms$habitat)
piska_df_gamms$habitat <- relevel(piska_df_gamms$habitat, c("control"))
piska_df_gamms$plot <- as.factor(piska_df_gamms$ID)

#GRAPHS FOR SPECIES & TOTAL ABUNDANCE & RICHNESS
            #(later in the publication we presented only plots for species, 
              #for which any difference between managed windthrow/unmanaged windthrow was confirmed by GAMMs)
      

plot_list <- list()

for (i in c(12:108)) { if(sum((piska_df[,i])) > 50) {
  g1<-gam(piska_df_gamms[,i] ~ s(plot,bs="re") + 
            habitat*year_since_dist + s(year_factor,bs='re'), data=piska_df_gamms, family=poisson, method="ML")
  g2<- gam(piska_df_gamms[,i] ~ s(plot,bs="re") + 
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
  contrasts <- contrast(emmeans(ref_grid(g_good, at = list(year_since_dist = seq(5, 19, length = 30)), regrid = "response"), "habitat", by = "year_since_dist"),
                        "trt.vs.ctrl", infer = c(TRUE, FALSE))
  contrasts <- as.data.frame(contrasts)
  levels(contrasts$contrast) <- c("managed", "unmanaged")
  control <- setNames(data.frame(matrix(ncol = length(colnames(contrasts)), nrow = 30)), c(colnames(contrasts)))
  control$contrast <- "control"
  control$year_since_dist <- seq(5, 19, length = 30)
  control[is.na(control)] <- 0
  contrasts <- rbind(contrasts,control)
  colorpal <- c("#DDCC77", "#6699CC", "#117733")
  plotplease <- ggplot(contrasts, aes(x = year_since_dist, y = estimate, colour = factor(contrast))) +
    theme_classic() +
    geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL, fill = factor(contrast)), alpha = 0.15, colour = NA) +
    scale_colour_manual(values=colorpal) +
    scale_fill_manual(values=colorpal) +
    geom_line(linewidth=1.5) +
    ggtitle(paste(colnames(piska_df)[i]), model_chosen) + xlab("year") +
    xlab("years after disturbance") +
    ylab("difference from control")
  j = i-13
  plot_list[[j]] = plotplease
}  else {
  print("less than 50 observations")}}

plot_list <- plot_list[!sapply(plot_list,is.null)]
plot_list
glist <- lapply(plot_list, ggplotGrob)
ggsave("effectplotsmegatest.pdf", marrangeGrob(glist, nrow = 2, ncol = 2), height=8, width=8, dpi = 300) 




#Additional pairwise comparisons of marginal means conducted for those 4 richness/abundance models 
                              #in which adding year*habitat interaction increased model parsimony 
                              #(i.e. forest abundance, forest richness, farmland abundance, farmland richness). 
                              #We added periods of significant difference to the graphs manually (CorelDRAW), using results from below

emms_total <- data.frame()

for (i in 109:112) {
g1<-gam(piska_df_gamms[,i] ~ s(plot,bs="re") + 
          habitat*year_since_dist + s(year_factor,bs='re'), data=piska_df_gamms, family=poisson, method="ML")
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

write.csv(emms_total, "emms_total.csv")

plot <- ggplot(emms_total, aes(x = year_since_dist, y = p.value, colour = factor(contrast))) +
  theme_classic() +
  geom_point() +
  geom_line() +
  geom_hline(yintercept=0.05)
plot + facet_wrap(~ pair)