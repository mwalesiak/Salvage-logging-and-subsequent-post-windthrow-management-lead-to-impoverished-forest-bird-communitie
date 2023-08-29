piska_df_mean <- read.table("clipboard",header=TRUE,sep="\t",na.strings="NA",dec=".",strip.white=TRUE)
library(ecotraj)
library(vegan)
library(dplyr)

species <- colnames(piska_df_mean[4:ncol(piska_df_mean)])
sites <- piska_df_mean$habitat
years <- piska_df_mean$year_since_dist
community <- select(piska_df_mean, all_of(species))
community_dist <- vegdist(community, method = "bray")

pdf("permanova_plot.pdf",width=8,height=6)
trajectoryPCoA(community_dist, sites, years, traj.colors = c("black","blue", "yellow"), lwd = 2,
               survey.labels = T)
dev.off()

trajectoryConvergence(community_dist, sites, years, symmetric = TRUE)
trajectoryLengths(community_dist, sites, years)

