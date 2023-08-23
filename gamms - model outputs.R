
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
  grid.table(output)
  while(!par('page')) plot.new()
}
dev.off()