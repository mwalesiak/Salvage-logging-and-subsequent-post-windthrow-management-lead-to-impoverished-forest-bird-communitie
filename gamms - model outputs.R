par(mfrow=c(1,2))
for(i in 12:14) { if(sum((piska_df[,i])) > 50){ 
  g1<-gam(piska_df[,i] ~ s(plot,bs="re") + #poissons should be fine judging on distribution (except richness & abu)
            srodowisko*lata_po_zab + s(rok_factor,bs='re'), data=piska_df, family=poisson, method="ML")
  g2<- gam(piska_df[,i] ~ s(plot,bs="re") + #poissons should be fine judging on distribution (except richness & abu)
             srodowisko + lata_po_zab + s(rok_factor,bs='re'), data=piska_df, family=poisson, method="ML")
  aicinter <- round(AIC(logLik(g1)), 2)
  aicnointer <- round(AIC(logLik(g2)), 2)
  if (AIC(logLik(g1)) > AIC(logLik(g2))) {
    g_good <- g2
  } else {
    g_good <- g1
  }
  summary(g_good)
  output<-capture.output(summary(g_good), file=NULL,append=FALSE)
  output[4] <- colnames(piska_df)[i]
  output[5] <- paste("AIC internakcja:", aicinter, "    ", "AIC no inter:", aicnointer)
  plot.new()
  grid.table(output)}
  else {
    plot.new() 
    output2 <- "malo ziomkow to jest ponizej 50"
    output3<- as.data.frame(output2, colnames(piska_df)[i])
    grid.table(output3)}
  while(!par('page')) plot.new()
}
dev.off()

