library(iNEXT)

piska_df <- read.table("clipboard",header=TRUE,sep="\t",na.strings="NA",dec=".",strip.white=TRUE)
View(piska_df)

mat1<-matrix(nrow=742,ncol=1)
for (i in 1:742){
  print(i)
  eg0<-as.numeric(piska_df[i,10:104])
  gga<-list(eg0=eg0)
  ut1<-iNEXT(gga,q=0,datatype="abundance",knots=50,nboot=100,endpoint=20)
  mat1[i,1]<-ut1$iNextEst$size_based[20,5]}
piska_df$richness.raref20<-as.numeric(mat1)

mat2<-matrix(nrow=742,ncol=1)
for (i in 1:742){
  print(i)
  eg0<-as.numeric(piska_df[i,forestbirds])# forestbirds
  gga<-list(eg0=eg0)
  ut1<-iNEXT(gga,q=0,datatype="abundance",knots=50,nboot=100,endpoint=20)
  mat2[i,1]<-ut1$iNextEst$size_based[20,5]}
piska_df$richness.raref.forestbirds20<-as.numeric(mat2)

mat3<-matrix(nrow=742,ncol=1)
for (i in 1:742){
  print(i)
  eg0<-as.numeric(piska_df[i,farmlandbirds])# all species
  if(sum(eg0>0)>0){
    gga<-list(eg0=eg0)
    ut1<-iNEXT(gga,q=0,datatype="abundance",knots=50,nboot=100,endpoint=20)
    mat3[i,1]<-ut1$iNextEst$size_based[20,5]
  }else{
    mat3[i,1]<-NA
  }}
piska_df$richness.raref.farmlandbirds20<-as.numeric(mat3)
