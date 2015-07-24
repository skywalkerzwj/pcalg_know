library(pcalg)
library(classGraph)
source("TEST111.R")

file="BR51-92m-translated.csv"
header<-readHeader(file)
num_miRNA<-92
miRNA<-header[1:num_miRNA]
mRNA<-header[-(1:num_miRNA)]

edgeTargetScan<-queryTargetFile(miRNA,mRNA,"Predicted/TargetScan-translated.csv")
edgeMirBase<-queryTargetFile(miRNA,mRNA,"Predicted/miRBase-translated.csv")
edgeTransMir<-queryTargetFile(miRNA,mRNA,"Predicted/TransmiR-translated.csv")
edgePPI<-queryPPIFile(miRNA,mRNA,"Predicted/PPI.csv")
# edgeTransMir<-t(edgeTransMir)

edgeTarTransPPI<-edgeTargetScan+edgeTransMir+edgePPI; edgeTarTransPPI<-edgeTarTransPPI!=0
edgeTarTrans<-edgeTargetScan+edgeTransMir;  edgeTarTrans<-edgeTarTrans!=0
edgeTarPPI<-edgeTargetScan+edgePPI;  edgeTarPPI<-edgeTarPPI!=0
edgeTransPPI<-edgeTransMir+edgePPI;  edgeTransPPI<-edgeTransPPI!=0


edgeMirTrans<-edgeMirBase+edgeTransMir;  edgeMirTrans<-edgeMirTrans!=0
edgeMirPPI<-edgeMirBase+edgePPI;  edgeMirPPI<-edgeMirPPI!=0
edgeMirPPITrans<-edgeMirTrans+edgePPI;  edegMirPPITrans<-edgeMirPPITrans!=0
edgeAll<-edgeTarTransPPI+edgeMirBase;  edgeAll<-edgeAll!=0
edgeTarMir<-edgeMirBase+edgeTargetScan;  edgeTarMir<-edgeTarMir!=0
edgeMirTarPPI<-edgeTarMir+edgePPI;  edgeMirTarPPI<-edgeMirTarPPI!=0
edgeMirTarTrans<-edgeTarMir+edgeTransMir;  edgeMirTarTrans<-edgeMirTarTrans!=0


dt<-Read(file)
stdData<-Standardise(dt)


resultNone<-CausaleffectsThuc(miRNA,mRNA,stdData,0.01,num_miRNA)

resultPPI<-CausaleffectsThuc(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgePPI)
resultTargetScan<-CausaleffectsThuc(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeTargetScan)
resultTransMir<-CausaleffectsThuc(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeTransMir)
resultMirBase<-CausaleffectsThuc(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeMirBase)

resultTarTrans<-CausaleffectsThuc(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeTarTrans)
resultMirTrans<-CausaleffectsThuc(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeMirTrans)
resultTarMir<-CausaleffectsThuc(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeTarMir)
resultMirPPI<-CausaleffectsThuc(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeMirPPI)
resultTransPPI<-CausaleffectsThuc(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeTransPPI)
resultTarPPI<-CausaleffectsThuc(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeTarPPI)


resultMirPPITrans<-CausaleffectsThuc(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeMirPPITrans)
resultMirTarPPI<-CausaleffectsThuc(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeMirTarPPI)
resultMirTarTrans<-CausaleffectsThuc(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeMirTarTrans)
resultTarTransPPI<-CausaleffectsThuc(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeTarTransPPI)

resultAll<-CausaleffectsThuc(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeAll)


save(resultNone,resultPPI,resultTargetScan,resultTransMir,resultMirBase,
     resultTarTrans,resultMirTrans,resultTarMir,resultMirPPI,resultTransPPI,resultTarPPI,
     resultMirPPITrans,resultMirTarPPI,resultMirTarTrans,resultTarTransPPI,
     resultAll,file="ResultBR92wPPI.Rdata")