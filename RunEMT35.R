library(pcalg)
library(classGraph)
source("TEST111.R")

file="EMT-35-translated.csv"
header<-readHeader(file)
num_miRNA<-35
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


resultNone<-CausaleffectsThuc(miRNA,mRNA,stdData,0.01,num_miRNA,skel.method="stable")

resultPPI<-CausaleffectsThuc(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgePPI,skel.method="stable")
resultTargetScan<-CausaleffectsThuc(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeTargetScan,skel.method="stable")
resultTransMir<-CausaleffectsThuc(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeTransMir,skel.method="stable")
resultMirBase<-CausaleffectsThuc(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeMirBase,skel.method="stable")

resultTarTrans<-CausaleffectsThuc(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeTarTrans,skel.method="stable")
resultMirTrans<-CausaleffectsThuc(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeMirTrans,skel.method="stable")
resultTarMir<-CausaleffectsThuc(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeTarMir,skel.method="stable")
resultMirPPI<-CausaleffectsThuc(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeMirPPI,skel.method="stable")
resultTransPPI<-CausaleffectsThuc(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeTransPPI,skel.method="stable")
resultTarPPI<-CausaleffectsThuc(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeTarPPI,skel.method="stable")


resultMirPPITrans<-CausaleffectsThuc(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeMirPPITrans,skel.method="stable")
resultMirTarPPI<-CausaleffectsThuc(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeMirTarPPI,skel.method="stable")
resultMirTarTrans<-CausaleffectsThuc(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeMirTarTrans,skel.method="stable")
resultTarTransPPI<-CausaleffectsThuc(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeTarTransPPI,skel.method="stable")

resultAll<-CausaleffectsThuc(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeAll,skel.method="stable")


save(resultNone,resultPPI,resultTargetScan,resultTransMir,resultMirBase,
     resultTarTrans,resultMirTrans,resultTarMir,resultMirPPI,resultTransPPI,resultTarPPI,
     resultMirPPITrans,resultMirTarPPI,resultMirTarTrans,resultTarTransPPI,
     resultAll,file="ResultEMT35wPPIStable.Rdata")
# save(resultTarTransThuc,resultMirTransThuc,resultAllThuc,file="")

# resultNonePC<-CausaleffectsPC(miRNA,mRNA,stdData,0.01,num_miRNA)
# resultMirBasePC<-CausaleffectsPC(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeMirBase)
# resultTargetScanPC<-CausaleffectsPC(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeTargetScan)
# resultTransMirPC<-CausaleffectsPC(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeTransMir)
# resultAllPC<-CausaleffectsPC(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeAll)
# resultTarTransPC<-CausaleffectsPC(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeTarTrans)
# resultMirTransPC<-CausaleffectsPC(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeMirTrans)
# resultTarMirsPC<-CausaleffectsPC(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeTarMir)

# ptm <- proc.time()
# tempRelaxedTRUE<-udag2pdagRelaxed(skel, solve.confl=TRUE,verbose=TRUE)
# tempRetry<-udag2pdagSpecial(skel)$pcObj
# proc.time() - ptm
