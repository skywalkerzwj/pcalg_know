library(pcalg)
library(classGraph)
source("TEST111.R")

#file="EMT-35-translated.csv"
file="data2.csv"
# header<-readHeader(file)
# num_miRNA<-35
# miRNA<-header[1:num_miRNA]
# mRNA<-header[-(1:num_miRNA)]
# 
# edge<-queryTargetFile(miRNA,mRNA,"Predicted/TargetScan-translated.csv"); edgeTargetScan<-edgeTargetScan+t(edgeTargetScan);edgeTargetScan<-edgeTargetScan!=0;
# 
# 
# edgeTargetScan<-queryTargetFile(miRNA,mRNA,"Predicted/TargetScan-translated.csv"); edgeTargetScan<-edgeTargetScan+t(edgeTargetScan);edgeTargetScan<-edgeTargetScan!=0;
# edgeMirBase<-queryTargetFile(miRNA,mRNA,"Predicted/miRBase-translated.csv"); edgeMirBase<-edgeMirBase+t(edgeMirBase);edgeMirBase<-edgeMirBase!=0;
# edgeTransMir<-queryTargetFile(miRNA,mRNA,"Predicted/TransmiR-translated.csv"); edgeTransMir<-edgeTransMir+t(edgeTransMir); edgeTransMir<-edgeTransMir!=0
# edgePPI<-queryPPIFile(miRNA,mRNA,"Predicted/PPI.csv");edgePPI<-edgePPI+t(edgePPI);edgePPI<-edgePPI!=0;
# # edgeTransMir<-t(edgeTransMir)
# 
# edgeTarTransPPI<-edgeTargetScan+edgeTransMir+edgePPI; edgeTarTransPPI<-edgeTarTransPPI!=0;
# edgeTarTrans<-edgeTargetScan+edgeTransMir;  edgeTarTrans<-edgeTarTrans!=0
# edgeTarPPI<-edgeTargetScan+edgePPI;  edgeTarPPI<-edgeTarPPI!=0
# edgeTransPPI<-edgeTransMir+edgePPI;  edgeTransPPI<-edgeTransPPI!=0
# 
# 
# edgeMirTrans<-edgeMirBase+edgeTransMir;  edgeMirTrans<-edgeMirTrans!=0
# edgeMirPPI<-edgeMirBase+edgePPI;  edgeMirPPI<-edgeMirPPI!=0
# edgeMirPPITrans<-edgeMirTrans+edgePPI;  edegMirPPITrans<-edgeMirPPITrans!=0
# edgeAll<-edgeTarTransPPI+edgeMirBase;  edgeAll<-edgeAll!=0
# edgeTarMir<-edgeMirBase+edgeTargetScan;  edgeTarMir<-edgeTarMir!=0
# edgeMirTarPPI<-edgeTarMir+edgePPI;  edgeMirTarPPI<-edgeMirTarPPI!=0
# edgeMirTarTrans<-edgeTarMir+edgeTransMir;  edgeMirTarTrans<-edgeMirTarTrans!=0


dt<-Read(file,sep="")
stdData<-Standardise(dt)


# resultNone<-CausaleffectsWeijia(miRNA,mRNA,stdData,0.01,num_miRNA)
# 
# resultPPI<-CausaleffectsWeijia(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgePPI)
# resultTargetScan<-CausaleffectsWeijia(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeTargetScan)
# resultTransMir<-CausaleffectsWeijia(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeTransMir)
# resultMirBase<-CausaleffectsWeijia(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeMirBase)
# 
# resultTarTrans<-CausaleffectsWeijia(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeTarTrans)
# resultMirTrans<-CausaleffectsWeijia(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeMirTrans)
# resultTarMir<-CausaleffectsWeijia(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeTarMir)
# resultMirPPI<-CausaleffectsWeijia(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeMirPPI)
# resultTransPPI<-CausaleffectsWeijia(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeTransPPI)
# resultTarPPI<-CausaleffectsWeijia(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeTarPPI)
# 
# 
# resultMirPPITrans<-CausaleffectsWeijia(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeMirPPITrans)
# resultMirTarPPI<-CausaleffectsWeijia(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeMirTarPPI)
# resultMirTarTrans<-CausaleffectsWeijia(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeMirTarTrans)
# resultTarTransPPI<-CausaleffectsWeijia(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeTarTransPPI)
# 
# resultAll<-CausaleffectsWeijia(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeAll)
# 
# 
# save(resultNone,resultPPI,resultTargetScan,resultTransMir,resultMirBase,
#      resultTarTrans,resultMirTrans,resultTarMir,resultMirPPI,resultTransPPI,resultTarPPI,
#      resultMirPPITrans,resultMirTarPPI,resultMirTarTrans,resultTarTransPPI,
#      resultAll,file="ResultEMT35wPPIDec.Rdata")
# # save(resultTarTransWeijia,resultMirTransWeijia,resultAllWeijia,file="")
# 
# # resultNonePC<-CausaleffectsPC(miRNA,mRNA,stdData,0.01,num_miRNA)
# # resultMirBasePC<-CausaleffectsPC(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeMirBase)
# # resultTargetScanPC<-CausaleffectsPC(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeTargetScan)
# # resultTransMirPC<-CausaleffectsPC(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeTransMir)
# # resultAllPC<-CausaleffectsPC(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeAll)
# # resultTarTransPC<-CausaleffectsPC(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeTarTrans)
# # resultMirTransPC<-CausaleffectsPC(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeMirTrans)
# # resultTarMirsPC<-CausaleffectsPC(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeTarMir)
# 
# # ptm <- proc.time()
# # tempRelaxedTRUE<-udag2pdagRelaxed(skel, solve.confl=TRUE,verbose=TRUE)
# # tempRetry<-udag2pdagSpecial(skel)$pcObj
# # proc.time() - ptm
