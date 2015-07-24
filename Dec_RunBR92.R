library(pcalg)
library(classGraph)
source("TEST111.R")

file="BR51-92m-translated.csv"
header<-readHeader(file)
num_miRNA<-92
miRNA<-header[1:num_miRNA]
mRNA<-header[-(1:num_miRNA)]

edgeCombined<-queryTargetFile(miRNA,mRNA,"Verified/combined-translated-verified.csv"); edgeCombined<-edgeCombined+t(edgeCombined);edgeCombined<-edgeCombined!=0;
edgeTargetScan<-queryTargetFile(miRNA,mRNA,"Predicted/TargetScan-translated.csv"); edgeTargetScan<-edgeTargetScan+t(edgeTargetScan);edgeTargetScan<-edgeTargetScan!=0;
edgeMirBase<-queryTargetFile(miRNA,mRNA,"Predicted/miRBase-translated.csv"); edgeMirBase<-edgeMirBase+t(edgeMirBase);edgeMirBase<-edgeMirBase!=0;
edgeTransMir<-queryTargetFile(miRNA,mRNA,"Predicted/TransmiR-translated.csv"); edgeTransMir<-edgeTransMir+t(edgeTransMir); edgeTransMir<-edgeTransMir!=0
edgePPI<-queryPPIFile(miRNA,mRNA,"Predicted/PPI.csv");edgePPI<-edgePPI+t(edgePPI);edgePPI<-edgePPI!=0;

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


resultNone<-CausaleffectsWeijia(miRNA,mRNA,stdData,0.01,num_miRNA)

resultPPI<-CausaleffectsWeijia(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgePPI)
resultTargetScan<-CausaleffectsWeijia(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeTargetScan)
resultTransMir<-CausaleffectsWeijia(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeTransMir)
resultMirBase<-CausaleffectsWeijia(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeMirBase)

resultTarTrans<-CausaleffectsWeijia(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeTarTrans)
resultMirTrans<-CausaleffectsWeijia(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeMirTrans)
resultTarMir<-CausaleffectsWeijia(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeTarMir)
resultMirPPI<-CausaleffectsWeijia(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeMirPPI)
resultTransPPI<-CausaleffectsWeijia(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeTransPPI)
resultTarPPI<-CausaleffectsWeijia(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeTarPPI)


resultMirPPITrans<-CausaleffectsWeijia(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeMirPPITrans)
resultMirTarPPI<-CausaleffectsWeijia(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeMirTarPPI)
resultMirTarTrans<-CausaleffectsWeijia(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeMirTarTrans)
resultTarTransPPI<-CausaleffectsWeijia(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeTarTransPPI)

resultAll<-CausaleffectsWeijia(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeAll)

save(resultNone,resultPPI,resultTargetScan,resultTransMir,resultMirBase,
     resultTarTrans,resultMirTrans,resultTarMir,resultMirPPI,resultTransPPI,resultTarPPI,
     resultMirPPITrans,resultMirTarPPI,resultMirTarTrans,resultTarTransPPI,
     resultAll,file="Dec_ResultBR92wPPI.Rdata")

countTargetScan<-pcAND(suffStat, indepTest=gaussCItest, p=ncol(stdData), alpha = 0.01,fixedEdges=edgeTargetScan)
countMirBase<-pcAND(suffStat, indepTest=gaussCItest, p=ncol(stdData), alpha = 0.01,fixedEdges=edgeMirBase)
countTransMir<-pcAND(suffStat, indepTest=gaussCItest, p=ncol(stdData), alpha = 0.01,fixedEdges=edgeTransMir)