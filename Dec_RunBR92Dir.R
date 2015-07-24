library(pcalg)
library(classGraph)
source("TEST111.R")

file="BR51-92m-translated.csv"
header<-readHeader(file)
num_miRNA<-92
miRNA<-header[1:num_miRNA]
mRNA<-header[-(1:num_miRNA)]

edgeTransMir<-queryTargetFile(miRNA,mRNA,"Predicted/TransmiR-translated.csv");
dirTransMir<-t(edgeTransMir); edgeTransMir<-edgeTransMir+t(edgeTransMir);
edgeTargetScan<-queryTargetFile(miRNA,mRNA,"Predicted/TargetScan-translated.csv"); edgeTargetScan<-edgeTargetScan+t(edgeTargetScan);
edgeMirBase<-queryTargetFile(miRNA,mRNA,"Predicted/miRBase-translated.csv"); edgeMirBase<-edgeMirBase+t(edgeMirBase)

edgeMirTrans<-edgeTransMir+edgeMirBase; 
edgeTarTrans<-edgeTargetScan+edgeTransMir;


edgeTargetScan<-edgeTargetScan!=0; edgeMirBase<-edgeMirBase!=0;edgeMirTrans<-edgeMirTrans!=0;edgeTarTrans<-edgeTarTrans!=0


dt<-Read(file)
stdData<-Standardise(dt)

resultTarDTrans<-CausaleffectsDirectionsWeijia(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeTarTrans,fixedDirections=dirTransMir,skel.method="stable")
resultMirDTrans<-CausaleffectsDirectionsWeijia(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeMirTrans,fixedDirections=dirTransMir,skel.method="stable")

# resultMirTarDTrans<-CausaleffectsDirectionsWeijia(miRNA,mRNA,stdData,0.01,num_miRNA,fixedEdge=edgeMirTar,fixedDirections=edgeTransMir)


save(resultTarDTrans,resultMirDTrans,file="Dec_ResultBR92DirStable.Rdata")
