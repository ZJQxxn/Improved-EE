source("improved_EE.r")
source('EE.r')
source('evaluation.r')
source('blkDiag.r')
library('glasso')
library('QUIC')
library("igraph")
library("matlab")
library('clime')

#thr<-0.08
count<-1

thrList<-c(0.05,0.075,0.1,0.125,0.15,0.175,0.2)
#thrList<-c(0.2)


graph_structure<-as.matrix(read.csv(sprintf('precision/grid_3000vars_20blks.csv'),header=FALSE))
data<-as.matrix(read.csv(sprintf('data/grid_3000vars_20blks_6000samples.csv'),header=FALSE))
print(sprintf('Finished reading'))

p=length(data[1,]) #Variable number
n=length(data[,1]) #Sample number
for (col in 1:p) {
   data[,col]<-data[,col]-mean(data[,col])
}
#Get sample covariance matrix and normalize it
S<-t(data)%*%(data)
rm(data)
S<-S/max(abs(S))
print('Finished normalizing')

for (thr in thrList){
print(sprintf('----------------------------%f thr-----------------',thr))

print('Improved EE:')
time<-proc.time()-proc.time()
for(i in 1:count){
     ptm<-proc.time()
    imp_ee_model<-improvedEE(S,thr,p,core_num=1)
     time<-time+(proc.time()-ptm)
}
time<-time/count
print(time)
var_seq<-imp_ee_model[[2]]
imp_ee_model<-imp_ee_model[[1]]
reorder_graph<-graph_structure[var_seq,]
reorder_graph<-reorder_graph[,var_seq]
imp_ee_model<-imp_ee_model/max(abs(imp_ee_model))
print(accuracy(reorder_graph,imp_ee_model))

"
print('EE:')
time<-proc.time()-proc.time()
for(i in 1:count){
    ptm<-proc.time()
    ee_model<-EE(S,thr,p)
    time<-time+(proc.time()-ptm)
}
time<-time/count
print(time)
ee_model<-ee_model/max(abs(ee_model))
print(accuracy(graph_structure,ee_model))


print('GLASSO:')
time<-proc.time()-proc.time()
for(i in 1:count){
    ptm<-proc.time()
    glasso_model<-glasso(S,thr)
    time<-time+(proc.time()-ptm)
}
time<-time/count
print(time)
glasso_model<-glasso_model[[2]]
glasso_model<-glasso_model/max(abs(glasso_model))
print(accuracy(graph_structure,glasso_model))


print('QUIC:')
time<-proc.time()-proc.time()
for(i in 1:count){
    ptm<-proc.time()
    quic_model<-QUIC(S,thr)
    time<-time+(proc.time()-ptm)
}
time<-time/count
print(time)
quic_model<-quic_model[[1]]
quic_model<-quic_model/max(abs(quic_model))
print(accuracy(graph_structure,quic_model))
"
}
