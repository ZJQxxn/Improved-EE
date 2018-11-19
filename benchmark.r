#Test on imbalanced graph
source("improved_EE.r")
source('splitEE.r')
source('EE.r')
source('evaluation.r')
source('blkDiag.r')
library('glasso')
library('QUIC')
library("igraph")
library("matlab")

thr<-0.23
count<-2


data<-as.matrix(read.csv('data/balance_3000vars_20blks_3000samples.csv',header=FALSE))
graph_structure<-as.matrix(read.csv('precision/balance_3000vars_20blks.csv',header=FALSE))

print('finished reading')

p=length(data[1,]) #Variable number
n=length(data[,1]) #Sample number
for (col in 1:p) {
   data[,col]<-data[,col]-mean(data[,col])
}
#Get sample covariance matrix and normalize it
S<-t(data)%*%(data)
rm(data)
S<-S/max(abs(S))

print('finished normalizing')

print('Split EE:')
time<-proc.time()-proc.time()
for(i in 1:count){
    ptm<-proc.time()
    spee_model<-splitEE(S,thr,p)
    time<-time+(proc.time()-ptm)
}
time<-time/count
print(time)
spee_model<-spee_model/max(abs(spee_model))
print(accuracy(graph_structure,spee_model))


print('Improved EE:')
print(system.time(imp_ee_model<-improvedEE(S,thr,p,core_num=1)))
var_seq<-imp_ee_model[[2]]
imp_ee_model<-imp_ee_model[[1]]
reorder_graph<-graph_structure[var_seq,]
reorder_graph<-reorder_graph[,var_seq]
imp_ee_model<-imp_ee_model/max(abs(imp_ee_model))
print(accuracy(reorder_graph,imp_ee_model))


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

"
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