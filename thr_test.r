source("improved_EE.r")
source('EE.r')
source('evaluation.r')
library("igraph")
library("matlab")

thr<-0.23

data<-as.matrix(read.csv('data/balance_2500vars_20blks_5000samples.csv',header=FALSE))
graph_structure<-as.matrix(read.csv('precision/balance_2500vars_20blks.csv',header=FALSE))
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
"
p<-5000
graph_structure<-as.matrix(read.csv('precision/balance_5000vars_20blks.csv',header=FALSE))
S<-as.matrix(read.csv('cov/S_5000.csv',header=FALSE))
print('finished reading')
"
print('Hard:')
imp_ee_model<-improvedEE(S,thr,p,core_num=1)
var_seq<-imp_ee_model[[2]]
imp_ee_model<-imp_ee_model[[1]]
reorder_graph<-graph_structure[var_seq,]
reorder_graph<-reorder_graph[,var_seq]
imp_ee_model<-imp_ee_model/max(abs(imp_ee_model))
print(accuracy(reorder_graph,imp_ee_model))

print('Soft:')
imp_ee_model<-improvedEE(S,thr,p,thr=softThreshold,core_num=1)
var_seq<-imp_ee_model[[2]]
imp_ee_model<-imp_ee_model[[1]]
reorder_graph<-graph_structure[var_seq,]
reorder_graph<-reorder_graph[,var_seq]
imp_ee_model<-imp_ee_model/max(abs(imp_ee_model))
print(accuracy(reorder_graph,imp_ee_model))

print('Adap. Lasso:')
imp_ee_model<-improvedEE(S,thr,p,thr=adapLasso,core_num=1)
var_seq<-imp_ee_model[[2]]
imp_ee_model<-imp_ee_model[[1]]
reorder_graph<-graph_structure[var_seq,]
reorder_graph<-reorder_graph[,var_seq]
imp_ee_model<-imp_ee_model/max(abs(imp_ee_model))
print(accuracy(reorder_graph,imp_ee_model))

print('SCAD:')
imp_ee_model<-improvedEE(S,thr,p,thr=SCAD,core_num=1)
var_seq<-imp_ee_model[[2]]
imp_ee_model<-imp_ee_model[[1]]
reorder_graph<-graph_structure[var_seq,]
reorder_graph<-reorder_graph[,var_seq]
imp_ee_model<-imp_ee_model/max(abs(imp_ee_model))
print(accuracy(reorder_graph,imp_ee_model))

