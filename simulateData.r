source('blkDiag.r')

p<-3000
blk<-20
n<-3000
#mat<-as.matrix(read.csv('precision/balance_3000vars_20blks.csv',header=FALSE))
for (coef in c(2,3,4,5)){

#mat<-gridGraph(vars,blk=20,row=10,col=15)
mat<-balancedBlk(p,blk,coef)
#Shuffle matrix
for(i in 1:vars){
    a<-floor(runif(1,min=1,max=vars))
    b<-floor(runif(1,min=1,max=vars))
    #Permute row
    tmp<-mat[a,]
    mat[a,]<-mat[b,]
    mat[b,]<-tmp
    #Permute column
    tmp<-mat[,a]
    mat[,a]<-mat[,b]
    mat[,b]<-tmp
}

data<-generateSample(mat,n)
print(sprintf('Finishd generating %d coef',coef))
write.table(data, file = sprintf("data/sparsity_3000vars_20blks_3000samples_%dsparse.csv",coef),row.names=FALSE,col.names=FALSE, sep=",")
write.table(mat, file = sprintf("precision/sparsity_3000vars_20blks_%dsparse.csv",coef),row.names=FALSE,col.names=FALSE, sep=",")
print(sprintf('Finished writing %d coef',coef))
}
