source('blkDiag.r')

blk<-20
for (p in c(20000)){
n<-p
mat<-balancedBlk(p,blk)
print(sprintf('Finishd generating mat'))
#Shuffle matrix
for(i in 1:p){
    a<-floor(runif(1,min=1,max=p))
    b<-floor(runif(1,min=1,max=p))
    #Permute row
    tmp<-mat[a,]
    mat[a,]<-mat[b,]
    mat[b,]<-tmp
    #Permute column
    tmp<-mat[,a]
    mat[,a]<-mat[,b]
    mat[,b]<-tmp
}
print(sprintf('Finishd shuffling mat'))
data<-generateSample(mat,n)
print(sprintf('Finishd generating samples'))
write.table(data, file = sprintf("data/balance_%dvars_20blks_%dsamples.csv",p,n),row.names=FALSE,col.names=FALSE, sep=",")
print('Finished writing samples')
write.table(mat, file = sprintf("precision/balance_%dvars_20blks.csv",p),row.names=FALSE,col.names=FALSE, sep=",")
print(sprintf('Finished writing precision'))
}
