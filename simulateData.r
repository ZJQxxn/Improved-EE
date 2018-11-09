source('blkDiag.r')

p<-3000
n<-p
for (b_num in c(2,5,10,50,100)){
mat<-balancedBlk(p,blk=b_num)
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
data<-generateSample(mat,n)
print(sprintf('Finishd generating %d blocks',b_num))
write.table(mat, file = sprintf("precision/balance_3000vars_%dblks.csv",b_num),row.names=FALSE,col.names=FALSE, sep=",")
write.table(data, file = sprintf("data/balance_3000vars_%dblks_3000samples.csv",b_num),row.names=FALSE,col.names=FALSE, sep=",")
print(sprintf('Finished writing %d blocks',b_num))
}
