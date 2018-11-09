source('blkDiag.r')
source('EE.r')
source('evaluation.r')
source('improved_EE.r')
library('magic')
library('matlab')



p<-10
thr<-0.275
mat<-adiag(sparseMatrix(5),sparseMatrix(5))

noise<-eye(p)
noise<-generateSample(noise,p)
noise<-noise/abs(max(abs(noise)))
noise<-noise*0.1
"
#Permuting
for(i in 1:5){
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

data<-generateSample(mat,2*p)
for (col in 1:p) {
   data[,col]<-data[,col]-mean(data[,col])
}
#Get sample covariance matrix and normalize it
S<-t(data)%*%(data)
S<-S/max(abs(S))

ee_model<-EE(S,thr,p)
imp_ee_model<-improvedEE(S,thr,p)
var_seq<-imp_ee_model[[2]]
imp_ee_model<-imp_ee_model[[1]]
row_ordered<-zeros(p)
for(i in 1:p){
    row_ordered[i,]<-mat[var_seq[i],]
}
col_ordered<-zeros(p)
for(j in 1:p){
    col_ordered[,j]<-row_ordered[,var_seq[j]]
}
reorder_mat<-col_ordered    

print('EE')
print(accuracy(mat,ee_model))
print('Improved EE')
print(accuracy(reorder_mat,imp_ee_model))
"