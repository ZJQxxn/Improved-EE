accuracy<-function(real_mat,est_mat){
    dim<-length(real_mat[1,])
    tn<-0
    fn<-0
    tp<-0
    fp<-0
    for(i in 1:dim){
        for(j in 1:dim){
            if (est_mat[i,j]!=0){
                if (real_mat[i,j]!=0){
                    tp<-tp+1
                }
                else{
                    fp<-fp+1
                }
            }
            else{
                if (real_mat[i,j]!=0){
                    fn<-fn+1
                }
                else{
                    tn<-tn+1
                }
            }
        }
    }
    tpr<-tp/(tp+fn)
    fpr<-fp/(fp+tn)
    frobenius<-norm(real_mat-est_mat,type="F")
    maximum<-norm(real_mat-est_mat,type="M")
    return (list(tpr=tpr,fpr=fpr,frobenius=frobenius,maximum=maximum))
}