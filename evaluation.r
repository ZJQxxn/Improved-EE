accuracy<-function(real_mat,est_mat){
    dim<-length(real_mat[1,])
    real_edge<-0
    real_nonedge<-0
    tp<-0
    fp<-0
    for(i in 1:dim){
        j<-i
        while(j<=dim){
            if(real_mat[i,j]!=0){
                real_edge<-real_edge+1
                if(est_mat[i,j]!=0){
                    tp<-tp+1
                }
            }
            else if(real_mat[i,j]==0){
               real_nonedge<-real_nonedge+1
                if(est_mat[i,j]!=0){
                    fp<-fp+1
                }
            }
            j<-j+1
        }
    }
    tpr<-tp/real_edge
    fpr<-fp/real_nonedge
    frobenius<-norm(real_mat-est_mat,type="F")
    maximum<-norm(real_mat-est_mat,type="M")
    return (list(tpr=tpr,fpr=fpr,frobenius=frobenius,maximum=maximum))
}