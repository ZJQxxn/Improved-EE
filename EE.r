source('thresholding.r')

EE<-function(S,lambda,p,thr_func=hardThreshold){
    S_lambda<-thr_func(S,lambda,p)
    return (solve(S_lambda))
}