#Hard-thresholding Function 
hardThreshold<-function(S,lambda,p){
    for (i in 1:p) {
       for(j in i:p){
           if(abs(S[i,j])<lambda){
               S[i,j] <- 0
               S[j,i] <- 0
           }
       }
    }
    return (S)
}

softThreshold<-function(S,lambda,p,eta=1){
    for (i in 1:p) {
	   for(j in i:p){
		   tmp<-0
           if(abs(S[i,j])>lambda){
               tmp<-sign(S[i,j])*(abs(S[i,j])-lambda)
           }
		   S[i,j]<-tmp
		   S[j,i]<-tmp
       }
    }
    return (S)
}

adapLasso<-function(S,lambda,p,eta=1){
    for (i in 1:p) {
	   for(j in i:p){
		   tmp<-0
		   delta<-(lambda^(eta+1))*(abs(S[i,j])^(-eta))
           if(abs(S[i,j])>delta){
               tmp<-sign(S[i,j])*(abs(S[i,j])-delta)
           }
		   S[i,j]<-tmp
		   S[j,i]<-tmp
       }
    }
    return (S)
}

SCAD<-function(S,lambda,p,alpha=3.7){
    for (i in 1:p) {
	   for(j in i:p){
		   tmp<-0
           if(abs(S[i,j])<=2*lambda){
			   if(abs(S[i,j])>lambda){
				  tmp<-sign(S[i,j])*(abs(S[i,j])-lambda)
			   }
           }
		   else if(abs(S[i,j])>alpha*lambda){
				tmp<-S[i,j]
		   }
		   else{
				tmp<-((alpha-1)*S[i,j]-sign(S[i,j])*alpha*lambda)/(alpha-2)
		   }
		   S[i,j]<-tmp
		   S[j,i]<-tmp
       }
    }
    return (S)
}