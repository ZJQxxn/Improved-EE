# The implementation of improved EE.
# Pay attention that we don't consider about storage cost of these functions. 
#-------------------------------------------------------------------------------
library("igraph")
library('magic')
library("matlab")

library('foreach')
library('doParallel')

source('thresholding.r')




#Improved Elementary Wstimator
improvedEE<-function(S,lambda,p,thr_func=hardThreshold,core_num=1){
    #Initialization
    Omega<-zeros(p)
    #Thresholding on covariance matrix
    S_lambda<-thr_func(S,lambda,p)
    #Get connected component
    graph<-graph_from_adjacency_matrix(S_lambda,mode="undirected",weighted=TRUE)
    comps<-components(graph)
    num<-comps[[3]]
    bins<-comps[[1]]
    csize<-comps[[2]]
    #Permuting the matrix so that it is block diagonal
    var_seq<-zeros(p,2)
    var_seq[,1]<-1:p
    var_seq[,2]<-bins
    var_seq<-var_seq[order(var_seq[,2]),]
    var_seq<-var_seq[,1]
    S_lambda<-S_lambda[var_seq,]
    S_lambda<-S_lambda[,var_seq]
    if(core_num>1){
        #Parallelizing
        each_inv<-function(pairs){
            library('magic')
            res<-list()
            for(i in 1:length(pairs)){
				pair<-pairs[[i]]
				start<-pair[[1]]
				end<-pair[[2]]
                res[i]<-list(solve(S_lambda[start:end,start:end]))
            }
            return (do.call('adiag',res))
        }       

        seq_pair<-list()
        tmp<-1
        for(i in 1:(length(csize))){
            seq_pair[[i]]<-list(tmp,csize[i]+tmp-1)
            tmp<-tmp+csize[i]
        }
        
        size<-length(seq_pair)/core_num
        node_seq<-list()
        tmp<-1
        for (i in 1:core_num){
            node_seq[i]<-list(seq_pair[tmp:(tmp+size-1)])
            tmp<-tmp+size
        }
        
        cl <- makeCluster(core_num)
		registerDoParallel(cl)
        inner_par<-list(S_lambda,node_seq)
		Omega <- foreach(each=1:core_num, .combine='adiag', .export='inner_par') %dopar% 
        {
            S_lambda<-inner_par[[1]]
            node_seq<-inner_par[[2]]
			each_inv(node_seq[[each]])
		}
        stopImplicitCluster()
    }
    else if(core_num==1){
        Omega<-zeros(p)
        tmp<-1
        for(i in 1:(length(csize))){
            Omega[tmp:(csize[i]+tmp-1),tmp:(csize[i]+tmp-1)]<-
                +solve(S_lambda[tmp:(csize[i]+tmp-1),tmp:(csize[i]+tmp-1)])
            tmp<-tmp+csize[i]
        }
    }
    return (list(soliution=Omega,var_seq=var_seq))
}