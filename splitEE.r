# The implementation of split EE.
# Pay attention that we haven't consider about storage cost of these functions. 
#-------------------------------------------------------------------------------
library("igraph")
library('magic')
library("matlab")
library('foreach')
library('doParallel')
source('thresholding.r')


splitAndInv<-function(blkMat){
    #======= TODO: The dim of matrix might be too samll and there might be more partitions =============
    dim<-length(blkMat[1,])
    size<-dim/2
    result<-zeros(dim)
    offset<-1
    partList<-c(blkMat[1:size,1:size],Matrix(blkMat[1:size,size+1:dim],sparse=TRUE),blkMat[size+1:dim,size+1:dim])
    result[1:size,1:size]<-solve(partList[1]-partList[2]%*%solve(partList[3])%*%t(partList[2]))
    result[size+1:dim,size+1:dim]<-solve(partList[3]-t(partList[2])%*%solve(partList[1])%*%partList[2])
    result[1:size,size+1:dim]<-(-solve(partList[3])%*%t(partList[2])%*%result[1:size,1:size]
    result[size+1:dim,size+1:dim]<-(-solve(partList[1])%*%partList[2]%*%result[size+1:dim,size+1:dim]
    return (result)
}

#Improved Elementary Wstimator
splitEE<-function(S,lambda,p,thr_func=hardThreshold,core_num=1){
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
    print(sprintf('%d components',num))
    #Permuting the matrix so that it is block diagonal
    var_seq<-zeros(p,2)
    var_seq[,1]<-1:p
    var_seq[,2]<-bins
    var_seq<-var_seq[order(var_seq[,2]),]
    var_seq<-var_seq[,1]
    S_lambda<-S_lambda[var_seq,]
    S_lambda<-S_lambda[,var_seq]
    #=============== TODO: Delte this part as we don't need parallel in our function===============
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
                +splitAndInv(S_lambda[tmp:(csize[i]+tmp-1),tmp:(csize[i]+tmp-1)])
            tmp<-tmp+csize[i]
        }
    }
    return (list(soliution=Omega,var_seq=var_seq))
}