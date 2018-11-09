sparseMatrix<-function(n){
    library('matlab')
    mat=eye(n)
    for (i in 1:n){
        j<-i+1
        while (j<=n){
            tmp<-0.7^abs(i-j)
            mat[i,j]<-tmp
            mat[j,i]<-tmp
            j<-j+1
        }
    }
    noise<-eye(n)
    noise<-generateSample(noise,n)
    noise<-noise/abs(max(abs(noise)))
    noise<-noise*0.001
    mat<-mat+noise  
    return (mat)
}

circleMatrix<-function(n){
    
}

balancedBlk<-function(p,blk=4){
    library('magic')
    library('MASS')
    n<-p/blk
    #mat<-randMatrix(n)
    blk_list<-list()
    for(i in 1:blk){
        blk_list[[i]]<-sparseMatrix(n)
    }
    mat<-do.call('adiag',blk_list)
    return (mat)
}

inbalancedBlk<-function(a,b,c,d){
    library('magic')
    mat1<-sparseMatrix(a)
    mat2<-sparseMatrix(b)
    mat3<-sparseMatrix(c)
    mat4<-sparseMatrix(d)
    return (adiag(mat1,mat2,mat3,mat4))
}

circleGraph<-function(p,blk=4){

}

generateSample<-function(Omega,size){
    library('matlab')
    library('MASS')
    samples<-mvrnorm(size,ones(1,length(Omega[1,])),solve(Omega))
    return (samples)
}