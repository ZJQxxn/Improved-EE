library('matlab')
library('magic')
library('MASS')

sparseMatrix<-function(n,coef=1){
    mat=eye(n)
    for (i in 1:n){
        j<-i+1
        while (j<=n){
            if(abs(i-j)%%coef==0){
                tmp<-0.7^abs(i-j)
                mat[i,j]<-tmp
                mat[j,i]<-tmp
            }
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
    mat<-eye(n)
    for(i in 1:n){
        mat[i,(i+1)%%n]<-0.3
        mat[(i+1)%%n,i]<-0.3
    }
    return (mat)
}

gridMatrix<-function(row,col){
    # row and col must be larger than 1
    mat<-eye(row*col)
    for(i in 1:(row-1)){
        for(j in 1:col){
            index=(i-1)*col+j
            mat[index,index+1]<-0.1
            mat[index+1,index]<-0.1
            mat[index,index+col]<-0.1
            mat[index+col,index]<-0.1
        }
    }
    return (mat)
}

balancedBlk<-function(p,blk=4,coef=1){
    n<-p/blk
    #mat<-randMatrix(n)
    blk_list<-list()
    for(i in 1:blk){
        blk_list[[i]]<-sparseMatrix(n,coef=coef)
    }
    mat<-do.call('adiag',blk_list)
    return (mat)
}

inbalancedBlk<-function(a,b,c,d){
    mat1<-sparseMatrix(a)
    mat2<-sparseMatrix(b)
    mat3<-sparseMatrix(c)
    mat4<-sparseMatrix(d)
    return (adiag(mat1,mat2,mat3,mat4))
}

circleGraph<-function(p,blk=4){
    n<-p/blk
    blk_list<-list()
    for(i in 1:blk){
        blk_list[[i]]<-circleMatrix(n)
    }
    mat<-do.call('adiag',blk_list)
    return (mat)
}

gridGraph<-function(p,blk,row,col){
    n<-p/blk
    blk_list<-list()
    for(i in 1:blk){
        blk_list[[i]]<-gridMatrix(row,col)
    }
    mat<-do.call('adiag',blk_list)
    return (mat)
}

generateSample<-function(Omega,size){
    samples<-mvrnorm(size,ones(1,length(Omega[1,])),solve(Omega))
    return (samples)
}