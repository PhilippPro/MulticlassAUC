

# Function sAUC.coeff
sAUC.coeff<-function(pred, truth, c1, c2, j=1) {
  pred.matrix=as.matrix(pred)
  n=nrow(pred.matrix)
  levels=levels(truth)
  temp=0
  for (u in c(1:n)) {
    for (v in c(1:n)) {
      if ((truth[u]==levels[c1])&&(truth[v]==levels[c2])) {
        if (pred.matrix[u,j]>pred.matrix[v,j]) {
          temp=temp+pred.matrix[u,j]-pred.matrix[v,j]
        } 
      }
    }
  }
  length.c1=length(which((truth==levels[c1])))
  length.c2=length(which((truth==levels[c2])))
  return(temp/(length.c1*length.c2))
}


binaryclass.Scoredauc.naive = function(pred, truth,  names=FALSE) {
  ptm <- proc.time()
  # Parameters
  K=dim(pred)[2]
  n=dim(pred)[1]
  levels=levels(truth)
  print(paste("levels",levels))
  
  # Making the matrix of probabilities p(i,j) is pred
  
  # function I
  I=function(a,b){
    res=0
    if (a>b){res=1}
    if (a==b){res=0.5}
    return(res)
  }
  
  
  # function makeMatrixI
  makeMatrixI.scored=function(j, pred){
    n=dim(pred)[1]
    matrix.I.temp=matrix(NA,nrow = n,ncol = n)
    for (i.i in c(1:n)) {
      for (i.t in c(1:n)) {
        matrix.I.temp[i.i,i.t]=I(pred[i.i,j],pred[i.t,j])*(pred[i.i,j]-pred[i.t,j])
      }
    }
    return(matrix.I.temp)
  }
  
  
  # Design Matrix f
  f= model.matrix( ~ 0 + truth, truth)
  f.summed=apply(f,2,sum)
  
  # matrix AUC
  AUC=matrix(rep(NA, K*K),K,K)
  
  #   for (j in c(1:K)) {
  #     ptmMakematrix <- proc.time()
  #     matrixI.tempj=makeMatrixI.scored(j,pred)
  #     print(proc.time()-ptmMakematrix)
  #     for (k in c(1:K)) {
  #       if (k!=j) {
  #         AUC[k,j]=t(f[,j])%*%matrixI.tempj%*%f[,k]/(f.summed[j]*f.summed[k])
  #       }
  #     }
  #   }
  
  for (j in c(1:K)) {
    for (k in c(1:K)) {
      if (k!=j) {
        temp=0
        for (u in c(1:n)) {
          for (v in c(1:n)) {
            if ((truth[u]==levels[j])&&(truth[v]==levels[k])) {
              if (pred[u,j]>pred[v,j]) {
                temp=temp+pred[u,j]-pred[v,j]
              } 
              else if (pred[u,j]==pred[v,j])
                temp=temp+(1/2)*(pred[u,j]-pred[v,j])
            }
          }
        }
        AUC[k,j]=temp/(f.summed[j]*f.summed[k])  
      }
    }
  }
  
  if (names==TRUE) {
    colnames(AUC)<-levels
    rownames(AUC)<-levels
  }
  print(proc.time()-ptm)
  return(AUC)
}



binaryclass.Scoredauc.naive.permutations = function(pred, truth,  names=FALSE) {
  ptm <- proc.time()
  # Parameters
  y<-truth
  X<-as.matrix(pred)
  ncol=ncol(X) #maj k
  n=nrow(X)
  levels.values=as.factor(levels(as.factor(truth)))
  levels.names=levels(as.factor(truth))
  K=length(levels.names)
  L = matrix(rep(levels.values, each = n), n, K)
  permutations = combs(1:K, 2)
  permutations = rbind(combs(1:K, 2),combs(K:1, 2))
  nP = nrow(permutations)
  Auc = matrix(0.5, nP, ncol) # Creates the AUC matrix

  
  # Making the matrix of probabilities p(i,j) is pred
  
  # function I
  I=function(a,b){
    res=0
    if (a>b){res=1}
    if (a==b){res=0.5}
    return(res)
  }
  
  
  # function makeMatrixI
  makeMatrixI.scored=function(j, pred){
    pred.matrix=as.matrix(pred)
    n=nrow(pred.matrix)
    matrix.I.temp=matrix(NA,nrow = n,ncol = n)
    for (i.i in c(1:n)) {
      for (i.t in c(1:n)) {
        matrix.I.temp[i.i,i.t]=I(pred.matrix[i.i,j],pred.matrix[i.t,j])*(pred.matrix[i.i,j]-pred.matrix[i.t,j])
      }
    }
    return(matrix.I.temp)
  }
  
  
  # Design Matrix f
  f= model.matrix( ~ 0 + truth, truth)
  f.summed=apply(f,2,sum)
  
  for (j in c(1:ncol)) {
    ptmMakematrix <- proc.time()
    matrixI.tempj=makeMatrixI.scored(j,pred)
    print(proc.time()-ptmMakematrix)
    
    for (i in 1:nP) {
      c1 = permutations[i, 1]
      c2 = permutations[i, 2]
      Auc[i,j]=t(f[,c1])%*%matrixI.tempj%*%f[,c2]/(f.summed[c1]*f.summed[c2])
    }
    
  }
  
  
  #Auc<-apply(Auc,c(1,2),function(x) return(max(x,1-x)))
  # Add the names
  if (names==TRUE) {
    rownames(Auc) = paste(levels.values[permutations[, 1]], " vs. ", levels.values[permutations[, 2]], sep = "")
    colnames(Auc) = colnames(X)
  }
  print(proc.time()-ptm)
  return(Auc)
}




## Scored AUC
binaryclass.Scoredauc = function(pred, truth,  names=FALSE) {
  ptm <- proc.time()
  # Parameters
  y<-truth
  X<-as.matrix(pred)
  ncol=ncol(X) #maj k
  n=nrow(X)
  levels.values=as.factor(levels(as.factor(truth)))
  levels.names=levels(as.factor(truth))
  K=length(levels.names)
  L = matrix(rep(levels.values, each = n), n, K)
  permutations = combs(1:K, 2)
  permutations = rbind(combs(1:K, 2),combs(K:1, 2))
  nP = nrow(permutations)
  Auc = matrix(0.5, nP, ncol) # Creates the AUC matrix
  count=0
  
  
  for (j in c(1:ncol)) {
    x = sort(X[, j], index = TRUE) # sorting with increasing probability for class j
    
    idx = y[x$ix] # indexes by increasing order of probabilities
    d = (matrix(rep(idx, K), n, K) == L) # Transforming the predicted vector in K columns of TRUE/FALSE
    d1=d
    d=apply(d,2,cumsum) # Cumulative sum of number in class
    dd = rbind(matrix(0, 1, K), d) # Add one line for integral computation
    nD = nrow(d) # New number of lines for dj
    d2=d1*1
    
    for (i in 1:nP) {
      c1 = permutations[i, 1]
      c2 = permutations[i, 2]
      number = d[nD, c1] * d[nD, c2]
      sum=0
      sumInt=0
      l<-NULL
      
      for (p in c(1:n)) {
        if (d2[p,c1]==1) {
          index=which(d2[c(1:p),c2]==1)
          if (length(index)>0) {
            for (q in c(1:length(index))) {
              indexq=index[q]
              diffProba = X[(x$ix[p]),j]-X[(x$ix[indexq]),j]
              #diffProba=1
              sum=sum + diffProba
              l<-c(l,diffProba)
              sumInt=sumInt+1
              count=count+1
            }
          }
        }
      }
      
      
      Auc[i, j] = sum/number # En gros le nombre de non variations de c2 par rapport a c1 qui est constant
    }
  }
  
  Auc<-abs(Auc)
  #Auc<-apply(Auc,c(1,2),function(x) return(max(x,1-x)))
  # Add the names
  if (names==TRUE) {
    rownames(Auc) = paste(levels.values[permutations[, 1]], " vs. ", levels.values[permutations[, 2]], sep = "")
    colnames(Auc) = colnames(X)
  }
  
  print(proc.time()-ptm)
  print(paste("count",count))
  return(Auc)
}






## Scored AUC
binaryclass.Scoredauc.final = function(pred, truth,  names=FALSE) {

  # Initializing
  y<-truth
  X<-as.matrix(pred)
  ncol=ncol(X)
  n=nrow(X)
  levels.values=as.factor(levels(as.factor(truth)))
  levels.names=levels(as.factor(truth))
  K=length(levels.names)
  L = matrix(rep(levels.values, each = n), n, K)
  permutations = rbind(combs(1:K, 2),combs(K:1, 2))
  nP = nrow(permutations)
  Auc = matrix(0.5, nP, ncol)
  
  
  for (j in c(1:ncol)) {                      # For all the classes probability predictions
    
    x = sort(X[, j], index = TRUE)              # sort with increasing probability for class j
    idx = y[x$ix]                               # get the resulting indexes
    d = (matrix(rep(idx, K), n, K) == L)        # Transform the predicted vector in K columns which have TRUE/FALSE values
    d1=d
    d=apply(d,2,cumsum)                         # Do the cumulative sum of element for each class
    nD = nrow(d)                                
    d2=d1*1                                     # Get the previously created matrix n*K with 1 and 0 instead of TRUE/FALSE
    
    for (i in 1:nP) {                           # For all the possible permutations
      c1 = permutations[i, 1]
      c2 = permutations[i, 2]
      number = d[nD, c1] * d[nD, c2]              # Number of elements in the cross product
      sum = 0
      
      for (p in c(1:n)) {                         # For all the observations
        if (d2[p,c1]==1) {                          # Get the 1 of the c1 column in the d2 matrix
          index=which(d2[c(1:p),c2]==1)             # count the number of elements in the c2 columns that are 1
          if (length(index)>0) {
            for (q in c(1:length(index))) {
              indexq=index[q]
              diffProba = X[(x$ix[p]),j]-X[(x$ix[indexq]),j]  # Compute the difference of probabilities for these points
              sum=sum + diffProba
            }
          }
        }
      }
      Auc[i, j] = sum/number 
    }
  }
  if (names==TRUE) {
    rownames(Auc) = paste(levels.values[permutations[, 1]], " vs. ", levels.values[permutations[, 2]], sep = "")
    colnames(Auc) = colnames(X)
  }
  return(Auc)
}




## Scored AUC
binaryclass.Scoredauc.final2 = function(pred, truth,  names=FALSE) {
  
  # Initializing
  y<-truth
  X<-as.matrix(pred)
  ncol=ncol(X)
  n=nrow(X)
  levels.values=as.factor(levels(as.factor(truth)))
  levels.names=levels(as.factor(truth))
  K=length(levels.names)
  L = matrix(rep(levels.values, each = n), n, K)
  permutations = rbind(combs(1:K, 2),combs(K:1, 2))
  nP = nrow(permutations)
  Auc = matrix(0.5, nP, ncol)
  
  
  for (j in c(1:ncol)) {                      # For all the classes probability predictions
    
    x = sort(X[, j], index = TRUE)              # sort with increasing probability for class j
    idx = y[x$ix]                               # get the resulting indexes
    d = (matrix(rep(idx, K), n, K) == L)        # Transform the predicted vector in K columns which have TRUE/FALSE values
    d2=d*1
    d=apply(d,2,cumsum)                         # Do the cumulative sum of element for each class
    nD = nrow(d)                                                                    # Get the previously created matrix n*K with 1 and 0 instead of TRUE/FALSE
    
    for (i in 1:nP) {                           # For all the possible permutations
      c1 = permutations[i, 1]
      c2 = permutations[i, 2]
      number = d[nD, c1] * d[nD, c2]              # Number of elements in the cross product
      sum = 0
      
      positiveProbability.number = d2[,c1]*d[,c2]
      negativeProbability.number = d2[,c2]*(d[nD, c1]-d[,c1])
      combined.number = positiveProbability.number-negativeProbability.number
      res = sum(x$x*combined.number)
      
      Auc[i, j] = res/number 
    }
  }
  if (names==TRUE) {
    rownames(Auc) = paste(levels.values[permutations[, 1]], " vs. ", levels.values[permutations[, 2]], sep = "")
    colnames(Auc) = colnames(X)
  }
  return(Auc)
}


