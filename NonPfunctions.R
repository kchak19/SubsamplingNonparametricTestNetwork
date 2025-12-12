library("foreach")
library("doParallel")
library("irlba")
library("distances")
library("matrixStats")
library("Rfast")
library("data.table")

# function to compute ASE of an adjacency matrix A
ASE = function(A, dim){
  if(nrow(A) >= 400){
    A.svd = irlba(A, nu = dim, nv = dim)
    A.svd.values = A.svd$d[1:dim]
    A.svd.vectors = A.svd$v[,1:dim]
    if(dim == 1)
      A.coords = sqrt(A.svd.values) * A.svd.vectors
    else
      A.coords = A.svd.vectors %*% diag(sqrt(A.svd.values))
  }
  else{
    A.svd = svd(A)
    if(dim == 1)
      A.coords = A.svd$v[,1] * sqrt(A.svd$d[1])
    else
      A.coords = A.svd$v[,1:dim] %*% diag(sqrt(A.svd$d[1:dim]))
  }
  
  return(A.coords)
}

#procrustes transformation between matrices X and Y
procr= function(X,Y){
  tmp = t(X) %*% Y
  tmp.svd = svd(tmp)
  W = tmp.svd$u %*% t(tmp.svd$v)
  return(list(d = norm(X%*%W - Y, type = "F"), W = W, X.tfrm = X%*%W))
}

#generate a graph from RDPG(X)
graph = function(X){
  n = nrow(X)
  # P = X %*% t(X)
  P = Rfast::Tcrossprod(X,X)
  A = matrix(0, nrow = n, ncol = n)
  A[col(A) > row(A)] = runif(n*(n-1)/2)
  A = (A + t(A))
  A = (A < P) + 0 ;
  diag(A) = 0
  return(A)
}

#subsampled estimation of ASE
estSS = function(A,o,k){
  n = nrow(A)
  b = (n-o)/k
  t1 = proc.time()[3]
  x = 1:n
  common = sample(x, size = o, replace = FALSE)
  x = setdiff(x,common)
  samp = c()
  for(i in 1:k)
  {
    samp = rbind(samp, sample(x, size = b, replace = FALSE))
    x = setdiff(x,samp)
  }
  # samp has all the partition indices
  t2 = proc.time()[3]
  
  ind1 = c(common,samp[1,])
  A.ref = A[ind1,ind1]
  X.hat.ref = ASE(A.ref,d)
  X.ref.0 = X.hat.ref[1:o,]
  X.ref.1 = X.hat.ref[(o+1):(o+b),]
  # Finding the reference latent vectors
  t3 = proc.time()[3]
  
  
  X.part = foreach(i = 2:k, .combine = rbind) %dopar%
    {
      A.sub = A[c(common,samp[i,]),c(common,samp[i,])]
      X.hat.sub = ASE(A.sub,d)
      X.sub.0 = X.hat.sub[1:o,]
      X.sub.i = X.hat.sub[(o+1):(o+b),]
      
      H = procr(X.sub.0, X.ref.0)$W
      X.sub.i.trans = X.sub.i %*% H
      X.sub.i.trans
    }
  
  t4 = proc.time()[3]
  
  # rearrange the vertices into the final matrix
  X.fin = matrix(0, n, d)
  
  for(j in 1:length(common)) {X.fin[common[j],] = X.ref.0[j,]}
  for(j in 1:length(samp[1,])) {X.fin[samp[1,j],] = X.ref.1[j,]}
  for(i in 2:k) for(j in 1:length(samp[i,])) {X.fin[samp[i,j],] = X.part[b*(i-2)+j,]}
  t5 = proc.time()[3]
  
  Times = data.frame(subsampling = t2-t1, refEst = t3-t2, partitionEst = t4-t3,
                     rearrange = t5-t4)
  
  return(X.fin)
  #return(list("Estimate"= X.fin, "Times" = Times))
}

#general estimation function with or without subsampling (COASub/Full Network)
est = function(A,o,k, subs = "SS"){
  if(subs == "ASE")
    return(ASE(A,d))
  if(subs == "SS")
    return(estSS(A,o,k))
  gc()
}

#rectangular distance between X and Y
rect.dist.fast = function(X,Y){
  n = nrow(X)
  m = nrow(Y)
  D = Rfast::Outer(rep(1, n), rowSums(Y^2)) - Rfast::Tcrossprod(Y,2*X) + 
    Rfast::Outer(rowSums(X^2), rep(1,m))
  return(D)
}

#median-matching
find.transform = function(X,Y){
  u = apply(X,2,median)
  v = apply(Y,2,median)
  if(ncol(X) == 1){
    T = sign(u/v)
  }
  if(ncol(X) > 1){
    T = sign(u/v)
  }
  nan = which(is.nan(T))
  n0 = which(T==0)
  T[nan] = rep(1, length(nan))
  T[n0] = sign(v[n0])
  return(diag(T))
}

#calculate the U-statistic between X and Y, using particular rotation method
UstatPar = function(X,Y,sig=0.5, method = "medmatch"){
  n = nrow(X)
  m = nrow(Y)
  if(method == "medmatch")
    TPP = find.transform(X, Y)
  if(method == "procr")
    TPP = procr(X, Y)$W
  if(method == "none")
    TPP = diag(1,d)
  
  X1 = X %*% TPP
  #KXX = exp(-(as.matrix(parDist(X1, threads = 6))^2)/(2*sig^2))
  #KXX = exp(-(as.matrix(dist(X1))^2)/(2*sig^2))
  
  # KXX = exp(-(as.matrix(dist(X1))^2)/(2*sig^2))
  # KYY = exp(-(as.matrix(dist(Y))^2)/(2*sig^2))
  
  KXX = exp(-(distance_matrix(distances::distances(X1)))^2/(2*sig^2))
  KYY = exp(-(distance_matrix(distances::distances(Y)))^2/(2*sig^2))
  KXY = exp(-(rect.dist.fast(X1,Y))/(2*sig^2))
  
  Usum = function(n, K) {n + 2*sum2(K)}
  stat = (m+n)*(Usum(n, KXX)/(n*(n-1)) + Usum(m, KYY)/(m*(m-1)) - 2*sum2(KXY)/(m*n))
  
  #stat = (m+n)*(sum2(KXX)/(n*(n-1)) + sum2(KYY)/(m*(m-1)) - 2*sum2(KXY)/(m*n))
  rm(KXX)
  rm(KYY)
  rm(KXY)
  return(stat)
}

#perform the test of equality between A and B with or without subsampling (COASub/Full Network)
sec.test = function(A, B, boot, ref=3, method = "medmatch", subs = "SS"){
  n = nrow(A)
  m = nrow(B)
  
  X = est(A,on,k,subs)
  Y = est(B,om,k,subs)
  
  stat = UstatPar(X, Y, sig, method)
  
  if(method == "medmatch")
    TPP = find.transform(X, Y)
  if(method == "procr")
    TPP = procr(X, Y)$W
  if(method == "none")
    TPP = diag(1,d)
  Xnew = X %*% TPP
  
  rowX = as.list(data.frame(t(Xnew)))
  rowY = as.list(data.frame(t(Y)))
  
  if(ref==1) FhatRows = rowX
  if(ref==2) FhatRows = rowY
  if(ref==3) FhatRows = c(rowX,rowY)
  
  Star = c()
  
  for(bt in 1:boot){
    #cat(bt)
    Xs = sample(FhatRows, nrow(A), replace = T)
    Xstar = t(cbind.data.frame(Xs))
    Ys = sample(FhatRows, nrow(B), replace = T)
    Ystar = t(cbind.data.frame(Ys))
    
    Astar = graph(Xstar)
    Bstar = graph(Ystar)

    Xhatstar = est(Astar, on, k, subs)
    Yhatstar = est(Bstar, om, k, subs)
    Star[bt] = UstatPar(Xhatstar, Yhatstar, sig, method)
    
    #Star[bt] = UstatPar(Xstar, Ystar, sig, method)
  }
  ls = list(main = stat, stats = Star, 
            pval = sum2(stat <= Star)/boot)
  return(ls)
}

#perform the test of equality between A and B by partition subsampling (PartSub)
test.split.pval = function(A, B, boot, k, ref=3, method = "medmatch"){
  n = nrow(A)
  m = nrow(B)
  bn = n/k
  bm = m/k

  x = 1:n
  sampn = c()
  for(i in 1:k){
    sampn = rbind(sampn, (1:bn)+bn*(i-1))
    x = setdiff(x,sampn)
  }
  y = 1:m
  sampm = c()
  for(i in 1:k){
    sampm = rbind(sampm, (1:bm)+bm*(i-1))
    y = setdiff(y,sampm)
  }

  ls = foreach(i = 1:k, .combine = c) %dopar%
  {
    VA = sampn[i,]
    VB = sampm[i,]
    A0 = A[VA,VA]
    B0 = B[VB,VB]
    output = sec.test(A0, B0, boot, ref, method, subs = "ASE")
    output$pval
  }

  S = min(ls)
  out = list(pval = round(S,3), ls = ls)
  return(out)
}

#perform the test of equality between A and B by pure subsampling (PureSub)
test.partial.pval = function(A, B, boot, rho, repl, ref=3, method = "medmatch"){
  n0 = rho(nrow(A))
  m0 = rho(nrow(B))

  ls = foreach(i = 1:repl, .combine = c) %dopar%
  {
    VA = sample.int(nrow(A), n0)
    VB = sample.int(nrow(B), m0)

    A0 = A[VA, VA]
    B0 = B[VB, VB]
    output = sec.test(A0, B0, boot, ref, method, subs = "ASE")
    output$pval
  }
  
  S = min(ls)
  out = list(pval = round(S,3), ls = ls)
  return(out)
}

#perform the test of scaling between A and B with or without subsampling (COASub/Full Network)
scale.test = function(A, B, boot, subs = "SS", c = "default", 
                      method = "medmatch"){
  n = nrow(A)
  m = nrow(B)
  
  X = est(A,on,k,subs)
  Y = est(B,om,k,subs)
   
  sX = norm(X, type = "F")/sqrt(n)
  sY = norm(Y, type = "F")/sqrt(m)
  
  stat = UstatPar(X/sX, Y/sY, sig, method)
  
  # Z = Y * sX/sY
  # TPP = find.transform(X, Z)
  # Xnew = X %*% TPP
  # 
  
  rowX = as.list(data.frame(t(X)))
  #rowX = as.list(data.frame(t(Xnew)))
  #rowY = as.list(data.frame(t(Z)))
  #Fhat = c(rowX, rowY)
  if(c == "default") c0 = sX/sY
  if (c != "default") c0 = c
  
  Star = c()
  
  for(bt in 1:boot){
    Xs = sample(rowX, nrow(A), replace = T)
    Xstar = t(cbind.data.frame(Xs))
    Ys = sample(rowX, nrow(B), replace = T)
    Ystar = t(cbind.data.frame(Ys))/c0
    
    Astar = graph(Xstar)
    Bstar = graph(Ystar)

    Xhatstar = est(Astar, on, k, subs)
    Yhatstar = est(Bstar, on, k, subs)

    sXstar = norm(Xhatstar, type = "F")/sqrt(n)
    sYstar = norm(Yhatstar, type = "F")/sqrt(m)
    Star[bt] = UstatPar(Xhatstar/sXstar, Yhatstar/sYstar, sig, method)
    
    # sXstar = norm(Xstar, type = "F")/sqrt(n)
    # sYstar = norm(Ystar, type = "F")/sqrt(m)
    # Star[bt] = UstatPar(Xstar/sXstar, Ystar/sYstar, sig, method)
  }
  #cat("\n")
  #ls = list(main = stat, stats = Star, pval = sum2(stat <= Star)/boot)
  ls = list(c = sX/sY, pval = sum2(stat <= Star)/boot)
  return(ls)
}

#perform the test of scaling between A and B by partition subsampling (PartSub)
scale.split.pval = function(A, B, boot, k, method = "medmatch"){
  n = nrow(A)
  m = nrow(B)
  bn = floor(n/k)
  bm = floor(m/k)
  
  x = 1:n
  sampn = c()
  for(i in 1:k){
    sampn = rbind(sampn, sample(x,bn))
    x = setdiff(x,sampn)
  }
  y = 1:m
  sampm = c()
  for(i in 1:k){
    sampm = rbind(sampm, sample(y,bm,replace = F))
    y = setdiff(y,sampm)
  }
  
  VA = sampn[1,]
  VB = sampm[1,]
  A0 = A[VA,VA]
  B0 = B[VB,VB]
  output0 = scale.test(A0, B0, boot, subs = "ASE")
  p0 = output0$pval
  c0 = output0$c
  
  ls = foreach(i = 2:k, .combine = c) %dopar%
    {
      VA = sampn[i,]
      VB = sampm[i,]
      A0 = A[VA,VA]
      B0 = B[VB,VB]
      output = scale.test(A0, B0, boot, subs = "ASE", c=c0)
      output$pval
    }
  p = min(c(p0, ls))
  out = list(pval = round(p,3), ls = c(p0, ls))
  return(out)
}

#perform the test of scaling between A and B by pure subsampling (PureSub)
scale.partial.pval = function(A, B, boot, rho, repl, method = "medmatch"){
  n0 = rho(nrow(A))
  m0 = rho(nrow(B))
  
  VA = sample.int(nrow(A), n0)
  VB = sample.int(nrow(B), m0)
  
  A0 = A[VA, VA]
  B0 = B[VB, VB]
  output0 = scale.test(A0, B0, boot, subs = "ASE")
  p0 = output0$pval
  c0 = output0$c
  
  ls = foreach(i = 1:repl, .combine = c) %dopar%
    {
      VA = sample.int(nrow(A), n0)
      VB = sample.int(nrow(B), m0)
      
      A0 = A[VA, VA]
      B0 = B[VB, VB]
      output = scale.test(A0, B0, boot, subs = "ASE", c=c0)
      output$pval
    }
  pval = min(c(p0, ls))
  out = list(pval = round(pval,3), ls = c(p0, ls))
  return(out)
}




