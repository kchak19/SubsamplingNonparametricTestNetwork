#setwd("~/Desktop/Nonparametric")
source("NonPfunctions.R")

#initialize simulation experiment parameters
n = 10000 #first network size
m = 8000 #second network size
on = 1000 #first network overlap region size
om = 1000 #second network overlap region size
k = 10 #number of subsamples
d = 2 #rank of RDPG/number of communities in SBM
c = 1.5
boot = 200 #number of bootstrap resamples
sig = 0.4 #kernel bandwidth
iter = 100 #number of iterations to find rejection rate

#create the SBM for experiment
J = matrix(0.5, d, d)
L = diag(0.1,d)
Ba = L + J
dim(Ba) = c(d,d)
I = diag(rep(1,d))
eps = 0
Be = Ba + eps*I
pi = sapply(1:(d-1), function(i) {1/d})
#pi = sapply(1:(d-1), function(i) {0.49^i})
pi = c(pi, 1-sum(pi))

#prepare to generate from corresponding RDPG
E = eigen(Ba)
Qx = E$vectors %*% diag(sqrt(E$values))
Ee = eigen(Be)
Qxe = Ee$vectors %*% diag(sqrt(Ee$values))

rowQx = as.list(data.frame(t(Qx)))
rowQxe = as.list(data.frame(t(Qxe)))

#######################

registerDoParallel(6) #use suitable number of parallel cores to improve computation

#this simulation experiment can perform the methods full network/COASub/PureSub/PartSub
#for test of equality/scaling. Choose which method to be used and comment the rest
subs = "SS"
rho = function(x) {100*floor(x^0.82/100)} #size for PureSub

pval1 = pval2 = rep(NA, iter)
fin.output = c()
#perform the test repeatedly to find rejection rate
for(i in 1:iter){
  
  Sx = sample(rowQx, n, prob = pi, replace = T)
  X = t(cbind.data.frame(Sx))
  Sy = sample(rowQxe, m, prob = pi, replace = T)
  Y = t(cbind.data.frame(Sy))/c
  
  A = graph(X)
  B = graph(Y)
  T1 = proc.time()[3]
  #output1 = sec.test(A, B, boot, 3, "medmatch", subs)
  #output1 = test.split.pval(A, B, boot, 5, 3, "medmatch")
  #output1 = test.partial.pval(A, B, boot, rho, 5, ref=3, method = "medmatch")
  #T12 = proc.time()[3]
  
  output1 = scale.test(A, B, boot, subs)
  #output1 = scale.partial.pval(A, B, boot, rho, 5)
  #output1 = scale.split.pval(A, B, boot, 5)
  
  T2 = proc.time()[3]
  
  fin.output = c(fin.output, list(output1))
  pval1[i] = output1$pval
  #pval1[i] = output1
  print(paste0("Iter ", i, " time ", round(T2-T1,2), " pval= ", pval1[i]))
  gc()
}
stopImplicitCluster()

#output the rejection rates here

#pvals = sapply(1:iter, function(i) {fin.output[i]})

#pvals = sapply(1:iter, function(i) {fin.output[[i]]$pval})
# mains = sapply(1:iter, function(i) {fin.output[[i]]$main})
#ls = t(cbind.data.frame(lapply(1:iter, function(i) {fin.output[[i]]$ls})))
#row.names(ls) = 1:iter

c(sum(pval1<0.05)/iter , T2-T1)
#c(sum(pval1<0.223)/iter , T2-T1)
#c(sum(pval1<0.05/5)/iter , T2-T1)

# sink(paste0("Output",subs,".txt"))
# c(sum(pval<0.05)/iter , T2-T1)
# sink()
# saveRDS(fin.output, "ptl5n5k10k.rds")

