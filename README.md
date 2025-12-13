# SubsamplingNonparametricTestNetwork

This is the README file for the R codes used in Subsampling Based Efficient Estimation and Nonparametric Two-Sample Testing for Large Networks by Chakraborty et. al.. The functions are available in NonPfunctions.R, and the simulation experiments can be recreated from NonPFinSim.R.

**Description of main functions:**

1) ASE(A,dim) : function to compute ASE of an adjacency matrix A at dimension dim (full network)
2) procr(X,Y) : procrustes transformation between matrices X and Y
3) graph(X) : generate a graph from RDPG(X)
4) estSS(A,o,k) : subsampled estimation of ASE of A with overlap size o and number of subsamples k (COASub)
5) est : estimation function that combines ASE and estSS, specify with subs
6) UstatPar(X,Y,sig,method) : compute the U-statistic between X and Y using kernel bandwidth sig and a rotation method for Y
7) sec.test(A,B,boot,subs,...) : test of equality between A and B with boot number of resamples, specify full network/COASub with subs
8) test.split.pval(A,B,boot,k,...) : test of equality between A and B with boot number of resamples, using PartSub (k partitions)
9) test.partial.pval(A,B,boot,rho,repl,...) : test of equality between A and B with boot number of resamples, using PureSub (size rho, repl number of replicates)
10) scale.test(A,B,boot,subs,c,...) : test of scaling between A and B with boot number of resamples, specify full network/COASub with subs. c is unknown by default, but option to use it as known.
11) scale.split.pval(A,B,boot,k,...) : test of scaling between A and B with boot number of resamples, using PartSub (k partitions)
12) scale.partial.pval(A,B,boot,rho,repl,...) : test of scaling between A and B with boot number of resamples, using PureSub (size rho, repl number of replicates)
