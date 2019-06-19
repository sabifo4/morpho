		seed   = -1
       seqfile = Ms_19sp_carnivores.aln
      treefile = tree.trees
       outfile = out.txt
      mcmcfile = mcmc.txt
	  
       seqtype = 0  * 0: nucleotides; 1:codons; 2:AAs;
         noisy = 3
       usedata = 1 * 0: no data; 1:seq like; 2:use in.BV; 3: out.BV

         ndata = 1   * 
         clock = 3   * 1: global clock; 2: independent rates; 3: correlated rates
       TipDate = 1 1  * TipDate (1) & time unit

        RootAge = B(37.3, 66.0, 0.025, 0.025)  * used if no fossil for root

         model = 0    * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

       BDparas = 1 1 0 0.001  * lambda, mu, rho, psi for birth-death-sampling model

   kappa_gamma = 2 1   * gamma prior for kappa
   alpha_gamma = 1 1   * gamma prior for alpha
   rgene_gamma = 2 5   * gamma prior for rate for genes
  sigma2_gamma = 2 2  * gamma prior for sigma^2  (for clock=2)

      finetune = 1: 0.1  0.1  0.1  0.01 .5  * auto (0 or 1) : times, musigma2, rates, mixing, paras, FossilErr

         print = 2
        burnin = 50000
      sampfreq = 50
       nsample = 20000
	   