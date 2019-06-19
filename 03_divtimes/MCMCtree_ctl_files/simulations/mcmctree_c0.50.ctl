          seed = -1
       seqfile = alignment.txt
      treefile = tree.trees
	  mcmcfile = mcmc.txt
       outfile = out.txt

         ndata = 1
       usedata = 1    * 0: no data; 1:seq like; 2:use in.BV; 3: out.BV
         clock = 1    * 1: global clock; 2: independent rates; 3: correlated rates
       RootAge = B(0.8, 1.2)  * safe constraint on root age, used if no fossil for root.
       TipDate = 1 1

         model = 0    * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
         alpha = 0    * alpha for gamma rates at sites
         ncatG = 5    * No. categories in discrete gamma

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

       BDparas = 1 1 0 0.001   * birth, death, sampling
   kappa_gamma = 1 1      * gamma prior for kappa
   alpha_gamma = 1 1      * gamma prior for alpha

   rgene_gamma = 2 1   * gamma prior for overall rates for genes
  sigma2_gamma = 1 1    * gamma prior for sigma^2     (for clock=2 or 3)

      finetune = 1: .0935  0.023  0.3  0.015 .05  * times, rates, mixing, paras, RateParas

         print = 1
        burnin = 5000
      sampfreq = 10
       nsample = 20000

*** Note: Make your window wider (100 columns) before running the program.
