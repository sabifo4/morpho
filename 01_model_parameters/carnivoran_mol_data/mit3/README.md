# RAxML analysis

The command used in `RAxML` to estimate the tree topology for the alignment consisting of one partition with the third codon positions is the following:

```
raxmlHPC -f a -m GTRGAMMA -p 12345 -# 100 -x 12345 -# 500 -s mit_3.aln -o Nan_bin,Ail_ful,Urs_ame,Par_her -n mt3g_GCP
```

Click [here](https://github.com/sabifo4/morpho/blob/master/01_model_parameters/carnivoran_mol_data/mit3/analysis_RAxML) to access the output files by `RAxML` for this analysis.

# Estimated tree topology 

The tree topology estimated by RAxML can be found [here](https://github.com/sabifo4/morpho/blob/master/01_model_parameters/carnivoran_mol_data/mit3/analysis_RAxML/RAxML_bestTree.mt3g_GCP.jpg)
and the newick file [here](https://github.com/sabifo4/morpho/blob/master/01_model_parameters/carnivoran_mol_data/mit3/analysis_RAxML/RAxML_bestTree.mt3g_GCP).

# Choosing the gamma prior

According to the estimated tree topology, we chose our gamma prior for the rate as it follows:  

* **Define parameters**:  
   Branch length: $b\sim 1$  
   Times for divergence times estimation (mean of gamma prior): $t_{div}=52$   
   Branch length ($b$) in terms of rate ($r$) and times ($t$): $b=r\times t$   
   Evolutionary mean rate as the mean of the gamma prior ($\mu_{\Gamma}=\alpha/\beta$). We can now find $\beta$:   
   
     - **a**. We know that $b=r\times t$.   
     - **b**. We know that the mean of the gamma prior is $\mu_{\Gamma}=\alpha/\beta$.
     - **c**. We want to equate the mean rate as the mean of the gamma prior, that is $r=\alpha/\beta$.   
     - **d**. We can isolate the rate $r$ in equations defined in point `a` and `c`:   
		
			r = b / t    
			r = alpha / beta
			=================			
			b / t = alpha / beta
	
     - **e**. We can now isolate beta, our unknown parameter: $\beta=b/t\times \alpha$.



* **Define the gamma prior**:  
   Shape parameter: $\alpha=2$  
   Scale parameter: $\beta=?$    
	- **For Bayes factors (BFs)**:  
	    - $\beta=\alpha\times t_{BF}/b=2\times 1/1=2$  
	    - Gamma prior: $\Gamma(2,2)$  
	    - The control files for `mcmc3r` \- which will create $n=16$ control files readable by `MCMCTree` per clock model with the corresponding $\beta$ values to later estimate the marginal likelihood for each model and then compute the BFs \- can be found [here](https://github.com/sabifo4/morpho/tree/master/02_bayes.model.sel/MCMCtree_ctl_files/carnivoran_mol_data/mit3).  
    - **For divergence times estimation**:  
	   - $\beta=\alpha\times t_{div}/b=2\times 52/1=104\sim 100$  
	   - Gamma prior: $\Gamma(2,100)$  
	   - The control file for `MCMCTree` can be found [here](https://github.com/sabifo4/morpho/tree/master/03_divtimes/MCMCtree_ctl_files/carnivoran_mol_data/mit3). Note that we only carried out the Bayesian inference with specifying the clock model to follow geomtric Brownian motion (autocorrelated-rates model) as it was the clock model supported by the BFs analysis.
