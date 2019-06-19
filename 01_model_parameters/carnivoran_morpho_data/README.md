# General description of the CONTML directory

Inside the **CONTML** directory you will find the analyses carried out when accounting for or ignoring within-lineage character correlation 
in the morphological carnivoran data sets. For each analysis, you can find the input and output files used by the `CONTML` software.
Click [here](https://github.com/sabifo4/morpho/blob/master/01_model_parameters/carnivoran_morpho_data/CONTML) to access the **CONTML** directory,
where more details about these files and how the software was run are given.

# Choosing the gamma prior on rates 

* **Define parameters**:  
   Branch length: $b\sim 20$ (approximated morphological distance from the tips to the root)  
   Times for the Bayes factors and divergence times estimation (mean of gamma prior): $t_{div}=52$   
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
   Scale parameter: $\beta = ?$    
    - For divergence times estimation:  
	   - $\beta=\alpha\times t_{div}/b=2\times 52/20=5.2\sim 5$  
	   - Gamma prior: $\Gamma(2,5)$  
	   - The control file for `MCMCTree` can be found [here](https://github.com/sabifo4/morpho/tree/master/03_divtimes/MCMCtree_ctl_files/carnivoran_morpho_mol_data). Note that we only carried out the Bayesian inference with specifying the clock model to follow and independent log-normal distribution as it was the clock model supported by the BFs analysis for both the morphological and molecular partitions.
