# Bayesian divergence times estimation

Once the best-fitting clock model was found in the Bayesian model selection for each data set, we inferred the divergence times using a Bayesian approach. 

**NOTE**: The path to the alignment and tree files in the control files are updated by a pipeline that runs `MCMCtree` in the [Apocrita HPC](https://docs.hpc.qmul.ac.uk/).

## Simulated data 

The control files used to run `MCMCtree` can be found [here](https://github.com/sabifo4/morpho/tree/master/03_divtimes/MCMCtree_ctl_files/simulations). 

The results obtained for the $R=1,000$ replicates simulated for each scenario were summarised as we detail below. Note that $i$ corresponds to node $i$ and $j$ to the replicate number:   
   * Mean posterior times, $\tilde{t}_{i,j}$ 
   * Mean 95% credibility intervals (CIs)   
   * Mean CI-width w (and relative CI-width $w_{r}=w/t_{i}$)   
   * Coverage (number of times the true node age falls within the 95% CI)   

Then, assuming that $i$ corresponds to the node number and $j$ to the replicate number, we calculated the bias and error as it follows:   
   
   * Mean bias, $b=\sum_{j=1}^{R}(\tilde{t}\_{i,j}\-t\_{i})/R$, and relative bias, $b_{r}=b/t_{i}$  
   * Mean squared error (MSE), $\epsilon=\sum_{j=1}^{R}(\tilde{t}\_{i,j}\-t_{i})^{2}/R$, and relative error, $\varepsilon_{r}=\varepsilon/t_{i}$ 

The directory with the scripts used to summarise the results obtained are [here](https://github.com/sabifo4/morpho/tree/master/03_divtimes/summary_simulations). Note that the working directory needs to be 
modified for each data set analysed.

## Carnivoran data 

### Molecule-only data set

The control files used to run `MCMCtree` for each partitioning scheme can be found [here](https://github.com/sabifo4/morpho/tree/master/03_divtimes/MCMCtree_ctl_files/carnivoran_mol_data).

### Morphology-only data set

The control files used to run `MCMCtree` accounting and ignoring the within-lineage character correlation can be found [here](https://github.com/sabifo4/morpho/tree/master/03_divtimes/MCMCtree_ctl_files/carnivoran_morpho_data).

### Joint data set with morphological data set

The control file used to run `MCMCtree` can be found [here](https://github.com/sabifo4/morpho/tree/master/03_divtimes/MCMCtree_ctl_files/carnivoran_morpho_mol_data).


