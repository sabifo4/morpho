# Bayesian model selection 

For the Bayesian clock model selection, we used the stepping-stone approach ([Xie et al. 2011](https://academic.oup.com/sysbio/article/60/2/150/2461669)) implemented in the `mcmc3r` R package, which prepares the 
directories with the control files needed for `MCMCtree` and uses its output to estimate the marginal likelihood and compute the Bayes factors and posterior probabilities
so the user can decide which model best fits the data. 

In these analyses, we used $n=16$ values of $\beta$, ranging from $\beta=0$ to $\beta=1$, so we could 
sample likelihood values from each corresponding power posterior. As the power posterior is proportional to 
the prior times the likelihood to the power of $\beta$, this procedure allows us to collect likelihood values 
in a path going from the prior (when $\beta=0$ we sampled from the prior
as $posterior\propto prior\times likelihood^{\beta}\propto prior\times likelihood^{0}\propto prior$)
to the posterior (when $\beta=1$ we sampled from the posterior
as $posterior\propto prior\times likelihood^{\beta}\propto prior\times likelihood^{1}\propto prior\times likelihood$). Once the MCMC chains in `MCMCtree` are completed, `mcmc3r` R package uses the output file with the sampled likelihood values to estimate the marginal likelihood values for each model, which later uses to compute the Bayes factors. 

For more information about the R commands needed to use the `mcmc3r` R package to estimate the marginal likelihood with 
`MCMCtree` output files, please follow the tutorial [here](https://dosreislab.github.io/2017/10/24/marginal-likelihood-mcmc3r.html). 

## Molecular data 

For the three molecular data sets with different partitioning schemes, we 
evaluated three different models regarding the clock:   

   * Strict clock (STR) model   
   * Geometric Brownian diffusion (GBM, also known as autocorrelated-rates) model   
   * Independent log-normal rate (ILN) model

The control files used by `mcmc3r` to generate the $n=16$ control files with the corresponding $\beta$ values
to later be used by `MCMCtree` can be found [here](https://github.com/sabifo4/morpho/tree/master/02_bayes.model.sel/MCMCtree_ctl_files/carnivoran_mol_data). 
Note that the path to the alignment and tree files in the control files are updated by a pipeline that runs `MCMCtree` in the [Apocrita HPC](https://docs.hpc.qmul.ac.uk/).

## Morphological data 
For the two morphological data sets, we evaluated six different models: 

   * Strict clock (STR) model + ignoring within-lineage character correlation ($\mathbf{R}=\mathbf{I}$)   
   * Strict clock (STR) model + accounting for within-lineage character correlation ($\mathbf{R}=\mathbf{R^{\*}}$)   
   * Geometric Brownian diffusion (GBM, also known as autocorrelated-rates) + ignoring within-lineage character correlation ($\mathbf{R}=\mathbf{I}$)   
   * Geometric Brownian diffusion (GBM, also known as autocorrelated-rates) + accounting for within-lineage character correlation ($\mathbf{R}=\mathbf{R^{\*}}$)   
   * Independent log-normal rate (ILN) models + ignoring within-lineage character correlation ($\mathbf{R}=\mathbf{I}$)
   * Independent log-normal rate (ILN) models + accounting for within-lineage character correlation ($\mathbf{R}=\mathbf{R^{\*}}$)

The control files used by `mcmc3r` to generate the $n=16$ control files with the corresponding $\beta$ values
to later be used by `MCMCtree` can be found [here](https://github.com/sabifo4/morpho/tree/master/02_bayes.model.sel/MCMCtree_ctl_files/carnivoran_morpho_data).


# Best-fitting models 

According to the posterior probabilities calculated for each model evaluated, the preferred models for each data set are those 
highlighted in bold in the table below: 

| Data<sup>**1**</sup>        | Model<sup>**2**</sup>                            | $\log L \pm\mathrm{S.E}$                      | $\mathrm{Pr}$                |
|--------------|-----------------------------------|---------------------------------------------|----------------------------|
| mit-3CP      | GBM                               | $\mathbf{-22,011.37\pm0.05}$                  | $\mathbf{0.74}$              |
|              | ILN                               | $-22,012.41\pm0.05$                           | $0.26$                       |
|              | STR                               | $-22,019.57\pm0.04$                           | $0.00$                       |
|              |                                   |                                             |                            |
| mit-12CP     | GBM                               | $-25,651.40\pm0.04$                           | $0.47$                       |
|              | ILN                               | $\boldsymbol{\mathbf{-25,651.28\pm0.04}}$     | $\boldsymbol{\mathbf{0.53}}$ |
|              | STR                               | $-25,657.82\pm0.03$                           | $0.00$                       |
|              |                                   |                                             |                            |
| mit-(12+3)CP | GBM                               | $-47,658.83\pm0.05$                           | $0.24$                       |
|              | ILN                               | $\mathbf{-47,657.71\pm0.05}$                  | $\mathbf{0.75}$              |
|              | STR                               | $-47,694.37\pm0.03$                           | $0.00$                       |
|              |                                   |                                             |                            |
| Morpho       | GBM - ($\mathbf{R}=\mathbf{R^{\*}}$) | $-4,097.41\pm0.04$                            | $0.00$                       |
|              | GBM - ($\mathbf{R}=\mathbf{I}$)     | $-4,221.13\pm0.04$                            | $0.00$                       |
|              | ILN - ($\mathbf{R}=\mathbf{R^{\*}}$) | $\mathbf{-4,085.03\pm0.02}$                   | $\mathbf{1.00}$              |
|              | ILN - ($\mathbf{R}=\mathbf{I}$)     | $\mathbf{-\mathrm{4,207.59}\pm\mathrm{0.02}}$ | $0.00$                       |
|              | STR - ($\mathbf{R}=\mathbf{R^{\*}}$) | $-4,158.38\pm0.01$                            | $0.00$                       |
|              | STR - ($\mathbf{R}=\mathbf{I}$)     | $-4,280.18\pm0.01$                            | $0.00$                       |

<sub><sup>**1**</sup>mit-12CP: 1 partition with the first and second codon positions (12CP) of the $12$ concatenated mitochondrial genes (12-mit genes); mit-3CP: $1$ partition with the third codon positions (3CP) of the 12-mit genes; mit-(12+3)CP: the two mitochondrial partitions analysed jointly; Morpho: $1$ partition with the morphological alignment of $87$ characters for the carnivoran data set</sub>  
<sub><sup>**2**</sup>STR: strict clock model, GBM: autocorrelated-rates model, ILN: independent-rates model, $\mathbf{R}=\mathbf{I}$: no correlation model, $\mathbf{R}=\mathbf{R^{\*}}$: correlation model. Note that, in all cases, population noise is explicitly accounted for in the models for morphological data.</sub>  




