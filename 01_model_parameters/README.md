# Simulated data sets 

We simulated quantitative characters on the phylogeny with $s=8$ species ($5$ extant and $3$ fossils) plotted below:

<p align="center">
  <img width="350" height="300" src="https://github.com/sabifo4/morpho/blob/master/figs/Fig2_Simulations_tree.jpeg">
</p>

When estimating the divergence times with the dating software `MCMCtree`, we used the following parameters:   

**Model of quantitative character evolution**: Brownian diffusion model as described in [Felsenstein, 1973](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1762641/). Accounting for within-lineage character correlation and character variation within populations is also 
integrated within the model and can be applied to the data sets. In the simulated data sets, the former was considered when studying the effect of population noise 
and both were considered when studying the effect of quantitative character correlation within lineages. We simulated our data sets to have a constant rate, hence we used the simplest morphological clock with a 
constant rate across the phylogeny, the equivalent to the strick clock (STR) for molecular data. The options used in the control file to set up this model for `MCMCtree` are `model = 0` and `clock = 1`.   

**Prior on node ages**: We used the birth-death-sequential-sampling (BDSS) process ([Stadler & Yang, 2013](https://academic.oup.com/sysbio/article/62/5/674/1684217)) for the prior on the node ages.
The BDSS parameters were chosen to generate a uniform density on the ages (figure below). These parameters are:   
   * Birth-rate: $\lambda=1$   
   * Death-rate: $\mu=1$   
   * Sampling fraction for extant species: $\rho=0$   
   * Rate of fossil sampling: $\psi=0.001$   

The options to enable this prior in the control file for `MCMCtree` are `TipDate = 1 1` and `BDparas = 1 1 0 0.001`. Note that $\rho=0$ should not be interpreted as a birth-death (BD) process without any species sampled. It is just a trick to obtain the uniform density. Note the
result here is a density similar to that of [Kishino et al., 2001](https://www.ncbi.nlm.nih.gov/pubmed/11230536); who use a Beta kernel. An example of the different 
prior densities we tested and the one with the parameters detailed above are given in the next figure:

<p align="center">
  <img width="500" height="500" src="https://github.com/sabifo4/morpho/blob/master/figs/Fig3_Test_prior.jpeg">
</p>

**Prior on the root**: The true root age was set to $t_{root}=1$ ($100$ myr time unit), hence the prior assigned was a uniform density with soft bounds between $0.8$ and $1.2$ 
(corresponding to a calibration of $80$ to $120$ Ma given our $100$ myr time unit). The option used in the control file for `MCMCtree` is `RootAge = B(0.8, 1.2)`.    

**Prior on the morphological rate**: Assuming that the true morphological rate was $r_{morpho}=1$, we used a diffuse gamma-Dirichlet distribution: $\Gamma(2,2)$, which 
in `MCMCtree` is set as `rgene_gamma = 2 2`.
Nevertheless, when we scaled the data sets when quantitative data were simulated to have population noise, the $\beta$ parameters changed for the gamma distributions 
used:   
  * When the population noise was $c=0.25$, the rate was scaled as $r_{morpho}/0.25=1/0.25=4$. Therefore, the new gamma prior was $\Gamma(2.0.5)$ (`rgene_gamma = 2 0.5`).   
  * When the population noise was $c=0.50$, the rate was scaled as $r_{morpho}/0.5=1/0.5=2$, thus the rate prior was set to $\Gamma(2,1)$ (`rgene_gamma = 2 2`).   

---
You can find the control files used in `MCMCtree` [here](https://github.com/sabifo4/morpho/tree/master/03_divtimes/MCMCtree_ctl_files/simulations) and the description of the simulations for each effect evaluated [here](https://github.com/sabifo4/morpho/tree/master/00_data/simulations).
  
# Carnivoran data: molecular data

After generating the molecular partitioned data sets, we carried out the Bayesian model selection and 
subsequent divergence times estimation analyses on three data sets:   

   * Data set with two partitions: 1st + 2nd CPs & 3rd CPs.   
   * Data set with one partition: 1st+2nd CPs.   
   * Data set with one partition: 3rd CPs.   
   
For each data set, the model parameters specified in the `MCMCtree` control file were the following:   

**Nucleotide substitution model**: HKY85 substitution model ([Hasegawa et al., 1984](https://www.jstage.jst.go.jp/article/pjab1977/60/4/60_4_95/_article) and [Hasegawa et al., 1985](https://link.springer.com/article/10.1007/BF02101694)) that allows for rate heterogeneity across nucleotide sites according to 
a discretized gamma distribution with $5$ categories (see [Yang, 1994](http://abacus.gene.ucl.ac.uk/ziheng/pdf/1994YangJMEv39p306.pdf) and [1995](http://abacus.gene.ucl.ac.uk/ziheng/pdf/1995YangGeneticsv139p993.pdf)). 
The options used in the control file for `MCMCtree` are `model = 4`, `alpha = 0.5`, and `ncatG = 5`.    

**Prior on node ages**: Prior constructed using the birth-death (BD) process ([Yang & Rannala, 2006](http://abacus.gene.ucl.ac.uk/ziheng/pdf/2006YangRannalaMBEv23p212.pdf)). The parameters 
used in `MCMCtree` are:   
   * Birth-rate: $\lambda=1$   
   * Death-rate: $\mu=1$   
   * Sampling frequency: $\rho=0.1$   

The option used to set this prior in the control file for `MCMCtree` is `BDparas = 1 1 0.1`.   

**Prior on the root**: For divergence times estimation analyses, the prior was set to a uniform fossil calibration with soft bounds between $37.3$ Ma and $66.0$ Ma ([Benton et al., 2015](https://palaeo-electronica.org/content/pdfs/424.pdf)).
The notation used in `MCMCtree` to set up this calibration in the control file is `RootAge = B(37.3, 66.0, 0.025, 0.025)`, and time unit is set to $1$ myr.
For Bayes factors analyses, the prior was $U(0.999,1.001)$ (`RootAge = '>0.999<1.001'`), fixing the age of the root to one. This is because we were not interested in divergence times estimation 
but in finding the best-fitting relaxed-clock model.   

**Prior on the molecular rate**: We used a gamma-Dirichlet prior ([dos Reis et al., 2014](https://academic.oup.com/sysbio/article/63/4/555/2848320)) so the mean of the prior rate 
(given by the mean of this distribution, $\alpha/\beta$) is close to empirical estimates 
based on the molecular branch lenghts of the phylogeny. We used `RAxML v8.2.10` ([Stamatakis, 2014](https://academic.oup.com/bioinformatics/article/30/9/1312/238053)) to estimate 
the branch lengths of the phylogeny with maximum likelihood. 
We then used the resulting maximum likelihood trees to obtain an approximation to the number of substitutions, so we could 
have a rough idea of mean molecular rates for each molecular data set. 
You can find the results of the analyses for each molecular data set and details
on the commands used in the links provided below:   

   * Data set with two partitions: [1st + 2nd CPs & 3rd CPs](https://github.com/sabifo4/morpho/blob/master/01_model_parameters/carnivoran_mol_data/mit12+3).   
   * Data set with one partition: [1st+2nd CPs](https://github.com/sabifo4/morpho/blob/master/01_model_parameters/carnivoran_mol_data/mit12).   
   * Data set with one partition: [3rd CPs](https://github.com/sabifo4/morpho/blob/master/01_model_parameters/carnivoran_mol_data/mit3).

The table below summarises the priors on the rates used for each of the three molecular data sets. 

| Analysis         | Data<sup>**1**</sup>                 | Prior on rates      | Prior on root age<sup>**2**</sup>           |
|------------------|----------------------|---------------------|------------------------------|
| Divergence times | mit-3CP              | $\Gamma(2,100)$        | $U(37.30,66.00)$               |
|                  | mit-12CP             | $\Gamma(2,1040)$       | $U(37.30,66.00)$               |
|                  | mit-(12+3)CP         | $\Gamma(2,100)$        | $U(37.30,66.00)$               |
|                  |                      |                     |                              |
| Bayes factors    | mit-3CP              | $\Gamma(2,2)$          | $U(0.999,1.001)$               |
|                  | mit-12CP             | $\Gamma(2,20)$         | $U(0.999,1.001)$               |
|                  | mit-(12+3)CP         | $\Gamma(2,2)$          | $U(0.999,1.001)$               |

<sub><sup>**1**</sup>mit-3CP: mitochondrial third codon positions; mit-12CP: mitochondrial first and second codon positions; mit-(12+3)CP: mitochondrial data with first and second codon positions in one partiton and third codon positions in another partition.</sub>  
<sub><sup>**2**</sup>Note that, in `MCMCtree`, uniform fossil calibrations have soft bounds, that is, there is a small tail probability ($p=2.5\%$ by default) that the time may lay outside each of the calibration bounds.</sub>  

The option used in the control file for `MCMCtree` is `rgene_gamma = alpha beta`, where the user can replace `alpha` and `beta` with the corresponding values 
of $\alpha$ and $\beta$. Note that the option `ndata = NDATA` is used to define the number of partitions, hence `NDATA` was set to `1` if there was $1$ molecular
partition or `2` for $2$ molecular partitions.   

**Clock model**: We ran all the molecular data sets under three clock models: the strict clock (STR), the geometric Brownian diffusion (GBM, also known as autocorrelated-rates, [Thorne et al., 1998](https://academic.oup.com/mbe/article/15/12/1647/963101); [Yang & Rannala, 2006](http://abacus.gene.ucl.ac.uk/ziheng/pdf/2006YangRannalaMBEv23p212.pdf)),
and independent log-normal rate (ILN) models ([Rannala & Yang, 2007](https://academic.oup.com/sysbio/article/56/3/453/1657118); [Lemey2010](https://academic.oup.com/mbe/article/27/8/1877/988944)). We then estimated the marginal likelihood in order to 
calculate Bayes factors and select the best-fitting clock model. After that, Bayesian divergence times estimation proceeded on the preferred clock model selected in the 
Bayesian model selection analysis. The option used in the control file for `MCMCtree` is `clock = CLOCK`, where `CLOCK` can be replaced with `1` for the STR
clock model, with `2` for the ILN model, or with `3` for the GBM model. 
Note that the gamma-Dirichlet prior on $\sigma_{i}^{2}$ for the GBM and ILN models was set to $\Gamma(2,2)$. This is set in the control file for `MCMCtree` as 
`sigma2_gamma = 2 2`.

---
You can find the control files for the Bayesian model selection analyses [here](https://github.com/sabifo4/morpho/tree/master/02_bayes.model.sel/MCMCtree_ctl_files/carnivoran_mol_data)
and for the Bayesian divergence times estimation analyses [here](https://github.com/sabifo4/morpho/tree/master/03_divtimes/MCMCtree_ctl_files/carnivoran_mol_data).

# Carnivoran data: morphological data 

When estimating the divergence times with the dating software `MCMCtree`, we used the following parameters:   

**Model of quantitative character evolution**: Brownian diffusion model as described in [Felsenstein, 1973](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1762641/). Accounting for within-lineage character correlation and character variation within populations is also 
integrated within the model and can be applied to the data sets. In the simulated data sets, the former was considered when studying the effect of population noise 
and both were considered when studying the effect of quantitative character correlation within lineages. The option used in the control file for `MCMCtree` is `model = 0`.   

**Prior on node ages**: We used the birth-death-sequential-sampling (BDSS) process ([Stadler & Yang, 2013](https://academic.oup.com/sysbio/article/62/5/674/1684217)) for the prior on the node ages.
The BDSS parameters were chosen to generate a uniform density on the ages, as it was done with the simulated data sets (see the corresponding section above 
for more details). These parameters are:   
   * Birth-rate: $\lambda=1$   
   * Death-rate: $\mu=1$   
   * Sampling fraction for extant species: $\rho=0$   
   * Rate of fossil sampling: $\psi=0.001$   

The options to enable this prior in the control file for `MCMCtree` are `TipDate = 1 1` and `BDparas = 1 1 0 0.001`.   

**Prior on the root**: Uniform fossil calibration with soft bounds between $37.3$ Ma and $66.0$ Ma ([Benton et al., 2015](https://palaeo-electronica.org/content/pdfs/424.pdf)).
The notation used in `MCMCtree` to set up this calibration in the control file is `RootAge = B(37.3, 66.0, 0.025, 0.025)`, and time unit is set to $1$ myr.   

**Prior on the morphological rate**: We used a gamma-Dirichlet prior ([dos Reis et al., 2014](https://academic.oup.com/sysbio/article/63/4/555/2848320)) so the mean of the prior rate 
(given by the mean of this distribution, $\alpha/\beta$) is close to empirical estimates 
based on the morphological branch lenghts of the phylogeny. We used `CONTML` (`PHYLIP` package, [Felsenstein, 1993](http://www0.nih.go.jp/~jun/research/phylip/main.html)) to estimate 
the branch lengths of the phylogeny using maximum likelihood, which were used to obtain an approximation to the units of morphological drift, so we could 
have an idea of mean morphological rates (see details about the software [here](https://github.com/sabifo4/morpho/tree/master/01_model_parameters/carnivoran_morpho_data/CONTML)). 
You can find the results of the analyses for the morphological data sets (ignoring or accounting for within-lineage character correlation) and details
on the commands used in the links provided below:   

   * [Accounting for within-lineage character correlation](https://github.com/sabifo4/morpho/tree/master/01_model_parameters/carnivoran_morpho_data/CONTML/with_correlation).   
   * [Ignoring within-lineage character correlation](https://github.com/sabifo4/morpho/tree/master/01_model_parameters/carnivoran_morpho_data/CONTML/without_correlation).   

The table below summarises the priors on the rates used for both morphological data sets. 

| Analysis         | Data<sup>**1**</sup>                 | Prior on rates      | Prior on root age<sup>**2**</sup>           |
|------------------|----------------------|---------------------|------------------------------|
| Divergence times | morpho               | $\Gamma(2,5)$          | $U(37.30,66.00)$               |
|                  |                      |                     |                              |
| Bayes factors    | morpho               | $\Gamma(2,5)$          | $U(37.30,66.00)$               |

<sub><sup>**1**</sup>morpho: morphological data.</sub>  
<sub><sup>**2**</sup>Note that in `MCMCtree`, uniform fossil calibrations have soft bounds, that is, there is a small tail probability ($p=2.5\%$ by default) that the time may lay outside each of the calibration bounds.</sub>  

**Clock model**: We ran all the morphological data sets under three clock models: the strict clock (STR), the geometric Brownian diffusion (GBM, also known as autocorrelated-rates, [Thorne et al., 1998](https://academic.oup.com/mbe/article/15/12/1647/963101); [Yang & Rannala, 2006](http://abacus.gene.ucl.ac.uk/ziheng/pdf/2006YangRannalaMBEv23p212.pdf)),
and independent log-normal rate (ILN) models ([Rannala & Yang, 2007](https://academic.oup.com/sysbio/article/56/3/453/1657118); [Lemey2010](https://academic.oup.com/mbe/article/27/8/1877/988944)). We then estimated the marginal likelihood in order to 
calculate Bayes factors and select the best-fitting clock model. After that, Bayesian divergence times estimation proceeded on the preferred clock model selected in the 
Bayesian model selection analysis. The option used in the control file for `MCMCtree` is `clock = CLOCK`, where `CLOCK` can be replaced with `1` for the STR
clock model, with `2` for the ILN model, or with `3` for the GBM model.  
Note that the gamma-Dirichlet prior on $\sigma_{i}^{2}$ for the GBM and ILN models was set to $\Gamma(2,2)$. This is set in the control file for `MCMCtree` as `sigma2_gamma = 2 2`.

---
You can find the control files for the Bayesian model selection analyses [here](https://github.com/sabifo4/morpho/tree/master/02_bayes.model.sel/MCMCtree_ctl_files/carnivoran_morpho_data)
and for the Bayesian divergence times estimation analyses [here](https://github.com/sabifo4/morpho/tree/master/03_divtimes/MCMCtree_ctl_files/carnivoran_morpho_data).


# Carnivoran data: combined data set (molecular + morphological data)

When estimating the divergence times with the dating software `MCMCtree`, we used the following parameters:   

**Model of morphological and molecular data evolution**:   
   * **Morphological partition**: Brownian diffusion model as described in [Felsenstein, 1973](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1762641/). Accounting for within-lineage character correlation and character variation within populations is also 
integrated within the model and can be applied to the data sets. In the simulated data sets, the former was considered when studying the effect of population noise 
and both were considered when studying the effect of quantitative character correlation within lineages. This model assumes the simplest morphological clock with a 
constant rate across the phylogeny, the equivalent to the strick clock (STR) for molecular data. The options used in the control file for `MCMCtree` are `model = 0` and `clock = 1`.    
   * **Molecular partition**: HKY85 substitution model ([Hasegawa et al., 1984](https://www.jstage.jst.go.jp/article/pjab1977/60/4/60_4_95/_article) and [Hasegawa et al., 1985](https://link.springer.com/article/10.1007/BF02101694)) that allows for rate heterogeneity across nucleotide sites according to 
a discretized gamma distribution with $5$ categories (see [Yang, 1994](http://abacus.gene.ucl.ac.uk/ziheng/pdf/1994YangJMEv39p306.pdf) and [1995](http://abacus.gene.ucl.ac.uk/ziheng/pdf/1995YangGeneticsv139p993.pdf)).
The options used in the control file for `MCMCtree` are `model = 4`, `alpha = 0.5` and `ncatG = 5`.   

Note that the option `ndata = NDATA` is used to define the number of partitions, hence `NDATA` was set to `3` as we have two molecular partitions (1st+2nd CPs & 3rd CPs) and one morphological partition.   

**Prior on node ages**: We used the birth-death-sequential-sampling (BDSS) process ([Stadler & Yang, 2013](https://academic.oup.com/sysbio/article/62/5/674/1684217)) for the prior on the node ages.
The BDSS parameters were chosen to generate a uniform density on the ages, as it was done with the simulated data sets. These parameters are:   
   * Birth-rate: $\lambda=1$   
   * Death-rate: $\mu=1$   
   * Sampling fraction for extant species: $\rho=0$   
   * Rate of fossil sampling: $\psi=0.001$

The options to enable this prior in the control file for `MCMCtree` are `TipDate = 1 1` and `BDparas = 1 1 0 0.001`.   

**Prior on the root**: Uniform fossil calibration with soft bounds between $37.3$ Ma and $66.0$ Ma ([Benton et al., 2015](https://palaeo-electronica.org/content/pdfs/424.pdf)).
The notation used in `MCMCtree` to set up this calibration in the control file is `RootAge = B(37.3, 66.0, 0.025, 0.025)`, and time unit is set to $1$ myr.   

**Prior on the rate**: We used a gamma-Dirichlet prior ([dos Reis et al., 2014](https://academic.oup.com/sysbio/article/63/4/555/2848320)) so the mean of the prior rate 
(given by the mean of this distribution, $\alpha/\beta$) is close to empirical estimates 
based on the morphological branch lenghts of the phylogeny. The mean rate for this combined data set was calculated as the average of the mean rates of the
morphological partition alone (previously estimated, $r_{morpho}=\alpha/\beta=2/5~0.4$)
and the molecular partitions (previously estimated, $r_{mit12+3}=\alpha/\beta=2/100=0.02$), that is $r_{joint}=(0.4+0.02)/2=0.21$. This means that the option used in the control
file for `MCMCtree` is `rgene_gamma = 2 10`.
More details about the procedure followed [here](https://github.com/sabifo4/morpho/tree/master/01_model_parameters/carnivoran_morpho_mol_data). 
The next table contains the prior on the rates used for divergence times estimation for the combined data set. Note that we did not carry out any 
Bayesian model selection as we used the preferred clock model found for the individual data sets: the independent log-normal rate (ILN) models ([Rannala & Yang, 2007](https://academic.oup.com/sysbio/article/56/3/453/1657118), [Lemey2010](https://academic.oup.com/mbe/article/27/8/1877/988944)). 

| Analysis         | Data<sup>**1**</sup>                 | Prior on rates      | Prior on root age<sup>**2**</sup>           |
|------------------|----------------------|---------------------|------------------------------|
| Divergence times | morpho+mit-(12+3)CP  | $\Gamma(2,10)$         | U(37.30,66.00)               |

<sub><sup>**1**</sup>morpho+mit-(12+3)CP: morphological and molecular data in three partitions.</sub>  
<sub><sup>**2**</sup>Note that in `MCMCtree`, uniform fossil calibrations have soft bounds, that is, there is a small tail probability ($p=2.5\%$ by default) that the time may lay outside each of the calibration bounds.</sub>  

**Clock model**: We ran all the combined data set under the preferred model for both morphological and molecular data set in previous analyses: 
the independent log-normal rate (ILN) models ([Rannala & Yang, 2007](https://academic.oup.com/sysbio/article/56/3/453/1657118), [Lemey2010](https://academic.oup.com/mbe/article/27/8/1877/988944)). 
Note that the gamma-Dirichlet prior on $\sigma_{i}^{2}$ for the GBM and ILN models was set to $\Gamma(2,2)$. The options in the control file for `MCMCtree` are 
`clock = 2` and `sigma2_gamma = 2 2`.

---
You can find the control file [here](https://github.com/sabifo4/morpho/tree/master/03_divtimes/MCMCtree_ctl_files/carnivoran_morpho_mol_data).




