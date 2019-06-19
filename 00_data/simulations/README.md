# Simulations

## Simulated data: an overview

We used the phylogeny shown below with $s=8$ species ($5$ extant and $3$ fossils) to simulate quantitative 
morphological traits:
 
<p align="center">
  <img width="400" height="300" src="https://github.com/sabifo4/morpho/blob/master/figs/Fig2_Simulations_tree.jpeg">
</p>

The morphological evolutionary rate is $r=1$ and constant along all the branches of the phylogeny. 
The simulated data matrices are generated under the Brownian diffusion model using the [`mcmc3r` R package](https://github.com/dosreislab/mcmc3r). 

In our simulations, we evaluated four different scenarios: effect of the number of characters, fossil age, population 
noise, and trait correlation. A total of $R=1,000$ replicates were simulated per setup.
An overview of the simulations is given in the next figure:

<p align="center">
  <img width="800" height="500" src="https://github.com/sabifo4/morpho/blob/master/figs/Fig3_Simulworkflow.png">
</p>

A more detailed explanation for each effect evaluated for this simulation analysis and links to 
the R code are given in the next section.

## Generating the simulated data sets

### 1. Effect of the number of characters 

We simulated data sets with $p=100$, $p=1,000$; and $p=10,000$ characters under the 
phylogeny detailed above ($s=8$ species). We assume independence among characters and no population noise ($c=0$).

The R script used for these simulations can be found [here](https://github.com/sabifo4/morpho/blob/master/00_data/simulations/R_scripts/00_simulate_num_chars.R). 

#### Summary of simulated matrices using `mcmc3r`
```
Matrices of size s x p
====================================================

  M1 ( size 8 x 100 )    | R = 1,000 replicates 
  M2 ( size 8 x 1000 )   | R = 1,000 replicates 
  M3 ( size 8 x 10,000 ) | R = 1,000 replicates 
```

### 2. Effect of fossil age 

We then varied the age of fossil $H$ so it ranged from younger to older: $t_{H}=0.7$, $0.5$, $0.3$, and $0.1$. Note that 
the rest of the fossil ages remained unchanged. For each tree topology with the different 
fossil age for species $H$, we simulated data sets with $p=100$, $p=1,000$; and $p=10,000$ characters. We assume 
independence among characters and no population noise ($c=0$).

The R script for these simulations can be found [here](https://github.com/sabifo4/morpho/blob/master/00_data/simulations/R_scripts/01_simulate_fossil_age.R). 

#### Summary of simulated matrices using `mcmc3r`
```
Matrices of size s x p per tree topology with changing fossil age for species H
================================================================================

A. Matrices when phylogeny with t_H = 0.7
  * M1.tH0.7 ( size 8 x 100 )    | R = 1,000 replicates 
  * M2.tH0.7 ( size 8 x 1000 )   | R = 1,000 replicates 
  * M3.tH0.7 ( size 8 x 10,000 ) | R = 1,000 replicates
  
B. Matrices when phylogeny with t_H = 0.5
  * M1.tH0.5 ( size 8 x 100 )    | R = 1,000 replicates 
  * M2.tH0.5 ( size 8 x 1000 )   | R = 1,000 replicates 
  * M3.tH0.5 ( size 8 x 10,000 ) | R = 1,000 replicates
  
C. Matrices when phylogeny with t_H = 0.3 [ These were already simulated for effect 1 ]
  * M1.tH0.3 ( size 8 x 100 )    | R = 1,000 replicates 
  * M2.tH0.3 ( size 8 x 1000 )   | R = 1,000 replicates 
  * M3.tH0.3 ( size 8 x 10,000 ) | R = 1,000 replicates
  
D. Matrices when phylogeny with t_H = 0.1
  * M1.tH0.1 ( size 8 x 100 )    | R = 1,000 replicates 
  * M2.tH0.1 ( size 8 x 1000 )   | R = 1,000 replicates 
  * M3.tH0.1 ( size 8 x 10,000 ) | R = 1,000 replicates
```

### 3. Effect of population noise 

We simulated data sets with low and high population noise ($c=0.25$ and $c=0.50$, respectively)
with $p=100$, $p=1,000$; and $p=10,000$ characters, which are assumed to evolve independently.

#### 3.1. How is population noise added?
We sampled $p$ random numbers from a normal distribution with mean $\mu=0$ and variance $c$ for $s=8$ 
species under the phylogeny detailed above. The resulting sampled vectors of size $p$ are then 
used to obtain the noise matrix of size $s\times p$. This noise is then added to the simulated morphological 
matrix of size $s\times p$. 

#### Summary of simulated matrices using `mcmc3r`
```
A. Generate noise matrices for low and high noise 
  * N1.c0.25 ( size 8 x 100 )    | R = 1,000 replicates 
  * N2.c0.25 ( size 8 x 1000 )   | R = 1,000 replicates 
  * N3.c0.25 ( size 8 x 10,000 ) | R = 1,000 replicates 
  
  * N1.c0.50 ( size 8 x 100 )    | R = 1,000 replicates 
  * N2.c0.50 ( size 8 x 1000 )   | R = 1,000 replicates 
  * N3.c0.50 ( size 8 x 10,000 ) | R = 1,000 replicates 
  
B. Generate morphological matrices
  * M1 ( size 8 x 100 )    | R = 1,000 replicates 
  * M2 ( size 8 x 1000 )   | R = 1,000 replicates 
  * M3 ( size 8 x 10,000 ) | R = 1,000 replicates 
  
C. Generate noisy matrices 
  * Mn1.c0.25 = M1 + N1.c0.25 ( size 8 x 100 )  | For each of the 1,000 replicates
  * Mn2.c0.25 = M2 + N2.c0.25 ( size 8 x 1000 ) | For each of the 1,000 replicates
  * Mn3.c0.25 = M3 + N3.c0.25 ( size 8 x 1000 ) | For each of the 1,000 replicates
  
  * Mn1.c0.50 = M1 + N1.c0.50 ( size 8 x 100 )  | For each of the 1,000 replicates
  * Mn2.c0.50 = M2 + N2.c0.50 ( size 8 x 1000 ) | For each of the 1,000 replicates
  * Mn3.c0.50 = M3 + N3.c0.50 ( size 8 x 1000 ) | For each of the 1,000 replicates
```

#### 3.2. How can we correct for population noise? 
For this part of the simulation, we use a sample of individuals for one of the 
species included in the phylogeny to get an estimate of the population noise for each character.
We assume that this population noise is more or less similar for the rest of the species in the phylogeny,
so we use the vector of estimated population variances to scale the matrix with the quantitative data
for all the species in the phylogeny.   

According to this assumption, we simulated a population sample of $n=20$ individuals to obtain a population 
matrix $\mathbf{P}$ of size $n\times p$, by sampling from a normal distribution with mean $\mu=0$ and 
variance $c$. 

Then, we obtained a vector with the estimated variances ($\hat{\mathbf{c}}$) for each character on matrix $\mathbf{P}$ and used 
it to scale the noisy matrix, $\mathbf{M}^{(n)}$, such as $\mathbf{M}^{(s)}=\mathbf{M}^{(n)}\mathrm{diag}\{1/\sqrt{\hat{\mathbf{c}}}\}$.
We also scaled the noisy matrix by the vector of true variances, $\mathbf{c}$, so we could have a a control test.

>###### **NOTE**: When $c=0.25$, the morphological rate for the scaled data is $r/0.25=1/0.25=4$. Similarly, when $c=0.50$, the morphological rate for the scaled data is $r/0.50=1/0.50=2$. This means that, during Bayesian inference, the rate priors change.
>###### We provide more information about the parameters for the control files used by `MCMCtree` for the simulated data sets [here](https://github.com/sabifo4/morpho/tree/master/01_model_parameters).

#### Summary of simulated matrices using `mcmc3r`
```
A. Generate population matrices for low and high noise 
  * P1.c0.25 ( size 8 x 100 )    | R = 1,000 replicates 
  * P2.c0.25 ( size 8 x 1000 )   | R = 1,000 replicates 
  * P3.c0.25 ( size 8 x 10,000 ) | R = 1,000 replicates 
  
  * P1.c0.50 ( size 8 x 100 )    | R = 1,000 replicates 
  * P2.c0.50 ( size 8 x 1000 )   | R = 1,000 replicates 
  * P3.c0.50 ( size 8 x 10,000 ) | R = 1,000 replicates 
  
B. Get vectors with estimated population noise
  * estc1.c0.25 ( size 100 )    | Obtained on P1.c0.25 | For each of the R = 1,000 replicates 
  * estc2.c0.25 ( size 1000 )   | Obtained on P2.c0.25 | For each of the R = 1,000 replicates 
  * estc3.c0.25 ( size 10,000 ) | Obtained on P3.c0.25 | For each of the R = 1,000 replicates 
  
  * estc1.c0.50 ( size 100 )    | Obtained on P1.c0.50 | For each of the R = 1,000 replicates 
  * estc2.c0.50 ( size 1000 )   | Obtained on P2.c0.50 | For each of the R = 1,000 replicates 
  * estc3.c0.50 ( size 10,000 ) | Obtained on P3.c0.50 | For each of the R = 1,000 replicates 
  
C. Scale matrices with vectors with estimated population noise
  * Ms1.estc0.25 = Mn1.c0.25 * diag{1/sqrt(estc1.c0.25)} ( size 8 x 100 )  | For each of the 1,000 replicates
  * Ms2.estc0.25 = Mn2.c0.25 * diag{1/sqrt(estc2.c0.25)} ( size 8 x 1000 ) | For each of the 1,000 replicates
  * Ms3.estc0.25 = Mn3.c0.25 * diag{1/sqrt(estc3.c0.25)} ( size 8 x 1000 ) | For each of the 1,000 replicates
  
  * Ms1.estc0.50 = Mn1.c0.50 * diag{1/sqrt(estc1.c0.50)} ( size 8 x 100 )  | For each of the 1,000 replicates
  * Ms2.estc0.50 = Mn2.c0.50 * diag{1/sqrt(estc2.c0.50)} ( size 8 x 1000 ) | For each of the 1,000 replicates
  * Ms3.estc0.50 = Mn3.c0.50 * diag{1/sqrt(estc3.c0.50)} ( size 8 x 1000 ) | For each of the 1,000 replicates

D. Scale matrices with vectors with true population noise
  * Ms1.c0.25 = Mn1.c0.25 * diag{1/sqrt(c0.25)} ( size 8 x 100 )  | For each of the 1,000 replicates
  * Ms2.c0.25 = Mn2.c0.25 * diag{1/sqrt(c0.25)} ( size 8 x 1000 ) | For each of the 1,000 replicates
  * Ms3.c0.25 = Mn3.c0.25 * diag{1/sqrt(c0.25)} ( size 8 x 1000 ) | For each of the 1,000 replicates
  
  * Ms1.c0.50 = Mn1.c0.50 * diag{1/sqrt(c0.50)} ( size 8 x 100 )  | For each of the 1,000 replicates
  * Ms2.c0.50 = Mn2.c0.50 * diag{1/sqrt(c0.50)} ( size 8 x 1000 ) | For each of the 1,000 replicates
  * Ms3.c0.50 = Mn3.c0.50 * diag{1/sqrt(c0.50)} ( size 8 x 1000 ) | For each of the 1,000 replicates
```

The R script for these simulations can be found [here](https://github.com/sabifo4/morpho/blob/master/00_data/simulations/R_scripts/02_simulate_popnoise.R). 

### 4. Effect of within-lineage character correlation

We simulated data sets using the constant correlation model, which means that all the 
within-lineage correlations in the correlation matrix $\mathbf{R}$ are equal to $\rho$.
We used different correlation values for $\rho$ ($\rho=0$, $0.25$, $0.35$, $0.5$, $0.70$, $0.80$, $0.9$)
and simulated data sets with $p=100$, $p=1,000$; and $p=10,000$ characters, which are assumed to evolve independently.
Note that the data are simulated on the tree, thus it already contains the among-lineage covariance induced by the
shared ancestry. We also fixed the population noise to be $c=0.25$.

#### 4.1. How is within-lineage character correlation added?
Once we had the correlation matrix $\mathbf{R}$ with within-lineage correlations equal to $\rho$, 
we obtained the lower triangular Cholesky decomposition, $\mathbf{L}$. We then added the 
within-lineage correlation to the matrix by multiplying the matrix with quantitative 
data $\mathbf{M}$ by the transposed of $\mathbf{L}$, such as $\mathbf{M}^{\mathrm{(R)}}=\mathbf{M}\mathbf{L^{\mathrm{T}}}$.

#### Summary of simulated matrices using `mcmc3r`
```
*NOTE: The same procedure is repeated for the 7 values of rho tested. The "X" after "rhoX"
corresponds to "rho0", "rho0.25", "rho0.35", "rho0.50", "rho0.70", "rho0.80", "rho0.90"

A. Generate noise matrices for low and high noise. 
  * R1.rhoX ( size 100 x 100 )       | R = 1,000 replicates 
  * R2.rhoX ( size 1000 x 1000 )     | R = 1,000 replicates 
  * R3.rhoX ( size 10,000 x 10,000 ) | R = 1,000 replicates 
  
B. Get lower triangular matrix using Cholesky decomposition
  * L1.rhoX ( p = 100 )    | For each of the R = 1,000 replicates 
  * L2.rhoX ( p = 1000 )   | For each of the R = 1,000 replicates 
  * L3.rhoX ( p = 10,000 ) | For each of the R = 1,000 replicates 
  
C. Generate matrices with within-lineage character correlation 
  * M1.rhoX = M1 * t(L1.rhoX) ( size 8 x 100 )  | For each of the 1,000 replicates
  * M2.rhoX = M2 * t(L2.rhoX) ( size 8 x 1000 ) | For each of the 1,000 replicates
  * M3.rhoX = M3 * t(L3.rhoX) ( size 8 x 1000 ) | For each of the 1,000 replicates
```

#### 4.2. How is population noise added?
As we did in the previous simulation setup, we sampled $p$ characters from a normal distribution
with mean $\mu=0$ and variance $c=0.25$ for each individual and obtain the noise matrix of size $s\times p$. 
Then, within-lineage correlation is also added to this noise matrix using the $\mathbf{L}$ matrix, such 
as $\mathbf{N}^{\mathrm{(R)}}=\mathbf{N}\mathbf{L^{\mathrm{T}}}$. This $\mathbf{N}^{\mathrm{(R)}}$ matrix is then added to the matrix with quantitative data to which 
correlation had been previously added, $\mathbf{M^{\mathrm{(R)}}}$, such as $\mathbf{M}^{(n)}=\mathbf{M^{\mathrm{(R)}}}+\mathbf{N^{\mathrm{(R)}}}$.

#### Summary of simulated matrices using `mcmc3r`
```
*NOTE: The same procedure is repeated for the 7 values of rho tested. The "X" after "rhoX"
corresponds to "rho0", "rho0.25", "rho0.35", "rho0.50", "rho0.70", "rho0.80", "rho0.90"

A. Generate noise matrices for a fixed population noise, c = 0.25 
  * N1.c0.25 ( size 8 x 100 )    | R = 1,000 replicates 
  * N2.c0.25 ( size 8 x 1000 )   | R = 1,000 replicates 
  * N3.c0.25 ( size 8 x 10,000 ) | R = 1,000 replicates 
  
B. Add within-lineage correlation to noise matrices using the lower triangular matrix previously computed
  * N1.c0.25.rhoX = N1.c0.25 * t(L1.rhoX) ( size 8 x 100 )    | For each of the R = 1,000 replicates 
  * N2.c0.25.rhoX = N2.c0.25 * t(L2.rhoX) ( size 8 x 1000 )   | For each of the R = 1,000 replicates 
  * N3.c0.25.rhoX = N3.c0.25 * t(L3.rhoX) ( size 8 x 10,000 ) | For each of the R = 1,000 replicates 
  
C. Generate matrices with within-lineage character correlation and noise added
  * Mn1.rhoX = M1.rhoX + N1.c0.25.rhoX ( size 8 x 100 )  | For each of the 1,000 replicates
  * Mn2.rhoX = M2.rhoX + N2.c0.25.rhoX ( size 8 x 1000 ) | For each of the 1,000 replicates
  * Mn3.rhoX = M3.rhoX + N3.c0.25.rhoX ( size 8 x 1000 ) | For each of the 1,000 replicates
```

#### 4.3. How can we correct for population noise and character correlation? 
Following the previous simulation setup, we simulated a population sample of $n=20$ individuals to obtain a population 
matrix $\mathbf{P}$ of size $n\times p$, by sampling from a normal distribution with mean $\mu=0$ and 
variance $c$. We then added within-lineage correlation by multiplying $\mathbf{P}$ by the transposed of 
the $\mathbf{L}$ matrix, such as $\mathbf{P}^{\mathrm{(R)}}=\mathbf{P}\mathbf{L^{\mathrm{T}}}$.

Then, we obtained a vector with the estimated variances for each character on matrix $\mathbf{P^{\mathrm{(R)}}}$ and used 
it to obtain the scaled data matrix, $\mathbf{M}^{(s)}$.  We also scaled the noisy matrix by the vector of true 
variances so we could have a control test.

We cannot use the unbiased correlation matrix to correct for within-lineage character correlation 
because these matrices tend to be singular. As we need the invert of $\mathbf{R}$ during the likelihood calculation
when traversing the phylogeny when Bayesian inference takes place, we need an invertible $\mathbf{R}$. Therefore,
we used a shrinkage approach to obtain an estimate of $\mathbf{R}$ with the `corpcor` R package. The algorithm 
used in this package generates a shrinkage correlation matrix $\mathbf{R}^{\*}$ based on a $\delta$ value that 
controls the level of shrinkage ranging from the identity matrix $\mathbf{I}$ to the unbiased correlation 
matrix $\mathbf{\hat{R}}$: $\mathbf{R}^{\*}=\delta\mathbf{I}+(1-\delta)\mathbf{\hat{R}}$. 
As $\delta$ has a strong impact on the final $\mathbf{R}^{\*}$, we decided to test two methods to obtain $\mathbf{R}^{\*}$.
We first used the automatic approach within this package, which finds the optimum value of 
$\delta$ to be used to generate the shrinkage matrix $\mathbf{R}^{\*}$. Then, we fixed the value of $\delta$ 
to $\delta=0.01$ so we could obtain an estimate very close to the unbiased correlation matrix $\mathbf{\hat{R}}$. 

After that, we wanted to find a matrix $\mathbf{A}$ such as $\mathbf{R}^{-1}=\mathbf{A}^{\mathrm{T}}\mathbf{A}$. This was done for the true correlation 
matrix $\mathbf{R}$ and for the shrinkage matrices obtained with the method described above, $\mathbf{R}\_{\delta_{0.01}}^{\*}$ (matrix estimated using $\delta$ fixed to $0.01$) 
and $\mathbf{R}\_{\delta_{DEF}}^{\*}$ (matrix estimated using the optimum $\delta$ value found by `corpcor`). Once $\mathbf{A}$ was found, then we 
could obtain the transformed data matrix in which correlation has been accounted for by computing
$\mathbf{Z}^{(s)}=\mathbf{M}^{(s)}\mathbf{A}^{\mathrm{T}}$.

#### Summary of simulated matrices using `mcmc3r`
```
*NOTE: The same procedure is repeated for the 7 values of rho tested. The "X" after "rhoX"
corresponds to "rho0", "rho0.25", "rho0.35", "rho0.50", "rho0.70", "rho0.80", "rho0.90"

A. Generate population matrices for low and high noise 
  * P1.c0.25 ( size 8 x 100 )    | R = 1,000 replicates 
  * P2.c0.25 ( size 8 x 1000 )   | R = 1,000 replicates 
  * P3.c0.25 ( size 8 x 10,000 ) | R = 1,000 replicates 

B. Add within-lineage correlation to population matrices
  * P1.c0.25.rhoX = P1.c0.25 * t(L1.rhoX) ( size 8 x 100 )    | For each of the R = 1,000 replicates 
  * P2.c0.25.rhoX = P2.c0.25 * t(L2.rhoX) ( size 8 x 1000 )   | For each of the R = 1,000 replicates 
  * P3.c0.25.rhoX = P3.c0.25 * t(L3.rhoX) ( size 8 x 10,000 ) | For each of the R = 1,000 replicates 
  
C. Get vectors with estimated population noise
  * estc1.c0.25.rhoX ( size 100 )    | Obtained on P1.c0.25.rhoX | For each of the R = 1,000 replicates 
  * estc2.c0.25.rhoX ( size 1000 )   | Obtained on P2.c0.25.rhoX | For each of the R = 1,000 replicates 
  * estc3.c0.25.rhoX ( size 10,000 ) | Obtained on P3.c0.25.rhoX | For each of the R = 1,000 replicates 
  
D. Scale matrices with vectors with estimated population noise
  * Ms1.estc0.25.rhoX = Mn1.rhoX * diag{1/sqrt(estc1.c0.25.rhoX)} ( size 8 x 100 )  | For each of the 1,000 replicates
  * Ms2.estc0.25.rhoX = Mn2.rhoX * diag{1/sqrt(estc2.c0.25.rhoX)} ( size 8 x 1000 ) | For each of the 1,000 replicates
  * Ms3.estc0.25.rhoX = Mn3.rhoX * diag{1/sqrt(estc3.c0.25.rhoX)} ( size 8 x 1000 ) | For each of the 1,000 replicates

E. Scale matrices with vectors with true population noise
  * Ms1.c0.25.rhoX = Mn1.rhoX * diag{1/sqrt(c0.25.rhoX)} ( size 8 x 100 )  | For each of the 1,000 replicates
  * Ms2.c0.25.rhoX = Mn2.rhoX * diag{1/sqrt(c0.25.rhoX)} ( size 8 x 1000 ) | For each of the 1,000 replicates
  * Ms3.c0.25.rhoX = Mn3.rhoX * diag{1/sqrt(c0.25.rhoX)} ( size 8 x 1000 ) | For each of the 1,000 replicates

F. Find a matrix "A" such as "R^(-1) = t(A) * A". This is done for the true R, for R* estimated using delta=0.01,
   and for R* estimated using the optimum delta found by the R package corpcor.
  * R1.rhoX = L1.rhoX * t(L1.rhoX) --> A1.trueR1 = L1.rhoX^(-1) ( size 100 x 100 )        | For each of the 1,000 replicates
  * R2.rhoX = L2.rhoX * t(L2.rhoX) --> A2.trueR2 = L2.rhoX^(-1) ( size 1000 x 1000 )      | For each of the 1,000 replicates
  * R3.rhoX = L3.rhoX * t(L3.rhoX) --> A3.trueR3 = L3.rhoX^(-1) ( size 10,000 x 10,000 )  | For each of the 1,000 replicates

  * R1*d0.01.rhoX = L1*d0.01.rhoX * t(L1*d0.01.rhoX) --> A1.R1*d0.01.rhoX = L1*d0.01.rhoX^(-1) ( size 100 x 100 )        | For each of the 1,000 replicates
  * R2*d0.01.rhoX = L2*d0.01.rhoX * t(L2*d0.01.rhoX) --> A2.R2*d0.01.rhoX = L2*d0.01.rhoX^(-1) ( size 1000 x 1000 )      | For each of the 1,000 replicates
  * R3*d0.01.rhoX = L3*d0.01.rhoX * t(L3*d0.01.rhoX) --> A3.R3*d0.01.rhoX = L3*d0.01.rhoX^(-1) ( size 10,000 x 10,000 )  | For each of the 1,000 replicates

  * R1*dDef.rhoX = L1*dDef.rhoX * t(L1*dDef.rhoX) --> A1.R1*dDef.rhoX = L1*dDef.rhoX^(-1) ( size 100 x 100 )        | For each of the 1,000 replicates
  * R2*dDef.rhoX = L2*dDef.rhoX * t(L2*dDef.rhoX) --> A2.R2*dDef.rhoX = L2*dDef.rhoX^(-1) ( size 1000 x 1000 )      | For each of the 1,000 replicates
  * R3*dDef.rhoX = L3*dDef.rhoX * t(L3*dDef.rhoX) --> A3.R3*dDef.rhoX = L3*dDef.rhoX^(-1) ( size 10,000 x 10,000 )  | For each of the 1,000 replicates

G. Use "A" to obtain the transformed data matrix in which correlation has been accounted for
  * Z1.c0.25.trueR1.rhoX = Ms1.estc0.25.rhoX * t(A1.trueR1) ( size 8 x 100 )  | For each of the 1,000 replicates
  * Z2.c0.25.trueR2.rhoX = Ms2.estc0.25.rhoX * t(A2.trueR2) ( size 8 x 1000 ) | For each of the 1,000 replicates
  * Z3.c0.25.trueR3.rhoX = Ms3.estc0.25.rhoX * t(A3.trueR3) ( size 8 x 1000 ) | For each of the 1,000 replicates

  * Z1.c0.25.R1*d0.01.rhoX = Ms1.estc0.25.rhoX * t(A1.R1*d0.01) ( size 8 x 100 )  | For each of the 1,000 replicates
  * Z2.c0.25.R2*d0.01.rhoX = Ms2.estc0.25.rhoX * t(A2.R2*d0.01) ( size 8 x 1000 ) | For each of the 1,000 replicates
  * Z3.c0.25.R3*d0.01.rhoX = Ms3.estc0.25.rhoX * t(A3.R3*d0.01) ( size 8 x 1000 ) | For each of the 1,000 replicates

  * Z1.c0.25.R1*dDef.rhoX = Ms1.estc0.25.rhoX * t(A1.R1*dDef) ( size 8 x 100 )  | For each of the 1,000 replicates
  * Z2.c0.25.R2*dDef.rhoX = Ms2.estc0.25.rhoX * t(A2.R2*dDef) ( size 8 x 1000 ) | For each of the 1,000 replicates
  * Z3.c0.25.R3*dDef.rhoX = Ms3.estc0.25.rhoX * t(A3.R3*dDef) ( size 8 x 1000 ) | For each of the 1,000 replicates
```

The R script for these simulations, ready to be run in an HPC as an openMP job, can be found [here](https://github.com/sabifo4/morpho/blob/master/00_data/simulations/R_scripts/03_simulate_char_correlation_cluster.R).
If you are simulating from your PC, you can alternatively use [this R script](https://github.com/sabifo4/morpho/blob/master/00_data/simulations/R_scripts/03_simulate_char_correlation_cluster.R), which does not require openMP.
Nevertheless, simulating $R=1,000$ data sets with $p=10,000$ characters can take several days, so we really 
recommend using the R script that has been written to be run as an openMP job.

