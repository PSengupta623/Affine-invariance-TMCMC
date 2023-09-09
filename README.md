# Affine-invariance-TMCMC
Two-stage Bayesian model updating_Stage-1-Model reduction_Stage-2-Affine-invariance-MDS factor accounting eigenfrequencies and mode shapes of different nature and dimensions-Tuning algorithm MDS based
Affine-invariance-based Transitional Markov Chain Monte Carlo (TMCMC) with the improved iterated model reduction technique using modal data
An affine-invariance sampling approach is proposed to estimate a multi-dimensional scaling (MDS) factor in the affine-transformed space that accounts for the eigenfrequencies and mode shapes which are of different nature and dimensions at the levels of likelihood, prior and posterior distributions. The TMCMC sampling adaptively utilizes the plausibility value from the proposed MDS-based tuning algorithm to improve its transition levels. 
This repository presents codes for the proposed TMCMC sampler (written in MATLAB) which seeks to demonstrate the implementation of the same along with its three illustrative examples. The details are as follows:
## 1) Illustrative Example: A simply supported steel beam
See: SS_steel_beam.m
An updating parameter theta (with i=1, 2….,10) is enforced on each stiffness parameter. Each element has a variable percentage of stiffness reductions relative to the nominal stiffness value to represent damage scenarios. The stiffness reduction percentages for various elements 
The modal data obtained from the nominal model is taken as the mean with the coefficient of variation (COV) of 0.1, 0.15 and 0.2 as Gaussian noise for the simulation of the modal responses. 
The measurements for corresponding DOFs are gradually decreased.
## 2) Example A three-dimensional ten-storied building:
See: 3D_FEA_building.jnl
A ten-storied RC building with the modelling details along with geometrical and material properties of the slab, beam, and column elements along with the contact elements between the element interfaces are provided in Sengupta and Chakraborty (2023) (https://doi.org/10.1016/j.ymssp.2022.109586). The total number of DOFs is 39,136. 
The proposed two-stage model updating algorithm is applied after subdividing the model into several substructures. Three conditions of substructuring are considered for the demonstration of the proposed model updating algorithm. Now, the varying stiffness reductions are considered. From the first to the fourth levels, the nominal stiffness of columns is reduced by 30 per cent and from the fifth to the tenth floors by 20 per cent. The slabs and beams from the first to the fourth levels are defined with a stiffness reduction of 30 per cent and from the fifth to the tenth floors, a reduction of 25 per cent is considered. 
The posterior mean values and COVs of the predicted parameters for different substructure conditions under different measurement configurations are obtained at 10 per cent, 15 per cent and 20 per cent noise levels with Nm=4.
## 3) Example Experimental data with a spring-mass model:
See: experimental_springmass_model.m
The system consists of 76.2 mm diameter and 25.4 mm thick aluminium discs with holes at the centre. A steel rod is placed in the holes for smooth sliding. Steel collars on each side of the discs are used to connect the collar springs. The springs have a constant stiffness of 56.7 kN/m and a mass of 559.3 gm, except the mass at the far left which is 419.4 gm. 
The mass-spring system is vibrated by an electrodynamic shaker. The acceleration measurements of each mass are taken. 

For more details on the experimental set-up, readers can also refer to the work by [Lyngdoh et al. (2019)] (https://doi.org/10.1061/(ASCE)EM.1943-7889.0001668).
# Reference(s):
* P. Sengupta, S. Chakraborty, An improved iterative model reduction technique to estimate the unknown responses using limited available responses, *Mechanical Systems and Signal Processing, 182 (2023) 109586. https://doi.org/10.1016/j.ymssp.2022.109586
* P. Sengupta, S. Chakraborty, An improved Bayesian model updating framework by enhanced iterative model reduction technique in the time domain, *Journal of Sound and Vibration, 549 (2023) 117589, https://doi.org/10.1016/j.jsv.2023.117589
* J. Coullon, R.J. Webber, Ensemble sampler for infinite-dimensional inverse problems, Stat. Comput. 31(3) (2021) 1-9. https://doi.org/10.1007/s11222-021-10004-y
* J.A. Vrugt, Markov chain Monte Carlo simulation using the DREAM software package: theory, concepts, and MATLAB implementation, *Environmental Modelling & Software, 75 (2016) 273–316. https://doi.org/10.1016/j.envsoft.2015.08.013
* A. Lye, A. Cicirello, E. Patelli, An efficient and robust sampler for Bayesian inference: Transitional Ensemble Markov Chain Monte Carlo. *Mechanical Systems and Signal Processing, 167 (2022) 108471. https://doi.org/10.1016/j.ymssp.2021.108471
* A. Lye, A. Cicirello, E. Patelli, Sampling methods for solving Bayesian model updating problems: A tutorial. *Mechanical Systems and Signal Processing, 159 (2021) 107760. https://doi.org/10.1016/j.ymssp.2021.107760
* D. Foreman-Mackey, D.W. Hogg, D. Lang, J. Goodman, emcee v3: A Python ensemble sampling toolkit for affine-invariant MCMC, arXiv preprint arXiv:1911.07688 (2019).
* J. Coullon, R.J. Webber, Ensemble sampler for infinite-dimensional inverse problems, Statistics and Computing 31(3) (2021) 1-9. https://doi.org/10.1007/s11222-021-10004-y
* R. Rocchetta, M. Broggi, Q. Huchet, E. Patelli, On-line Bayesian model updating for structural health monitoring. *Mechanical Systems and Signal Processing, 103 (2018) 174-195. https://doi.org/10.1016/j.ymssp.2017.10.015
* Z. Yuan, P. Liang, T. Silva, K. Yu, J. E. Mottershead, Parameter selection for model updating with global sensitivity analysis. *Mechanical Systems and Signal Processing, 115 92019) 483-496. https://doi.org/10.1016/j.ymssp.2018.05.048
* S. Bansal, Bayesian Model updating using modal data based on dynamic condensation, J. Eng. Mech. 146(2) (2020) 04019123. https://doi.org/10.1061/(ASCE)EM.1943-7889.0001714.

# Author:
* Name: Partha Sengupta
* Contact: senguptapartha91@gmail.com
* Affiliation: Indian Institute of Engineering Science and Technology Shibpur
