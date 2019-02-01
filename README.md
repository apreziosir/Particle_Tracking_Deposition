# Numerical Particle Tracking Model to Estimate Fine Sediments' Deposition in an Idealized Sand Dune

This folder provides the code necessary to run a numerical particle-tracking model in an homogeneous sandbed with periodical boundary conditions and perform a basic results popost processing. Moreover, this folder is the Supplementary Information provided for the paper "Fine Sediment Deposition and Filtration Under Losing and Gaining Conditions: A Particle-Tracking Model Approach" 
(this information has to be completed when the paper gets accepted by the journal)

The code is based in the numerical Particle Tracking model by Li et al. 2017 (1) that estimates nitrogen uptake in arctic rivers.
to estimate fine particles' deposition in an idealized homogeneous bed with steady free surface flow conditions, similar to the ones proposed by Elliott and Brooks in 1997 (2) and Packman et al. 2000 (3). 

## Code readme and instructions

### Particle tracking code
The main part of the numerical particle tracking model is located in the script PTFunct.m. This function runs in Octave and Matlab and its inputs are the groundwater velocities and the filtration coefficient of the bed. Other physical parameters can be changed inside the code. 

The aforementioned code generates a .mat file that sotres particles' coordinates and various parameters that can be found in the final part of the PTFunct.m script. 

### Counting and grouping particles




## References
(1) Li, A., Aubeneau, A. F., Bolster, D., Tank, J. L., & Packman, A. I.(2017, aug). Covariation in patterns of turbulence-driven hyporheic flow and denitrification enhances reach-scale nitrogen removal.Water Resources Research,53(8),6927–6944. doi: 10.1002/2016WR019949
(2) Elliott, A. H., & Brooks, N. H.(1997a, jan).Transfer of nonsorbing solutes to a streambed with bed forms: Laboratory experiments. Water Resources Research,33(1), 137–151.  doi:  10.1029/96WR02783
(3) Packman, A. I., Brooks, N. H., & Morgan, J. J. (2000). A physicochemical model for colloid exchange between a stream and a sand streambed with bed forms. Water Resources Research,36(8), 2351. doi: 10.1029/2000WR900059
