# README #

Copyright (C) 2025 Oscar Bouverot-Dupuis <oscar.bouverot-dupuis@ipht.fr>, Alberto Rosso <alberto.rosso@cnrs.fr> and Laura Foini <laura.foini@ipht.fr>.

This program is free software; you can redistribute it and/or modify it under the terms of either the CeCILL license or the GNU General Public license, as included with the software package.

### Summary of the numerical project ###

This is the public version of the code presented in "The generic Mott transition in the sine-Gordon model through an embedded worm algorithm", Oscar Bouverot-Dupuis, Laura Foini and Alberto Rosso.

See the arXiv preprint:

This package performs Monte Carlo simulations on the tilted sine-Gordon model, using C++. It incorporates the Worm algorithm (Wo) and the smooth worm algorithm (SmoWo).

### Organization of the code ###
- Constants.h is the file where the user declares all the simulation parameters. It is the only file the user needs to edit.

- main.cpp and Imports.h respectively contain the common sampling routine for all algorithms and necessary package imports.

- PhysSystem\_\*.h and Sampler\_\*.h (\* = Wo, SmoWo) contain the algorithm-specific parts of the code. Sampler\_\*.h performs a high-level sampling routine by calling lower-level functions grouped in PhysSystem\_\*.h.


### Set up ###

The version of C++ used was C++17.
Choose the simulation parameters by editing the Constants.h file. The parameters to update from are
- #include "Sampler\_\*.h" : choose algorithm to use by replacing \* by Wo or SmoWo.
- Ls : all system sizes $L$. Format : vector of strictly positive ints.
- Bs : all system inverse temperatures $\beta=1/T$. Format : vector of strictly positive ints.
- Ks : all Luttinger parameters $K$. Format : vector of strictly positive doubles. 
- mus : all chemical potentials $\mu$. Format : vector of positive doubles.
- g : Umklapp strength $g$. Format : positive double.
- lambda_w : worm event rate $\lambda_\mathrm{w}$. Format : strictly postive double.
- lambda_r_prefactor : sets the refreshment rate as $\lambda_\mathrm{r}=$ lambda_r_prefactor $/(\beta L)$. Format : strictly postive double.
- path : path to store the data. Format : string.
- OBSERVABLES : observables (can do one or several) to compute and save. "N_space" and "N_time" return the winding numbers $N_x$ and $N_\tau$, "kappa" and "rho_s" are the compressibility $\kappa$ and superfluid stiffness $\rho_s$, "C_2kf" is the amplitude of $2k_F$ modulation of the density, "C_theta" is an array containing the times the worm head has spent at a distance $x$ and time $\tau$ from its tail, "Ct_phi" and "Cx_phi" are $C_\varphi(x,0)$ and $C_\varphi(0,\tau)$, "algotime" is the algorithmic time in number of operations (not sweeps) and "field" the actual field configuration. Format : vector of strings.
- PERCENT\_SAVE : save observables after every PERCENT\_SAVE % of the simulation. Format : strictly positive int.
- TOTAL\_NUMBER\_SAMPLE : total number of samples for each observables. Format : strictly positive double.
- OUTPUT\_DISTANCE\_SWEEPS : samples are ouputted after OUTPUT\_DISTANCE\_SWEEPS*N^2 operations. Format : strictly positive int.
- data_block_avg : the array-formatted observables are averaged over data_blcok_avg samples before being saved. This prevents the output files from being too large. Format: strictly positive int.

### Usage ###
Update the Constants.h file as explained in the above. Execute the file main.cpp using the g++ compiler command or a C++ editor (Visual Studio, etc.).


### Funding ###
This work was made possible thanks to the support of the French ANR under the grant ANR-22-CMAS-0001 (_QuanTEdu-France_ project), and the support of the French government through the France 2030 program (PhOM – Graduate School of Physics), under reference ANR-11-IDEX-0003 (Project Mascotte, L. Foini).
