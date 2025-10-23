# README #

Copyright (C) 2025 Oscar Bouverot-Dupuis <oscar.bouverot-dupuis@ipht.fr>, Alberto Rosso <alberto.rosso@cnrs.fr> and Manon Michel <manon.michel@cnrs.fr>.

This program is free software; you can redistribute it and/or modify it under the terms of either the CeCILL license or the GNU General Public license, as included with the software package.

### Summary of the numerical project ###

This is the public version of the code presented in "Bosonized 1D quantum systems through enhanced Event-Chain Monte Carlo", Oscar Bouverot-Dupuis, Alberto Rosso, and Manon Michel. 

See the arXiv preprint: https://arxiv.org/abs/2503.11577.

This package performs Monte Carlo simulations on the bosonized dissipative XXZ spin chain, using C++. It incorporates the Metropolis algorithm (Met), the Metropolis algorithm enhanced with cluster moves (Clu-Met), the Event-Chain Monte Carlo algorithm (ECMC), and the ECMC enhanced with cluster moves (Clu-EC).

### Organization of the code ###
- Constants.h is the file where the user declares all the simulation parameters. It is the only file the user needs to edit.

- main.cpp and Imports.h respectively contain the common sampling routine for all algorithms and necessary package imports.

- PhysSystem\_\*.h and Sampler\_\*.h (\* = Met, Clu\_Met, ECMC, Clu\_EC) contain the algorithm-specific parts of the code. Sampler\_\*.h performs a high-level sampling routine by calling lower-level functions grouped in PhysSystem\_\*.h.


### Set up ###

The version of C++ used was C++17.
Choose the simulation parameters by editing the Constants.h file. The parameters to update from are
- #include "Sampler\_\*.h" : choose algorithm to use by replacing \* by Met, Clu\_Met, ECMC, or Clu\_EC.
- Ns : all system sizes N (the 2D system is of size NxN). Format : vector of strictly positive ints.
- Ks : all Luttinger parameters K. Format : vector of strictly positive doubles. 
- alphas : all dissipation strengths alpha. Format : vector of positive doubles.
- g : Umklapp strength g. Format : positive double.
- s : bath exponent s. Format : positive double.
- path : path to store the data. Format : string.
- OBSERVABLES : observables (can do one or several) to compute and save. "orderparameter" returns the staggered magnetization, "algotime" the algorithmic time in number of operations (not sweeps) and "field" the actual field configuration. Format : vector of strings.
- PERCENT\_SAVE : save observables after every PERCENT\_SAVE % of the simulation. Format : strictly positive int.
- TOTAL\_NUMBER\_SAMPLE : total number of samples for each observables. Format : strictly positive double.
- OUTPUT\_DISTANCE\_SWEEPS : samples are ouputted after OUTPUT\_DISTANCE\_SWEEPS*N^2 operations. Format : strictly positive int.


### Usage ###
Update the Constants.h file as explained in the above. Execute the file main.cpp using the g++ compiler command or a C++ editor (Visual Studio, etc.).


### Funding ###

This work was made possible thanks to Institut Pascal at Université Paris-Saclay with the support of the program _Investissements d’avenir_ ANR-11-IDEX-0003-01 and thanks to the support of the French ANR under the grant ANR-22-CMAS-0001 (_QuanTEdu-France_ project), the grant ANR-20-CE46-0007 (_SuSa_ project) and the grant ANR-23-CE30-0031-04 (_DISCREEP_ project).
