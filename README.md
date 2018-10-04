## Internship Summer 2018

Author: Oriol CHANDRE-VILA.

Advisors: Joseph MORLIER (ISAE-SUPAERO) and Sylvain DUBREUIL (ONERA). 

ONERA, 2018. 

Département Traitement de l'Information et Systèmes (DTIS).

### Synthesis

This internship focuses the use of Proper Orthogonal Decomposition in an aeroelastic problem. Given a code "aero_struct" created by ONERA, the POD approach has been used in order to perform a offline/online processus.
Thanks to this process, super fast simulations can be achieved in online mode, but a robust expensive offline phase is needed.
The code is implemented with Python 2.7.

### Methodology

- OFFLINE PHASE: In this stage, a Greedy Algorithm has been used to create the Reduced Base POD for each parameter of interest; circulation in Aerodynamics and displacement in Structures. In each iteration (total number of iterations is predetermined as a variable), a finite number of points (second predetermined variable) is performed. The goal of the Greedy algorithm is to perform the reduced approach an compute an Error Indicator. Then, there where this EI is the highest, run a complete analysis and add the results to the corresponding RB.
- ONLINE PHASE: Here, a kriging interpolation is used in order to run only a specific case. With the kriging, the values of displacement and circulation are calculated, and then the Von Mises Stress, the error and others characteristics can be computed too.


### In this GitHub...

1. Codes: One can find all the codes needed to perform the study. Just two things to be commented:

    - The Python version used is the 2.7

    - A mesh creator (and viewer) is used: GMSH 3.0.6, please note that it is not the last version available!

2. Files: One can find attached both the EXCEL and the images. The EXCEL file is used to campare the results between the POD and the original code.

3. Presentation: the support for the presentation can be found inside.

4. Report: The document where the internship's results are explained is attached.

### Conclusions

(From the internship report)

First of all, it can easily be proved that a real time process saves a great amount of computation time and efforts once an offline process is performed before. However, it will be interesting to be applied in a problem if only a lot of cases must be simulated. If not, it does not worth to perform such an expensive offline calculus.

Secondly, the parameter that drives the accuracy of the RB is the number of samples used to create it (as it was expected). If more samples are computed, the approximation will use more modes of the problem, so the difference between the approximation and the real solution is smoother. Consequently, the number of candidates increases the computation time and adds more information to the kriging (more points to interpolate). Furthermore, because of the kriging process, they must be more candidates than samples (nCandidates > nSamples).

Finally, the kriging procedure is giving an acceptable variance values, but the real solution is not localized between the limits (in the case of the displacement). This can only be explained by a lag of accuracy of the POD-RB. A source of error could be the construction of the RB. Since it is a coupled problem, the real solution added to the RB could not reduced the error to zero. For sure, the error in that point is reduced, but the coupling of the problem can produce a non-zero error indicator. In the original Greedy algorithm, this problem does not exist because the same point can be added several times.

In future works, the idea of working on the different source of error to achieve more accurate results should be studied by improving the creation of the reduced base and by possibly adding points to the kriging functions.
