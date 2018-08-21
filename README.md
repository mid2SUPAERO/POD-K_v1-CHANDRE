## Stage Summer 2018

Author: Oriol CHANDRE-VILA.

Advisors: Joseph MORLIER (ISAE-SUPAERO) and Sylvain DUBREUIL (ONERA). 

ONERA, 2018. 

Département Traitement de l'Information et Systèmes (DTIS).

### Synthese

This stage focuses the use of Proper Orthogonal Decomposition in an aeroelastic problem. Given a code "aero_struct" created by ONERA, the POD approach has been used in order to perform a offline/online processus.
Thanks to this process, super fast simulations can be achieved in online mode, but a robust expensive offline phase is needed.
The code is implemented with Python 2.7.

### Methodology

- OFFLINE PHASE: In this stage, a Greedy Algorithm has been used to create the Reduced Base POD for each parameter of interest; circulation in Aerodynamics and displacement in Structures. In each iteration (total number of iterations is predetermined as a variable), a finite number of points (second predetermined variable) is performed. The goal of the Greedy algorithm is to perform the reduced approach an compute an Error Indicator. Then, there where this EI is the highest, run a complete analysis and add the results to the corresponding RB.
- ONLINE PHASE: Here, a kriging interpolation is used in order to run only a specific case. With the kriging, the values of displacement and circulation are calculated, and then the Von Mises Stress, the error and others characteristics can be computed too.


### In this GitHub...

1. Codes: all the live codes can be found inside.

MATLAB tutorial on [PGD for Poisson Eq. 2D](http://htmlpreview.github.io/?https://github.com/mid2SUPAERO/PIR_CHANDRE_ROM/blob/master/Codes/html_Poisson2D_PGD/main.html)
  
    
MATLAB tutorial on [POD for Steady Advection-Diffusion Problem](http://htmlpreview.github.io/?https://github.com/mid2SUPAERO/PIR_CHANDRE_ROM/blob/master/Codes/html_AdvDiff_POD/AdvectionDiffusion.html)
    
    
MATLAB tutorial on [PGD for Steady Advection-Diffusion Problem](http://htmlpreview.github.io/?https://github.com/mid2SUPAERO/PIR_CHANDRE_ROM/blob/master/Codes/html_AdvDiff_PGD_Online/Online.html) 
    
    But of course before need to do the offline learning
    
MATLAB tutorial on [PGD for Steady Advection-Diffusion Problem](http://htmlpreview.github.io/?https://github.com/mid2SUPAERO/PIR_CHANDRE_ROM/blob/master/Codes/html_AdvDiff_PGD_Offline/Offline.html) 
    
2. Presentation: the support for the presentation can be found inside.

3. Report: the editable and the PDF file of the written report can be found inside.
