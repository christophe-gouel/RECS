**RECS**

A solver for nonlinear, dynamic, stochastic, rational expectations equilibrium
models.

Christophe Gouel <christophe.gouel@grignon.inra.fr>

To use the solver you have to

1. Install the CompEcon toolbox (see Dependencies), see CompEcon installation
   instructions and do not forget to create the mex files.
2. Add the solver folder to the Matlab path.
3. Run `recsInstall.m` to copy a mex file from CompEcon to the solver folder.

When the installation is finished, you can test it by running `recsdemos.m`, which
launches all demonstration files.

DEPENDENCIES
============

RECS programs depend on the following packages, which have to be installed:

* CompEcon (http://www4.ncsu.edu/~pfackler/compecon/).

OPTIONAL INSTALLATION
---------------------

* Path solver for Matlab (http://www.cs.wisc.edu/cpnet/cpnetsoftware/). This is
  the reference solver for mixed complementarity problems. It is highly recommended
  to install if difficult complementarity problems need to be solved.
* Matlab Statistics Toolbox. Useful to simulate models in which shocks follow
  distribution other than normal.
* Sundials Toolbox (https://computation.llnl.gov/casc/sundials/main.html), which
  provides a compiled Newton-Krylov solver for solving the rational expectations
  equilibrium.
* Matlab Optimization Toolbox. The solver fsolve can be used to solve both the
  equilibrium equations and the rational expectations equilibrium.
