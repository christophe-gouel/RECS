**RECS**

A Matlab solver for nonlinear, dynamic, stochastic, rational expectations
equilibrium models.

[Christophe Gouel](http://www.christophegouel.com) (<christophe.gouel@grignon.inra.fr>)

INSTALL
=======

If you use git, you can get the latest version with

`git clone https://github.com/christophe-gouel/RECS.git recs`

Without git, you can download an archive of the latest version by clicking on the
**Downloads** button on [RECS webpage](https://github.com/christophe-gouel/RECS).

To use the solver you have to

1. Install the CompEcon toolbox (see Dependencies), follow CompEcon installation
   instructions and do not forget to create the mex files.
2. Add the solver folder (avoid folder name with space, even for parent folders)
   to the Matlab path.
3. Run `recsInstall.m` to finalize RECS installation.

When the installation is finished, you can test it by running `recsdemos.m`, which
launches all demonstration files.

DEPENDENCIES
============

In addition to Matlab (from 7.9 (R2009b) to 7.13 (R2011b)), RECS programs depend
on the following packages, which have to be installed:

* CompEcon (<http://www4.ncsu.edu/~pfackler/compecon/>).

OPTIONAL INSTALLATION
---------------------

* Path solver for Matlab (<http://pages.cs.wisc.edu/~ferris/path.html>). This is
  the reference solver for mixed complementarity problems. It is highly recommended
  to install if difficult complementarity problems need to be solved.
* Matlab Optimization Toolbox. The solver `fsolve` can be used to solve both the
  equilibrium equations and the rational expectations equilibrium.
* Sundials Toolbox (<https://computation.llnl.gov/casc/sundials/main.html>), which
  provides a compiled Newton-Krylov solver for solving the rational expectations
  equilibrium.
* Matlab Statistics Toolbox. Useful to simulate models in which shocks follow
  distribution other than normal.

Supported platforms
===================

RECS Matlab files work on any platform supported by Matlab. However, the parser
that calculates automatically the Jacobians has been written in
Python. Executables have been prepared for Windows (32- and 64-bit) and Linux
(32-bit). For other platforms, you can use dolo-recs from its Python sources:
<https://github.com/albop/dynare-python>. Send me an email if you need help on
this!
