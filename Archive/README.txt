**RECS**

A Matlab solver for nonlinear, dynamic, stochastic, rational expectations
equilibrium models.

[Christophe Gouel](http://www.christophegouel.com) (<christophe.gouel@grignon.inra.fr>)

INSTALL
=======

To use the solver you have to

1. Extract the archive to a folder (avoid folder name with space, even for
   parent folders).
2. Install the CompEcon toolbox (see Dependencies), follow CompEcon installation
   instructions and do not forget to create the mex files.
3. Add RECS folder to the Matlab path.

When the installation is finished, you can test it by running `recsdemos.m`, which
launches all demonstration files. RECS documentation should also be available in
Matlab when you type `doc`.

DEPENDENCIES
============

In addition to Matlab (from 7.9 (R2009b) to 7.13 (R2011b), and likely beyond),
RECS programs depend on the following packages, which have to be installed:

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
that calculates automatically the Jacobians (used in recsmodelinit.m) has been
written in Python. Executables have been prepared for Windows (32- and
64-bit). For other platforms, you can use dolo-recs from its Python sources:
<https://github.com/albop/dolo>. Send me an email if you need help on this!

Sources
=======

RECS sources can be found on the following git repository:
`https://github.com/christophe-gouel/RECS.git`

