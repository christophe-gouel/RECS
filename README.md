**RECS**

A Matlab solver for nonlinear, dynamic, stochastic, rational expectations
equilibrium models. RECS stands for Rational Expectations Complementarity
Solver, which emphasizes the fact that RECS is specifically developed to solve
models that include complementarity equations, also known as models with
occasionally binding constraints.

[Christophe Gouel](http://www.christophegouel.com) (<christophe.gouel@grignon.inra.fr>)

Download
========

RECS Toolbox zip archives are available at
<https://github.com/christophe-gouel/RECS/downloads>.

Why is this archive as large as 12 MB? Most of it is due to an executable for
Windows. This executable includes a complete Python distribution necessary to
parse RECS model files.

Dependencies
============

* Matlab R2009b or later.
* [CompEcon toolbox](http://www4.ncsu.edu/~pfackler/compecon/). RECS depends on
  the CompEcon toolbox for many programs (especially with respect to
  interpolation). Please follow CompEcon installation instructions; in
  particular, do not forget to create the mex files if you want your models
  solved in a reasonable time.

Optional dependencies
---------------------

* [Path solver for Matlab](http://pages.cs.wisc.edu/~ferris/path.html). This is
  the reference solver for mixed complementarity problems. It is highly
  recommended to install if difficult complementarity problems need to be
  solved.
* MATLAB Optimization Toolbox. The solver fsolve can be used to solve both the
  equilibrium equations and the rational expectations equilibrium.
* [Sundials Toolbox](https://computation.llnl.gov/casc/sundials/main.html),
  which provides a compiled Newton-Krylov solver for solving the rational
  expectations equilibrium.
* Matlab Statistics Toolbox. Useful to simulate models in which shocks follow
  distributions other than normal.

Installation instructions
=========================

1. [Download the latest RECS
   archive](https://github.com/christophe-gouel/RECS/downloads) and unzip it
   into some folder (avoid folder name with space, even for parent folders),
   called here `recsfolder`.
2. Install the CompEcon toolbox: 
   1. [Download CompEcon toolbox archive](http://www4.ncsu.edu/~pfackler/compecon/); 
   2. Unzip the archive into some folder, called here `compeconfolder`; 
   3. Add CompEcon to the Matlab path: `addpath('compeconfolder/CEtools','compeconfolder/CEdemos')`; 
   4. Type `mexall` in Matlab prompt to create all CompEcon mex files.
3. (optional) Install other dependencies.
4. Add RECS folder to the Matlab path: `addpath('recsfolder')`.
5. On Windows, you are all set. On other architectures, you have to to install
   some Python packages. see instructions below.
6. You can test your installation by running RECS demonstration files by typing
   `recsdemos`. You can also access RECS documentation in Matlab by typing `doc`.

Install on Linux
----------------

Python 2.7.X is required. On Debian/Ubuntu, to install the necessary packages
type in a terminal:
 
    sudo apt-get install python-yaml python-sympy python-scipy

Install on Mac
--------------

In this case, you are on your own. You have to install 

* [Python 2.7.X](http://www.python.org/download/). Python is preinstalled on
  Mac, but is usually too old to be useful.
* [PyYaml](http://pyyaml.org/wiki/PyYAML).
* [SymPy](http://sympy.org).
* [SciPy](http://www.scipy.org/Download).

One solution could be to install a scientific Python distribution such as
[EPD](http://www.enthought.com/).

Let me know if it works or not.

Source
=======

RECS source can be found on the following git repository:
<https://github.com/christophe-gouel/RECS.git>.

License
=======

Unless stated otherwise, all files in the RECS toolbox are licensed using the
Expat license, a permissive free software license. Please see the [software
license](https://github.com/christophe-gouel/RECS/blob/master/LICENSE.txt) for
more information.

