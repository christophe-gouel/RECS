**RECS toolbox Version 0.5-git**

A MATLAB solver for nonlinear, dynamic, stochastic, rational expectations
equilibrium models. RECS stands for "Rational Expectations Complementarity
Solver". This name emphasizes that RECS has been developed specifically to solve
models that include complementarity equations, also known as models with
occasionally binding constraints.

[Christophe Gouel](http://www.christophegouel.com) (<christophe.gouel@grignon.inra.fr>)

**Main page is at**: [www.recs-solver.org](http://www.recs-solver.org).

Download
========

RECS Toolbox zip archives are available at
[code.google.com/p/recs/](http://code.google.com/p/recs/).

Why is this archive 12 MB? Much of this size is due to an executable for
Windows. The executable file includes a complete Python distribution necessary
to parse RECS model files.

Dependencies
============

* MATLAB R2010a or later.
* [CompEcon toolbox](http://www4.ncsu.edu/~pfackler/compecon/). RECS depends on
  the CompEcon toolbox for many programs (especially with respect to
  interpolation). Please follow CompEcon installation instructions and do not
  forget to create the mex files if you want your models solved in a reasonable
  time.

Optional dependencies
---------------------

* [Path solver for MATLAB](http://pages.cs.wisc.edu/~ferris/path.html). Path is
  the reference solver for mixed complementarity problems. Its installation is
  highly recommended if difficult complementarity problems need to be solved.
* MATLAB Optimization Toolbox. The solver fsolve can be used to solve both the
  equilibrium equations and the rational expectations equilibrium.
* MATLAB Parallel Computing Toolbox. This toolbox allows many RECS programs to
  be run in parallel to speed-up computation.
* [Sundials Toolbox](https://computation.llnl.gov/casc/sundials/main.html),
  which provides a compiled Newton-Krylov solver for solving the rational
  expectations equilibrium.
* MATLAB Statistics Toolbox. Useful to simulate models in which shocks follow
  other distributions than the normal.

Installation instructions
=========================

1. [Download the latest RECS archive](http://code.google.com/p/recs/) and unzip
   it into a folder, called here `recsfolder` (avoid folder names that include
   spaces, even for parent folders).
2. Install the CompEcon toolbox:
   1. [Download the CompEcon toolbox archive](http://www4.ncsu.edu/~pfackler/compecon/);
   2. Unzip the archive into a folder, called here `compeconfolder`;
   3. Add CompEcon to the MATLAB path: `addpath('compeconfolder/CEtools','compeconfolder/CEdemos')`;
   4. Type `mexall` in MATLAB prompt to create all CompEcon mex files.
3. (optional) Install other dependencies.
4. Add the RECS folder to the MATLAB path: `addpath('recsfolder')`.
5. On Windows, you are all set. On other architectures, you will have to install
   some Python packages. see instructions below.
6. You can test your installation by running RECS demonstration files by typing
   `recsdemos`. You can also access RECS documentation in MATLAB by typing `doc`.

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

One solution might be to install a scientific Python distribution such as
[EPD](http://www.enthought.com/).

Let me know whether or not it works.

Source
=======

RECS source can be found on the following git repository:
<https://github.com/christophe-gouel/RECS.git>.

License
=======

Unless stated otherwise, all files in the RECS toolbox are licensed using the
Expat license, a permissive free software license. Please see the [software
license](https://raw.github.com/christophe-gouel/RECS/master/LICENSE.txt) for
more information.

