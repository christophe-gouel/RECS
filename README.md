RECS toolbox Version 0.7-beta
=============================

A MATLAB solver for nonlinear, dynamic, stochastic, rational expectations
equilibrium models. RECS stands for "Rational Expectations Complementarity
Solver". This name emphasizes that RECS has been developed specifically to solve
models that include complementarity equations, also known as models with
occasionally binding constraints.

[Christophe Gouel](http://www.christophegouel.com) (<christophe.gouel@grignon.inra.fr>)

**Main page is at**: [www.recs-solver.org](http://www.recs-solver.org).

## Download

RECS Toolbox zip archives are available at
[https://github.com/christophe-gouel/RECS/releases](https://github.com/christophe-gouel/RECS/releases).

Why is this archive 12 MB? Much of this size is due to an executable for
Windows. The executable file includes a complete Python distribution necessary
to parse RECS model files.

## Dependencies

* MATLAB R2010a or later.
* [CompEcon toolbox](http://www4.ncsu.edu/~pfackler/compecon/). RECS depends on
  the CompEcon toolbox for many programs (especially with respect to
  interpolation). Please follow CompEcon installation instructions and do not
  forget to create the mex files if you want your models solved in a reasonable
  time.

### Optional dependencies

* [Path solver for MATLAB](http://pages.cs.wisc.edu/~ferris/path.html). Path is
  the reference solver for mixed complementarity problems. Its installation is
  highly recommended if difficult complementarity problems need to be solved.
* MATLAB Optimization Toolbox. The solver fsolve can be used to solve both the
  equilibrium equations and the rational expectations equilibrium.
* MATLAB Parallel Computing Toolbox. This toolbox allows many RECS programs to
  be run in parallel to speed-up computation.

## Installation instructions

1. [Download the latest RECS archive](https://github.com/christophe-gouel/RECS/releases) and unzip
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

### Install on Linux and Mac

Python 2.7.X and the following packages are required and have to be accessible
on your PATH:

* [Python 2.7.X](http://www.python.org/download/). Python is preinstalled on
  Linux and Mac, but is often too old to be useful.
* [PyYaml](http://pyyaml.org/wiki/PyYAML).
* [SymPy](http://sympy.org).
* [SciPy](http://www.scipy.org/Download).

On Debian/Ubuntu, to install the necessary packages type in a terminal:

    sudo apt-get install python-yaml python-sympy python-scipy

One solution might be to install a scientific Python distribution that would
contain all the necessary packages such as
[Anaconda](https://store.continuum.io/cshop/anaconda/).

## Installation from source

If you want to work with the bleeding edge version of RECS, which may be
unstable, or if you want to contribute to RECS development, you need to install
RECS from source. The installation requires [Git](http://git-scm.com/).

When installing from source, all platforms (Linux, Mac, and Windows) require
[Python 2.7.X](http://www.python.org/download/), along with
[PyYaml](http://pyyaml.org/wiki/PyYAML), [SymPy](http://sympy.org), and
[SciPy](http://www.scipy.org/Download). One solution is to install everything at
once by installing a scientific Python distribution such as
[Anaconda](https://store.continuum.io/cshop/anaconda/).

Optionnaly, under Windows, to be able to generate a binary to run to the solver
without a Python installation, one can also install
[PyInstaller](http://www.pyinstaller.org/) and make its folder available in
Windows Path.

1. Download the latest version of RECS from the git repository by typing in a
   command line: `git clone https://github.com/christophe-gouel/RECS.git recs`
2. From RECS folder (`cd recs`), download recs submodules with two commands:
   `git submodule init` and `git submodule update`.
3. Install the CompEcon toolbox:
   1. [Download the CompEcon toolbox archive](http://www4.ncsu.edu/~pfackler/compecon/);
   2. Unzip the archive into a folder, called here `compeconfolder`;
   3. Add CompEcon to the MATLAB path: `addpath('compeconfolder/CEtools','compeconfolder/CEdemos')`;
   4. Type `mexall` in MATLAB prompt to create all CompEcon mex files.
4. (optional) Install other dependencies.
5. Add the RECS folder to the MATLAB path: `addpath('recsfolder')`.
6. Finalizes RECS installation from source by running in MATLAB the function
   `recsInstall`.
7. You can test your installation by running RECS demonstration files by typing
   `recsdemos`. You can also access RECS documentation in MATLAB by typing `doc`.

## Source

RECS source can be found on the following git repository:
<https://github.com/christophe-gouel/RECS.git>.

## License

Unless stated otherwise, all files in the RECS toolbox are licensed using the
Expat license, a permissive free software license. Please see the [software
license](https://raw.github.com/christophe-gouel/RECS/master/LICENSE.txt) for
more information.

