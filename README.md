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

* MATLAB R2011b or later.
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

Python 2.7.X and the following packages are required:

* [Python 2.7.X](http://www.python.org/download/). Python is preinstalled on
  Linux and Mac, but you should check the version number.
* [NumPy](http://www.numpy.org/).
* [PyYaml](http://pyyaml.org/wiki/PyYAML).
* [SymPy](http://sympy.org), version 0.7.2.

To make the Python programs available to RECS, you have two options:

*   Install them in a virtual Python environment. By default, RECS looks for a
    folder PythonVirtualEnv inside the Python folder. To do the installation,
    from RECS folder type in a terminal

        cd Python
        virtualenv PythonVirtualEnv
        source PythonVirtualEnv/bin/activate
        pip install numpy PyYAML sympy==0.7.2
        deactivate

    If your default Python installation is not the version 2.7.X, replace the
	second command by

        virtualenv PythonVirtualEnv -p /usr/bin/python2.7

    where `/usr/bin/python2.7` should be replaced by the address of your Python
    2.7.X interpreter.

* Install them in your default Python installation (recommended only if you do
  not use Python otherwise). In this case, just type

        pip install numpy PyYAML sympy==0.7.2

## Installation from source

If you want to work with the bleeding edge version of RECS, which may be
unstable, or if you want to contribute to RECS development, you need to install
RECS from source. The installation requires [Git](http://git-scm.com/).

When installing from source, all platforms (Linux, Mac, and Windows) require
[Python 2.7.X](http://www.python.org/download/), along with
[NumPy](http://www.numpy.org/), [PyYaml](http://pyyaml.org/wiki/PyYAML), and
[SymPy](http://sympy.org) version 0.7.2. See above for instructions.

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

