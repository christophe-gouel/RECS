%% Installation

%% Download
% RECS Toolbox zip archives are available at <https://github.com/christophe-gouel/RECS/releases>.
%
% Why is this archive 12 MB? Much of this size is due to an executable for
% Windows. The executable file includes a complete Python distribution
% necessary to parse RECS model files.

%% Dependencies
% * MATLAB R2010a or later.
% * <http://www4.ncsu.edu/~pfackler/compecon/ CompEcon toolbox>. RECS depends on
%   the CompEcon toolbox for many programs (especially with respect to
%   interpolation). Please follow CompEcon installation instructions and do not
%   forget to create the mex files if you want your models solved in a
%   reasonable time.

%% Optional dependencies
% * <http://pages.cs.wisc.edu/~ferris/path.html Path solver for MATLAB>. Path is
%   the reference solver for <MCP.html mixed complementarity problems>. Its
%   installation is highly recommended if difficult complementarity problems
%   need to be solved.
% * MATLAB Optimization Toolbox. The solver |fsolve| can be used to solve both the
%   equilibrium equations and the rational expectations equilibrium.
% * MATLAB Parallel Computing Toolbox. This toolbox allows many RECS programs to
%   be run in parallel to speed-up computation.

%% Installation instructions
% # <https://github.com/christophe-gouel/RECS/releases Download the latest RECS archive> and unzip
% it into a folder, called here |recsfolder| (avoid folder names that include
% spaces, even for parent folders).
% # Install the CompEcon toolbox: (i) <http://www4.ncsu.edu/~pfackler/compecon/
% Download the CompEcon toolbox archive>; (ii) Unzip the archive into a folder,
% called here |compeconfolder|; (iii) Add CompEcon to the MATLAB path:
% |addpath('compeconfolder/CEtools','compeconfolder/CEdemos')|; (iv) Type
% |mexall| in MATLAB prompt to create all CompEcon mex files.
% # (optional) Install other dependencies.
% # Add the RECS folder to the MATLAB path: |addpath('recsfolder')|.
% # On Windows, you are all set. On other architectures, you will have to
% install some Python packages. see instructions below.
% # You can test your installation by running RECS demonstration files by typing
% |recsdemos|. You can also access RECS documentation in MATLAB by typing |doc|.
%
%% Install Python packages (for Linux and Mac, or Windows when installing from source)
%
% Python 2.7.X and the following packages are required:
%
% * <http://www.python.org/download/ Python 2.7.X>. Python is preinstalled on
%   Linux and Mac, but you should check the version number.
% * <http://www.numpy.org/ NumPy>.
% * <http://pyyaml.org/wiki/PyYAML PyYaml>.
% * <http://sympy.org SymPy>, version 0.7.2.
%
% To make the Python programs available to RECS, you have two options:
%
% *Install them in a virtual Python environment.*
%
% By default, RECS looks for a folder PythonVirtualEnv inside the Python
% folder. To do the installation, from RECS folder type in a terminal
%
%  cd Python
%  virtualenv PythonVirtualEnv
%  source PythonVirtualEnv/bin/activate
%  pip install numpy PyYAML sympy==0.7.2
%  deactivate
%
% If your default Python installation is not the version 2.7.X, replace the
% second command by
%
%  virtualenv PythonVirtualEnv -p /usr/bin/python2.7
%
% where |/usr/bin/python2.7| should be replaced by the address of your Python
% 2.7.X interpreter.
%
% *Install them in your default Python installation.*
%
% Recommended only if you do not use Python otherwise. In this case, just type
%
%  pip install numpy PyYAML sympy==0.7.2

%% Installation from source
% If you want to work with the bleeding edge version of RECS, which may be
% unstable, or if you want to contribute to RECS development, you need to
% install RECS from source. The installation requires <http://git-scm.com/ Git>.
%
% When installing from source, all platforms (Linux, Mac, and Windows) require
% <http://www.python.org/download/ Python 2.7.X>, along with
% <http://www.numpy.org/ NumPy>, <http://pyyaml.org/wiki/PyYAML PyYaml>, and
% <http://sympy.org SymPy> version 0.7.2.  See above for instructions.
%
% Optionnaly, under Windows, to be able to generate a binary to run to the
% solver without a Python installation, one can also install
% <http://www.pyinstaller.org/ PyInstaller> and make its folder available in
% Windows Path.
%
% # Download the latest version of RECS from the git repository by typing in a
% command line: |git clone https://github.com/christophe-gouel/RECS.git recs|
% # From RECS folder (|cd recs|), download recs submodules with two commands:
% |git submodule init| and |git submodule update|.
% # Install the CompEcon toolbox: (i) <http://www4.ncsu.edu/~pfackler/compecon/
% Download the CompEcon toolbox archive>; (ii) Unzip the archive into a folder,
% called here |compeconfolder|; (iii) Add CompEcon to the MATLAB path:
% |addpath('compeconfolder/CEtools','compeconfolder/CEdemos')|; (iv) Type
% |mexall| in MATLAB prompt to create all CompEcon mex files.
% # (optional) Install other dependencies.
% # Add the RECS folder to the MATLAB path: |addpath('recsfolder')|.
% # Finalizes RECS installation from source by running in MATLAB the function
% |recsInstall|.
% # You can test your installation by running RECS demonstration files by typing
% |recsdemos|. You can also access RECS documentation in MATLAB by typing |doc|.

