%% Installation

%% Download
% RECS Toolbox zip archives are available at <http://code.google.com/p/recs/>.
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
% # <http://code.google.com/p/recs/ Download the latest RECS archive> and unzip
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
% *Install on Linux*
%
% Python 2.7.X is required. On Debian/Ubuntu, to install the necessary packages
% type in a terminal:
%
%  sudo apt-get install python-yaml python-sympy python-scipy
%
% *Install on Mac*
%
% In this case, you are on your own. You have to install
%
% * <http://www.python.org/download/ Python 2.7.X>. Python is preinstalled on
% Mac, but is usually too old to be useful.
% * <http://pyyaml.org/wiki/PyYAML PyYaml>.
% * <http://sympy.org SymPy>.
% * <http://www.scipy.org/Download SciPy>.
%
%
% One solution might be to install a scientific Python distribution such as
% <http://www.enthought.com/ EPD>.
%
% Let me know whether or not it works.

%% Installation from source
% If you want to work with the bleeding edge version of RECS, which may be
% unstable, or if you want to contribute to RECS development, you need to
% install RECS from source. The installation requires <http://git-scm.com/ Git>.
%
% When installing from source, all platforms (Linux, Mac, and Windows) require
% <http://www.python.org/download/ Python 2.7.X>, along with
% <http://pyyaml.org/wiki/PyYAML PyYaml>, <http://sympy.org SymPy>, and
% <http://www.scipy.org/Download SciPy>. Under Windows, it is also necessary to
% install <http://www.pyinstaller.org/ PyInstaller> and to make its folder
% available in Windows Path.
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

