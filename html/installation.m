%% Installation

%% Download
% RECS Toolbox zip archives are available at
% <https://github.com/christophe-gouel/RECS/downloads>.
%
% Why is this archive as large as 12 MB? Most of it is due to an executable for
% Windows. This executable includes a complete Python distribution necessary to
% parse RECS model files.

%% Dependencies
% * Matlab R2009b or later.
% * <http://www4.ncsu.edu/~pfackler/compecon/ CompEcon toolbox>. RECS depends on
%   the CompEcon toolbox for many programs (especially with respect to
%   interpolation). Please follow CompEcon installation instructions; in
%   particular, do not forget to create the mex files if you want your models
%   solved in a reasonable time.

%% Optional dependencies
% * <http://pages.cs.wisc.edu/~ferris/path.html Path solver for Matlab>. This is
%   the reference solver for <MCP.html mixed complementarity problems>. It is
%   highly recommended to install if difficult complementarity problems need to
%   be solved.
% * MATLAB Optimization Toolbox. The solver fsolve can be used to solve both the
%   equilibrium equations and the rational expectations equilibrium.
% * <https://computation.llnl.gov/casc/sundials/main.html Sundials Toolbox>,
%   which provides a compiled Newton-Krylov solver for solving the rational
%   expectations equilibrium.
% * Matlab Statistics Toolbox. Useful to simulate models in which shocks follow
%   distributions other than normal.

%% Installation instructions
% # <https://github.com/christophe-gouel/RECS/downloads Download the latest RECS
% archive> and unzip it into some folder (avoid folder name with space, even for
% parent folders), called here |recsfolder|.
% # Install the CompEcon toolbox: (i) <http://www4.ncsu.edu/~pfackler/compecon/
% Download CompEcon toolbox archive>; (ii) Unzip the archive into some folder,
% called here |compeconfolder|; (iii) Add CompEcon to the Matlab path:
% |addpath('compeconfolder/CEtools','compeconfolder/CEdemos')|; (iv) Type
% |mexall| in Matlab prompt to create all CompEcon mex files.
% # (optional) Install other dependencies.
% # Add RECS folder to the Matlab path: |addpath('recsfolder')|.
% # On Windows, you are all set. On other architectures, you have to to install
% some Python packages. see instructions below.
% # You can test your installation by running RECS demonstration files by typing
% |recsdemos|. You can also access RECS documentation in Matlab by typing |doc|.
%
% *Install on Linux*
%
% On Debian/Ubuntu, type in a terminal:
% 
%  sudo apt-get install python-setuptools python-yaml python-scipy python-sympy python-matplotlib
%
% *Install on Mac*
%
% In a terminal type:
%
%  sudo port install py27-setuptools
%  sudo port install py27-yaml
%  sudo port install py27-scipy
%  sudo port install py27-sympy
%  sudo port install py27-matplotlib
