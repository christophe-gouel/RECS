%% Installation

%% How to get RECS latest version?
% * With <http://git-scm.com/ git>, you can get the latest version with
%   |git clone https://github.com/christophe-gouel/RECS.git recs|.
% * Without git, you can download an archive of the latest version by clicking on
%   the *Downloads* button on
%   <matlab:web('https://github.com/christophe-gouel/RECS','-browser') RECS webpage>.

%% Installation instructions
% # Install the CompEcon toolbox.
% # (optional) Install other dependencies
% # Add RECS folder (avoid folder name with space, even for parent folders) to
%   the Matlab path.
% # Run |recsInstall| to complete RECS installation. It takes 2-3 minutes,
%   because it creates RECS documentation by simulating all demonstration files.

%% CompEcon toolbox
% RECS depends on the CompEcon toolbox for many programs (especially with respect
% to interpolation and multiplication of multidimensional arrays). CompEcon is
% available at <http://www4.ncsu.edu/~pfackler/compecon/>. Please follow CompEcon
% installation instructions; in particular, do not forget to create the mex files!

%% Optional dependencies
% * Path solver for Matlab (<http://pages.cs.wisc.edu/~ferris/path.html>). This is
%   the reference solver for mixed complementarity problems. It is highly
%   recommended to install if difficult complementarity problems need to be solved.
% * MATLAB Optimization Toolbox. The solver fsolve can be used to solve both the
%   equilibrium equations and the rational expectations equilibrium.
% * Sundials Toolbox (<https://computation.llnl.gov/casc/sundials/main.html>),
%   which provides a compiled Newton-Krylov solver for solving the rational
%   expectations equilibrium.
% * Matlab Statistics Toolbox. Useful to simulate models in which shocks follow
%   distribution other than normal.

%% Supported platforms
% RECS Matlab files work on any platform supported by Matlab. However, the
% parser that calculates automatically the Jacobians has been written in
% Python. Executables have been prepared for Windows (32- and 64-bit) and Linux
% (32-bit). For other platforms, you can use dolo-recs from its Python sources:
% <https://github.com/albop/dynare-python>. Send me an email if you need help on
% this!
