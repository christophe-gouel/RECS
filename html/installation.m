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
% # Add RECS folder to the Matlab path.
% # Run |recsInstall| to copy a mex file from CompEcon to the solver folder.
%
% When the installation is finished, you can test it by running |recsdemos|, which
% launches all demonstration files.

%% CompEcon toolbox
% RECS depends on the CompEcon toolbox for many programs (especially with respect
% to interpolation and multiplication of multidimensional arrays). CompEcon is
% available at <http://www4.ncsu.edu/~pfackler/compecon/>. Please follow CompEcon
% installation instructions; in particular, do not forget to create the mex files!

%% Optional dependencies
% * Path solver for Matlab (<http://www.cs.wisc.edu/cpnet/cpnetsoftware/>). This is
%   the reference solver for mixed complementarity problems. It is highly
%   recommended to install if difficult complementarity problems need to be solved.
% * MATLAB Optimization Toolbox. The solver fsolve can be used to solve both the
%   equilibrium equations and the rational expectations equilibrium.
% * Sundials Toolbox (<https://computation.llnl.gov/casc/sundials/main.html>),
%   which provides a compiled Newton-Krylov solver for solving the rational
%   expectations equilibrium.
% * Matlab Statistics Toolbox. Useful to simulate models in which shocks follow
%   distribution other than normal.

%% Issues with 64-bit architecture
% * CompEcon programs cannot be compiled on 64-bit
%   architecture. You can aks me for how to modify them to achieve
%   compilation or you can use uncompiled files, which may, however, lead to very slow
%   programs.
% * RECS calls Python binaries that are not yet prepared for
%   64-bit computers. 

%%
% Copyright (C) 2011 Christophe Gouel
%
% Licensed under the Expat license, see <matlab:licensetohelp LICENSE.txt>
