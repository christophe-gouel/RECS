%% Solvers for equilibrium equations

%% Overview of solvers for equilibrium equations
% There are two types of solver for solving equilibrium equations: solvers for
% system of nonlinear equations and solvers for system of mixed complementarity
% constraints (MCP) equations. As MCP equations include as particular cases
% nonlinear actions, the latter type is perfectly able to solve systems of
% traditional nonlinear equations.
%
% For each solver, options can be defined which govern their precision, solution
% method or output display. There is no common pattern for these options, so
% please see solvers' documentation for how to manage options.

%% MCP solvers
% Three MCP solvers can be used with RECS of which two are available by
% default. |lmmcp| is RECS default solver (for both MCP and traditional
% problems). It is included in RECS files. |ncpsolve| comes with the CompEcon
% toolbox, which has to be installed for RECS to run. Path has to be installed
% separately. To install Path, see
% <http://pages.cs.wisc.edu/~ferris/path.html>. To distinguish Path from the
% MATLAB internal command, the name of the Path MATLAB file is |recspathmcp.m|.

%% Nonlinear equations solver
% |fsolve| is the only nonlinear equations solver interfaced with RECS. You need
% to have MATLAB Optimization Toolbox to use this solver.

%% See also
% <matlab:doc('fsolve') |fsolve|> <matlab:doc('lmmcp') |lmmcp|>
% <matlab:doc('ncpsolve') |ncpsolve|> <matlab:doc('recspathmcp') |recspathmcp|>
%
% <html>
% <hr/>
% </html>
%
%%
% Copyright (C) 2011-2012 Christophe Gouel
%
% Licensed under the Expat license, see <LICENSE.txt>
