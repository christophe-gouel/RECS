%% Solvers for systems of nonlinear equations and for mixed complementarity problems

%% Overview of solvers for equilibrium equations
% There are two types of solver for solving equilibrium equations: solvers for
% system of nonlinear equations and solvers for system of <MCP.html mixed
% complementarity constraints (MCP)> equations. As MCP equations include as
% particular cases nonlinear equations, MCP solvers are perfectly able to solve
% also systems of traditional nonlinear equations.
%
% For each solver, options can be defined which govern their precision, solution
% method or output display. There is no common pattern for these options, as
% RECS is only providing an interface to existing solvers, so please see solvers'
% documentation for how to manage options.

%% MCP solvers
% Four MCP solvers can be used with RECS of which three are available by
% default.
%
% *|lmmcp|*
%
% |lmmcp| is RECS default solver (for both MCP and traditional
% nonlinear problems). It is included in RECS files.
%
% *|mcpsolve|*
%
% |mcpsolve| is a faster implementation of |ncpsolve| from the CompEcon
% toolbox. It is a simple MCP solver that will be very fast on simple problems
% but that may lack robustness.
%
% *|ncpsolve|*
%
% |ncpsolve| comes with the CompEcon toolbox, which has to be installed for RECS
% to run.
%
% *Path*
%
% Path is the reference solver for mixed complementarity problems. Path has to
% be installed separately. It is highly recommended to install if difficult
% complementarity problems need to be solved. To install Path, see
% <http://pages.cs.wisc.edu/~ferris/path.html>. Path is called by the MATLAB
% file |recspathmcp.m|.

%% Nonlinear equations solver
% *|fsolve|*
%
% |fsolve| is the only nonlinear equations solver interfaced with RECS. You need
% to have MATLAB Optimization Toolbox to use this solver.

%% See also
% <matlab:isfilepresent('fsolve') |fsolve|>, <matlab:doc('lmmcp') |lmmcp|>,
% <matlab:doc('mcpsolve') |mcpsolve|>, <matlab:doc('ncpsolve') |ncpsolve|>,
% <matlab:doc('recspathmcp') |recspathmcp|>.
