% RECS Toolbox
% Version 0.7 (R2017b) 5-Jan-2018
%
% RECS functions
%   recsAccuracy              - Evaluate accuracy of a RECS solution.
%   recsAccuracySP            - Evaluate accuracy of a RECS solution for models with subperiods.
%   recsAuxiliary             - Calculates auxiliary variables not included in the core model.
%   recsCheck                 - Check analytical derivatives against numerical ones.
%   recsConvert               - Convert the interpolation structure of a model to another form.
%   recsDecisionRules         - Plot a model decision rules.
%   recsdemos                 - Run all RECS demonstration files.
%   recsFirstGuess            - Find a first guess using the perfect foresight solution or the first-order approximation of the model.
%   recsFirstGuessSP          - Find a first guess for models with subperiods.
%   recsinterpinit            - Prepare a RECS interpolation structure.
%   recsmodel                 - Prepare a RECS model object.
%   recsmodelsp               - Prepare a RECS model object for models with subperiods.
%   recsSimul                 - Simulate a model.
%   recsSimulSP               - Simulate a model with subperiods.
%   recsSolveDeterministicPb  - Solve a perfect foresight problem.
%   recsFirstGuessSP          - Find a first guess for models with subperiods.
%   recsSolveLocal            - Calculate the first-order perturbation solution.
%   recsSolveREE              - Find the rational expectations equilibrium (REE) of a model.
%   recsSolveREESP            - Find the rational expectations equilibrium (REE) of a model with subperiods.
%   recsSolveREEFiniteHorizon - Solve a finite horizon rational expectations problem.
%   recsSS                    - Solve for the deterministic steady state of a model.
%   recsSSSP                  - Solve for the deterministic steady state of a model with subperiods.
%
% Nonlinear equations and MCP solvers included in RECS
%   lmmcp       - Solve mixed complementarity problems.
%   mcpsolve    - Solve mixed complementarity problems.
%   nsoli       - Solve system of nonlinear equations by a Jacobian-free Newton-Krylov solver.
%   recspathmcp - Solve mixed complementarity problems (require the installation of Path).
%   SA          - Solve a system of equations by successive approximation.
%   SCP         - Solve a problem through simple continuation method (homotopy).

% Copyright (C) 2011-2018 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt
