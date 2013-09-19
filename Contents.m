% RECS Toolbox
% Version 0.-beta (R2013b) 19-Sep-2013
%
% RECS functions
%   recsAccuracy              - Evaluate accuracy of a RECS solution.
%   recsCheck                 - Check analytical derivatives against numerical ones.
%   recsConvert               - Convert the interpolation structure of a model to another form.
%   recsDecisionRules         - Plot a model decision rules
%   recsdemos                 - Run all RECS demonstration files.
%   recsFirstGuess            - Find a first guess using the perfect foresight solution or the first-order approximation of the model
%   recsinterpinit            - Prepare a RECS interpolation structure.
%   recsmodel                 - Prepare a RECS model object.
%   recsSimul                 - Simulate a model.
%   recsSolveDeterministicPb  - Solve a perfect foresight problem.
%   recsSolveLocal            - Calculate the first-order perturbation solution
%   recsSolveREE              - Find the rational expectations equilibrium (REE) of a model.
%   recsSolveREEFiniteHorizon - Solve a finite horizon rational expectations problem.
%   recsSS                    - Solve for the deterministic steady state of a model.
%
% Nonlinear equations and MCP solvers included in RECS
%   lmmcp       - Solve mixed complementarity problems.
%   nsoli       - Solve system of nonlinear equations by a Jacobian-free Newton-Krylov solver.
%   recspathmcp - Solve mixed complementarity problems (require the installation of Path).
%   SA          - Solve a system of equations by successive approximation.
%   SCP         - Solve a problem through simple continuation method (homotopy)

% Copyright (C) 2011-2013 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt
