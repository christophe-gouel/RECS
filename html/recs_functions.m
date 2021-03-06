%% Function Reference

%% RECS functions
% * <matlab:doc('recsAccuracy') |recsAccuracy|> - Evaluate accuracy
% of a RECS solution
% * <matlab:doc('recsAccuracySP') |recsAccuracySP|> - Evaluate accuracy
% of a RECS solution for models with subperiods
% * <matlab:doc('recsAuxiliary') |recsAuxiliary|> - Calculates auxiliary variables not included in the core model
% of a RECS solution for models with subperiods
% * <matlab:doc('recsCheck') |recsCheck|> - Check analytical derivatives against
% numerical ones
% * <matlab:doc('recsConvert') |recsConvert|> - Convert the
% interpolation structure of a model to another form
% * <matlab:doc('recsDecisionRules') |recsDecisionRules|> - Plot a model
% decision rules
% * <matlab:doc('recsdemos') |recsdemos|> - Run all RECS demonstration files
% * <matlab:doc('recsFirstGuess') |recsFirstGuess|> - Find a first guess using the
% perfect foresight solution or the first-order approximation of the model
% * <matlab:doc('recsFirstGuessSP') |recsFirstGuessSP|> - Find a first guess for
% models with subperiods
% * <matlab:doc('recsinterpinit') |recsinterpinit|> - Prepare a RECS
% interpolation structure
% * <matlab:doc('recsmodel.recsmodel') |recsmodel|> - Prepare a RECS model object
% * <matlab:doc('recsmodelsp.recsmodelsp') |recsmodelsp|> - Prepare a RECS model object for
% models with subperiods
% * <matlab:doc('recsSimul') |recsSimul|> - Simulate a model
% * <matlab:doc('recsSimulSP') |recsSimulSP|> - Simulate a model with subperiods
% * <matlab:doc('recsSolveDeterministicPb') |recsSolveDeterministicPb|> - Solve a
% perfect foresight problem
% * <matlab:doc('recsSolveDeterministicPbSP') |recsSolveDeterministicPbSP|> - Solve a
% perfect foresight problem with subperiods
% * <matlab:doc('recsSolveREE') |recsSolveREE|> - Find the rational expectations
%   equilibrium (REE) of a model
% * <matlab:doc('recsSolveREESP') |recsSolveREESP|> - Find the rational expectations
%   equilibrium (REE) of a model with subperiods
% * <matlab:doc('recsSolveREEFiniteHorizon') |recsSolveREEFiniteHorizon|> -
% Solve a finite horizon rational expectations problem
% * <matlab:doc('recsSS') |recsSS|> - Solve for the deterministic
% steady state of a model
% * <matlab:doc('recsSSSP') |recsSSSP|> - Solve for the deterministic
% steady state of a model with subperiods

%% Nonlinear equations and MCP solvers that can be used with RECS
% * <matlab:isfilepresent('fsolve') |fsolve|> - Solve system of nonlinear
%   equations (from MATLAB Optimization Toolbox)
% * <matlab:doc('lmmcp') |lmmcp|> - Solve mixed complementarity problems
% * <matlab:doc('mcpsolve') |mcpsolve|> - Solve mixed complementarity problems
% * <matlab:doc('ncpsolve') |ncpsolve|> - Solve mixed complementarity problems
%   (from CompEcon)
% * <matlab:doc('nsoli') |nsoli|> - Solve system of nonlinear
%   equations by a Jacobian-free Newton-Krylov solver
% * <matlab:doc('recspathmcp') |recspathmcp|> - Solve mixed complementarity problems
% * <matlab:doc('SA') |SA|> - Solve a system of equations by successive approximation
% * <matlab:doc('SCP') |SCP|> - Solve a problem through simple continuation
%   method (homotopy)