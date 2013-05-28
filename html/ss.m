%% Steady state
% In contrast to pertubation approaches, the deterministic steady state does not
% play an important role in finding the rational expectations equilibrium with
% collocation methods, and thus with RECS. However, finding the deterministic
% steady state proved very useful in the overall model building for
%
% * Calibration
% * Definition of the approximation space
% * Calculation of a first guess using the corresponding perfect foresight problem
% * Checking model structure

%% Steady state definition
% The deterministic steady state is the state reached in the absence of shocks
% and ignoring future shocks. Following the convention adopted in RECS (see
% <def_sre.html Definition of a stochastic rational expectations problem>), the
% deterministic steady state is the set $\left\{s,x,z\right\}$ of state,
% response and expectations variables that solves the following system of
% equations
%
% $\underline{x}(s) \le x \le \overline{x}(s) \perp f(s,x,z)$
%
% $z = h(s,x,\mathrm{E}\left(e\right),s,x)$
%
% $s = g(s,x,\mathrm{E}\left(e\right))$

%% Finding the steady state with RECS
% *Automatically when initializing model structure*
%
% When writing a model file (see <ug_model_files.html Writing RECS model
% files>), it is possible, at the end of the Yaml file in the |calibration|
% block, to define an initial guess for finding the steady state. When the model
% structure is created by |recsmodel|, if the definition of the shocks is
% provided to |recsmodel|, a Newton-type solver will attempt to find the steady
% state starting from the initial guess provided in the model file. If a steady
% state is found, it is displayed in MATLAB command window.
%
% *Manually*
%
% Otherwise, the steady state can be found manually by feeding the function
% <matlab:doc('recsSS') |recsSS|> with the model and an initial guess for the
% steady state.
%
% Both approaches rely on a Newton-type solver to find the steady state. See
% <ug_solvers_eq.html Solvers for systems of nonlinear equations and for mixed
% complementarity problems> for details on solver choice.

%% Uses of the deterministic steady state with RECS
%
% *Calibration*
%
% Compared to the stochastic rational expectations problem, the deterministic
% steady state is easy to find. It does not even require the definition of an
% interpolation structure. Since in practice it is often close to the long-run
% average values of the stochastic problem, the steady state is useful for
% calibrating the model so that, on the asymptotic distribution, the model can
% reproduce the desired long-run average.
%
% *Define the approximation space*
%
% As the values of the state variables in the stochastic problem are often
% located around the deterministic steady state, the steady state serves as a
% good reference point around which the state space can be define. See
% <ug_interpolation.html Defining the interpolation structure> for more.
%
% *First guess calculation for the stochastic problem*
%
% The rational expectations solver requires a good first guess to ensure
% convergence to the solution. RECS proposes to provide a first guess by
% calculating the perfect foresight problem corresponding to the stochastic
% problem. The perfect foresight problem assumes that the model converges in the
% long-run to its deterministic steady state. See <first_guess.html> for more.
%
% *Checking model structure*
%
% Before solving the stochastic problem, even if there are no motives for using
% the steady state, it is good practice to find it in order to ensure that the
% model written in the Yaml file behaves as expected in this simple case.
