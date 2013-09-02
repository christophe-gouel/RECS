%% Solution methods
% This page describes the default algorithm used in RECS to solve for a model's
% rational expectations equilibrium. It is possible to choose other methods by
% changing the <ug_options.html options>.
%
% The numerical algorithm used here is inspired by Fackler (2005) and Miranda
% and Fackler (1995). It is a projection method with a collocation approach
% solved by time iteration and approximation of the behavior of response
% variables.
%
%% Default algorithm
% The method will attempt to find a function that approximates well the behavior
% of response variables,
%
% $$x \approx\mathcal{X}\left(s,c_{x}\right),$$
%
% where $c_{x}$ are the parameters defining the approximation. To calculate this
% approximation, we discretize the state space, and the spline has to hold
% exactly for all points of the grid. By default, RECS uses a spline
% approximation of response variables as a function of state variables,
%
% The expectations operator is replaced by a sum over possible realization of
% shocks, $e_{l}$, to which are associated the probability $w_{l}$. If shocks
% are normal, the pairs $\left\{e_{l},w_{l}\right\}$ are calculated by RECS a
% Gaussian quadrature. Using this discretization, we can express the equilibrium
% equation as
%
% $$\underline{x}(s) \le x \le \overline{x}(s) \perp f\left(s,x,\sum_{l} w_{l} h\left(s,x,e_{l},g\left(s,x,e_{l}\right),\mathcal{X} \left(g\left(s,x,e_{l}\right),c_{x}\right)\right)\right).$$
%
% For a given approximation, $c_{x}$, and a given $s$, this equation is a
% function of $x$ only and can be solved using a mixed complementarity solver.
%
% Once all the above elements are defined, we can proceed to the algorithm,
% which runs as follows:
%
% # Initialize the spline approximation, $c_{x}^{0}$, based on a
% <first_guess.html first guess>, $x^{0}$.
% # For each point of the grid of state variables, $s_{i}$, solve for $x_{i}$
% using an <ug_solvers_eq.html MCP solver> the following equation:
% $\underline{x}(s_{i}) \le x_{i} \le \overline{x}(s_{i}) \perp f\left(s_{i},x_{i},\sum_{l}w_{l}h\left(s_{i},x_{i},e_{l},g\left(s_{i},x_{i},e_{l}\right),\mathcal{X}\left(g\left(s_{i},x_{i},e_{l}\right),c_{x}^{n}\right)\right)\right).$
% # Update the approximation using the new values of response variables,
% $x=\mathcal{X}\left(s,c_{x}^{n+1}\right)$.
% # If $|| c_{x}^{n+1}-c_{x}^{n}||_{2}\ge \epsilon$ then increment $n$
% to $n+1$ and go to step 2.
%
%% References
% <http://dx.doi.org/10.1007/s10614-005-1784-z Fackler, P. L. (2005). A MATLAB
% Solver for Nonlinear Rational Expectations Models. _Computational Economics_,
% 26(2), 173-181.>
%
% Miranda, M. J. and Fackler, P. L. (2002). _Applied Computational Economics and
% Finance_. Cambridge: MIT Press.

