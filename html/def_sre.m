%% Definition of a stochastic rational expectations problem

%% Rational expectation models
% We adopt here the framework proposed in Fackler [1], and used also in Winschel et
% Krätzig [2]. A model can be defined by the following three equations, where
% next-period variables are indicated with a dot on top of the character
%
% $\underline{x}(s) \le x \le \overline{x}(s) \perp f(s,x,z)$, where
% $f:R^{d+m+p}\rightarrow R^{m}$,
%
% $z = \mathrm{E}_{\dot e} \left[h(s,x,\dot e,\dot s,\dot x)\right]$, where
% $h:R^{d+m+q+d+m}\rightarrow R^{p}$,
%
% $\dot s = g(s,x,\dot e)$, where $g:R^{d+m+q}\rightarrow R^{d}$.
%
% Variables have been partitioned into state variables, $s$, response variables,
% $x$, and shocks, $e$. Response variables can have lower and upper bounds,
% $\underline{x}$ and $\overline{x}$, which can themselves be function of state
% variables. Expectations variables, denoted by $z$, are defined in addition since
% they are necessary for solving the model.
%
% The first equation is the equilibrium equation. It characterizes the behavior of
% response variables given state variables and expectations about next period. For
% generality, it is expressed as a mixed complementarity problem (MCP, if you
% ignore what MCP means, see <MCP.html Introduction to mixed complementarity
% problems>). In cases where response variables have no lower and upper bounds,
% it simplifies to a traditional equation:
%
% $0=f(s,x,z)$.
%
% The second equation defines the expectations. The last equation is the
% transition equation.

%% Solving a rational expectations model
% The problem created by rational expectations models is that the second
% equation, defining the expectations, is not a traditional algebraic
% equation. It is an equation that expresses the consistency between agents'
% expectations, their information set and realized outcomes. Solving a rational
% expectations model amounts to find an approximation or an algebraic
% representation of the expectations terms so that the model can be reduced to
% something easier to solve.
%
% For example, if it is possible to find an approximation of the relationship
% between expectations and current-period state variables (i.e., the parameterized
% expectations approach [3]), the equilibrium equation can be simplified to
%
% $$\underline{x}(s) \le x \le \overline{x}(s) \perp f\left(s,x,\phi(s)\right),$$
%
% where $z \approx \phi(s)$. This equation can be solved for x with any
% <ug_solvers_eq.html MCP solvers>.
%
% The RECS solver implements various methods to solve rational expectations
% models and find approximation of expectations terms. See <ug_methods.html
% Solution methods> for further information.

%% References
%
% [1] <http://dx.doi.org/10.1007/s10614-005-1784-z Fackler, P. L. (2005). A
% MATLAB Solver for Nonlinear Rational Expectations Models. _Computational
% Economics_, 26(2), 173-181.>
%
% [2] <http://dx.doi.org/10.3982/ECTA6297 Winschel, V. and Krätzig,
% M. (2010). Solving, Estimating, and Selecting Nonlinear Dynamic Models Without
% the Curse of Dimensionality. _Econometrica_, 78(2), 803-821.>
%
% [3] <http://www.jstor.org/stable/1391746 den Haan, W. J. and Marcet, A.
% (1990). Solving the Stochastic Growth Model by Parameterizing
% Expectations. _Journal of Business & Economic Statistics_, 8(1), 31-34.>

%%
% Copyright (C) 2011-2012 Christophe Gouel
%
% Licensed under the Expat license, see <LICENSE.txt>
