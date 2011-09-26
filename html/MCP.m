%% Introduction to mixed complementarity problems
% All equilibrium equations are represented in RECS as mixed complementarity
% problems (MCP), so we present here a short introduction to this kind of
% problems. For more information, see Rutherford [1], and Ferris and Pang [2].

%% Definition of a mixed complementarity problem
% Complementarity problem is an extension of square system of nonlinear
% equations. Many economic problems can be expressed as complementarity
% problems. A MCP is defined as follows (adapted from Munson [3]):
%
% *Definition (Mixed Complementarity Problem)* Given a continuously
% differentiable function $F: R^n\rightarrow R^n$, and lower and upper
% bounds
%
% $$l \in \{R \cup \{-\infty\}\}^n,$$
%
% $$u \in \{R \cup \{+\infty\}\}^n.$$
%
% The mixed complementarity problem $l \le x \le u \perp F(x)$ is to find a $x
% \in R^n$ such that one of the following holds for each $i \in \{1,\ldots,n\}:$
%
% $F_i(x)=0$ and $l_i\le x_i \le u_i$,
%
% $F_i(x)>0$ and $x_i=l_i$,
%
% $F_i(x)<0$ and $x_i=u_i$.
%
% It is easy to see that with infinite lower and upper bounds, a MCP problem is
% equivalent to solving a square system of nonlinear equations:
%
% $$-\infty\le x\le +\infty \perp F(x) \Leftrightarrow F(x)=0.$$

%% When do complementarity problems arise?
% Complementarity problems appear naturally in economics. Without pretending
% to be exhaustive, here are a few situations that lead to complementarity
% problems:
%
% * *Karush-Kuhn-Tucker conditions* of a constrained nonlinear program. The
% first-order conditions of the following nonlinear programming problem
% $\min_x f(x)$ subject to $g(x)\le 0$ and $l\le x \le u$ can be written as
% a system of complementarity equations: $l\le x\le u \perp f'(x)-\lambda
% g'(x)$ and $\lambda\ge 0 \perp -g(x)\ge 0$, where $\lambda$ is the Lagrange
% multiplier on the first inequality constraint.
% * Natural representation of *regime-switching behaviors*. For example, a
% system of intervention prices backed by public storage can be represented
% as $S\ge 0 \perp P-P^I\ge 0$, where $S$, $P$ and $P^I$ are, respectively,
% the storage, price and intervention price.
% * *A Walrasian equilibrium* can be formulated as a complementarity problem
% (Mathiesen [4]).

%% Conventions of notations adopted for representing and solving complementarity problems
% For solving complementarity problems, RECS uses several solvers listed in
% <ug_solvers_eq.html MCP solvers>. The convention adopted in most MCP
% solvers and used by RECS is the one used above in MCP definition: superior
% or equal inequalities are associated with the lower bounds and inferior or
% equal inequalities are associated with the upper bounds. So, when defining
% your model, be careful to respect this convention.

%% References
%
% [1] <http://dx.doi.org/10.1016/0165-1889(94)00831-2 Rutherford,
% T. F. (1995). Extension of GAMS for complementarity problems arising in
% applied economic analysis. _Journal of Economic Dynamics and Control_, 19(8),
% 1299-1324.>
%
% [2] <http://www.jstor.org/pss/2132696 Ferris, M. C. and Pang,
% J. S. (1997). Engineering and economic applications of complementarity
% problems. _Siam Review_, 39(4), 669-713.>
%
% [3] Munson, T. (2002). Algorithms and Environments for
% Complementarity. Unpublished PhD thesis from University of Wisconsin, Madison.
%
% [4] <http://dx.doi.org/10.1007/BF02591680 Mathiesen, L. (1987). An algorithm
% based on a sequence of linear complementarity problems applied to a Walrasian
% equilibrium model: An example. _Mathematical Programming_, 37, 1-18.>

%%
% Copyright (C) 2011 Christophe Gouel
%
% Licensed under the Expat license, see <matlab:filetohelp('LICENSE.txt') LICENSE.txt>
