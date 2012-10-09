%% Introduction to mixed complementarity problems
% To be able to reliably solve models with occasionally binding constraints, all
% equilibrium equations are represented in RECS as mixed complementarity
% problems (MCP), so we present here a short introduction to this kind of
% problems. For more information, see Rutherford (1995), and Ferris and Pang
% (1997).

%% Definition of a mixed complementarity problem
% Complementarity problem can be seen as an extension of square system of
% nonlinear equations that incorporates a mixture of equations and
% inequalities. Many economic problems can be expressed as complementarity
% problems. An MCP is defined as follows (adapted from Munson, 2002):
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
% To sum up, an MCP problem associates each variable, $x_i$, to a lower bound,
% $l_i$, an upper bound, $u_i$, and an equation, $F_i(x)$. The solution is such
% that if $x_i$ is between its bounds then $F_i(x)=0$. If $x_i$ is equal to its
% lower (upper) bound then $F_i(x)$ is positive (negative).
%
% This format encompasses several cases. In particular, it is easy to see that
% with infinite lower and upper bounds, solving an MCP problem is equivalent to
% solving a square system of nonlinear equations:
%
% $$-\infty\le x\le +\infty \perp F(x) \Leftrightarrow F(x)=0.$$

%% A simple example of mixed complementarity
% We consider here the traditional consumption/saving problem with borrowing
% constraint (Deaton, 1991). This problem is solved as demonstration, see
% <cs1.html CS1>. A consumer with utility $\sum_{t=0}^{\infty}
% u(C_t)/(1+\delta)^t$ has a stochastic income, $Y_t$ and has to choose each
% period how much to consume and how much to save. He is prevented from
% borrowing, but can save at the rate $r$. Without the borrowing constraint, his
% problem consists in choosing $C_t$ such that it satisfies the standard Euler equation:
%
% $$u'\left(C_{t}\right)=\frac{1+r}{1+\delta}\mathrm{E}_{t}\left[u'\left(C_{t+1}\right)\right].$$
%
% The borrowing constraint implies the inequality $C_t \le X_t$, where $X_t$ is
% the cash on hand available in period $t$ and defined as
%
% $$X_{t}=\left(1+r\right)\left(X_{t-1}-C_{t-1}\right)+Y_{t}.$$
%
% When the constraint is binding (i.e., $C_t=X_t$), the Euler equation does no
% longer hold. The consumer would like to borrow but cannot, so its marginal
% utility of immediate consumption is higher than its discounted marginal utility
% over future consumption:
%
% $$u'\left(C_{t}\right) \ge \frac{1+r}{1+\delta}\mathrm{E}_{t}\left[u'\left(C_{t+1}\right)\right].$$
%
% Since the Euler equation holds when consumption is not constrained, this
% behavior can be summarized be the following complementarity equation
%
% $$C_{t}\le X_{t} \quad \perp \quad \frac{1+r}{1+\delta}\mathrm{E}_{t}\left[u'\left(C_{t+1}\right)\right]\le u'\left(C_{t}\right).$$

%% When do complementarity problems arise?
% In addition to encompass nonlinear systems of equations, complementarity
% problems appear naturally in economics. Without pretending to be exhaustive,
% here are a few situations that lead to complementarity problems:
%
% * *Karush-Kuhn-Tucker conditions* of a constrained nonlinear program. The
% first-order conditions of the following nonlinear programming problem
% $\min_x f(x)$ subject to $g(x)\le 0$ and $l\le x \le u$ can be written as
% a system of complementarity equations: $l\le x\le u \perp f'(x)-\lambda
% g'(x)$ and $\lambda\ge 0 \perp -g(x)\ge 0$, where $\lambda$ is the Lagrange
% multiplier on the first inequality constraint.
% * Natural representation of *regime-switching behaviors*. Let's consider two
% examples. First, a variable $y$ that is determined as a function of another
% variable $x$ as long as $y$ is superior to a lower bound $a$. Put simply:
% $y=\max\left[a,\lambda\left(x\right)\right]$. In MCP, we would write this as
% $y\ge a \perp y\ge \lambda\left(x\right)$. Second, a system of intervention
% price in which a public stock is accumulated when prices decrease below an
% intervention price and sold when they exceed the intervention price (see also
% <sto2.html STO2>) can be represented as $S\ge 0 \perp P-P^I\ge 0$, where $S$,
% $P$ and $P^I$ are, respectively, the public storage, the price and the
% intervention price.
% * *A Walrasian equilibrium* can be formulated as a complementarity problem
% (Mathiesen, 1987) by associating three sets of variables to three sets of
% inequalities: non-negative levels of activity are associated to zero-profit
% conditions, non-negative prices are associated to market clearance, and income
% levels are associated to income balance.

%% Conventions of notations adopted for representing and solving complementarity problems
% For solving complementarity problems, RECS uses several solvers listed in
% <ug_solvers_eq.html MCP solvers>. The convention adopted in most MCP
% solvers and used by RECS is the one used above in MCP definition: superior
% or equal inequalities are associated with the lower bounds and inferior or
% equal inequalities are associated with the upper bounds. So, when defining
% your model, be careful to respect this convention.
%
% In addition, in Yaml files inequalities should not be written
% explicitly. As an example, let's consider these 5 equations:
%
% $$ F(x) = 0,$$
%
% $$ y\ge a \quad \perp \quad G\left(y\right)\ge 0,$$
%
% $$ z\le b \quad \perp \quad H\left(z\right)\le 0,$$
%
% $$ a\le k \le b \quad \perp \quad J\left(k\right),$$
%
% $$ l\ge a \quad \perp \quad M\left(l\right)\le 0.$$
%
% They would be written in a Yaml file as
%
%    - F(x)=0 | -inf<=x<=inf
%    - G(y)   |    a<=y<=inf
%    - H(z)   | -inf<=z<=b
%    - J(k)   |    a<=k<=b
%    - -M(l)  |    a<=l<=inf
%
% As you can note, one needs to associate lower and upper bounds with every
% variables. So if there are no finite bounds, one can associate infinite
% bounds. The last equation was not respecting the convention that associates
% lower bounds on variables with superior or equal inequality, so one should be
% careful when writing it in the Yaml file to reverse the sign of the equation.

%% References
% <http://www.jstor.org/stable/2938366 Deaton, A. (1991). Saving and liquidity
% constraints. _Econometrica_, 59(5), 1221-1248.>
%
% <http://dx.doi.org/10.1016/0165-1889(94)00831-2 Rutherford,
% T. F. (1995). Extension of GAMS for complementarity problems arising in
% applied economic analysis. _Journal of Economic Dynamics and Control_, 19(8),
% 1299-1324.>
%
% <http://www.jstor.org/pss/2132696 Ferris, M. C. and Pang,
% J. S. (1997). Engineering and economic applications of complementarity
% problems. _Siam Review_, 39(4), 669-713.>
%
% Munson, T. (2002). Algorithms and Environments for
% Complementarity. Unpublished PhD thesis from University of Wisconsin, Madison.
%
% <http://dx.doi.org/10.1007/BF02591680 Mathiesen, L. (1987). An algorithm based
% on a sequence of linear complementarity problems applied to a Walrasian
% equilibrium model: An example. _Mathematical Programming_, 37, 1-18.>
%
% <html>
% <hr/>
% </html>
%
%%
% Copyright (C) 2011-2012 Christophe Gouel
%
% Licensed under the Expat license, see <LICENSE.txt>
