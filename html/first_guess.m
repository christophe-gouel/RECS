%% First guess
% The solution methods used in RECS to find the rational expectations
% equilibrium of a problem rely all on nonlinear equation solvers, with which
% convergence can be achieved only if the starting point is not too far from the
% solution. Hence, the importance of providing a good first guess.
%
%% The perfect foresight problem as first guess
% RECS can attempt to calculate for you a first guess that is good in most
% situations. It does it by calculating the perfect foresight solution of the
% deterministic problem in which the shocks in the stochastic problem have been
% substituted by their expectations:
%
% $$\underline{x}(s) \le x \le \overline{x}(s) \perp f(s,x,z),$$
%
% $$z =  h(s,x,\mathrm{E}\left(e\right),s_{+},x_{+}),$$
%
% $$s = g(s_{-},x_{-},\mathrm{E}\left(e\right)).$$
%
% The solver for perfect foresight problems (see <deterministic.html
% Deterministic problems>) assumes that the problem converges to its
% deterministic steady state after T periods (by default T=50). The perfect
% foresight problem is solved for each grid point and the resultant behavior is
% used as a first guess for the stochastic problem. This is done by the
% following call:
% 
%  [interp,x] = recsFirstGuess(interp,model,s,sss,xss);
%
% This function updates the interpolation structure, |interp|, with the solution
% of the perfect foresight problem and outputs the value of response variables
% on the grid, |x|.
%
% *Solving time*
% 
% As the perfect foresight problem is solved for each point of the grid and for
% T periods, this step can take some time. In many cases, it should be expected
% that it may take more time to find a first guess through the perfect foresight
% solution than to solve the stochastic problem from this first guess.
%
% *How good is this first guess?*
%
% It depends on the model and the solution horizon. For most models, if the
% state space has been properly defined (i.e., not too small or too large), this
% first guess is good enough to allow the stochastic solver to converge.
%
% The quality of the perfect foresight solution as a first guess depends also on
% the model nonlinearity. For models with behavior close to linear, the first
% guess can be extremely good. For example, with <gro1.html the stochastic
% growth model> the first guess leads to an initial deviation from rational
% expectations of 1E-5. So after a few iterations, the solver converges to the
% solution. However, with <gro2.html the stochastic growth model with
% irreversible investment>, which is much more nonlinear, the residual when
% starting the stochastic solver from the first guess is 1E-2. This is
% sufficient to achieve convergence, but it is still far from the true solution.
% 
%% User-provided first guess
% It is not necessary to use the deterministic problem as a first guess. A first
% guess can be provided directly by the user. In this case there is no need to
% call a function, one just need to define a n-by-m matrix |x| that provides the
% values of the m response variables for the n points on the grid of state
% variable. The matrix |x| should then be supplied to |recsSolveREE| as would be
% the one created by |recsFirstGuess| (see the next step <solve_REE.html Solve
% the rational expectations problem> for more information).
%
% For example of this, see in demos <cs1.html CS1>, <sto1.html STO1>, <sto4.html
% STO4>, <sto5.html STO5>, and <sto6.html STO6>.




