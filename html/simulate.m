%% Simulate

% to speed-up: possible to start from several initial points at once

% Random number from funrand (or without it from random draw from discretize shocks based on their probability).
% reset random number generator

% Long-run statistics (exclude first 20 obs): mean of autocorrelation

% Variables organization

%% Choice of simulation techniques
% They are two principal approaches to simulate the model once an approximated
% rational expectations equilibrium has been found.
%
% *Use approximated decision rule*
%
% The first, and most commonly used, method consists in simulating the model by
% applying the approximated decision rules recursively. 
%
% Starting from a known $s_1$, for $t=1:T$
%
% $$x_t=\mathcal{X}\left(s_t,c_{x}\right)$$
%
% Draw a shock realization $e_{t+1} from its distribution and update the state
% variable:
%
% $$s_{t+1}=g\left(s_t,x_t,e_{t+1}\right)$$
%
% *Solve equilibrium equations using approximated expectations*
%
% The second uses the function being approximated (the decisions rules, the
% conditional expectations, or the expectations function as shown in
% <ug_methods.html Solution methods>) to approximate next-period expectations
% and solve the equilibrium equations to find the current decisions (in a
% time-iteration approach). For example, if the approximated expectations
% ($z\approx\mathcal{Z}\left(s,c_z\right)$) are used to replace the conditional
% expectations in the equilibrium equation, the algorithm runs as follows:
%
% Starting from a known $s_1$, for $t=1:T$, find $x_t$ by solving the following
% MCP equation using as first guess $x_t=\mathcal{X}\left(s_t,c_{x}\right)$
%
% $$\underline{x}(s_t) \le x_t \le \overline{x}(s_t) \perp f\left(s_t,x_t,\mathcal{Z}\left(s_t,c_z\right)\right)$$
%
% Draw a shock realization $e_{t+1} from its distribution and update the state
% variable:
%
% $$s_{t+1}=g\left(s_t,x_t,e_{t+1}\right)$$
%
% The two approaches differ in speed and in precision. Simulating a model
% through the recursive application of approximated decision rules is very fast
% because it does not require any nonlinear solve. The second approach is much
% slower because each period requires solving a system of complementarity
% equations. However, its precision is much higher, because the approximation is
% used only for the expectations of next-period conditions. The precision gain
% is all the more important when decision rules have kinks, which makes them
% difficult to approximate. Wright and Williams (1984), who first proposed the
% parameterized expectations algorithm, suggest using the second method for
% simulating the storage model.
%
% *When does this distinction matter?*
%
% Most of the time, you should use the default approach of applying recursively
% the approximated decision rules. Simulating the model by solving the equations
% should be considered in the following two cases:
%
% * If the model has been solved with a low precision, using the approximated
% decision rules can lead to large errors, while simulating the model by solving
% the equilibrium equations using only the approximation in expectations will
% yield a much more precise solution. It is the case in the following example
% inspired from <sto2.html STO2>: a price-floor policy backed by public
% storage. If the problem is approximated by spline with a small grid of 6
% points, it is not possible to expect a good approximation of the highly
% nonlinear storage rules. However, if the approximation is used only in
% expectations, it yields a very precise solution:
sto2simu(1)
%%
% * If you are interested in the percentage of time spent in various regimes,
% solving the equilibrium equations is the only way to get a precise
% estimate. Even when solved with a high precision, a simulation using the
% approximated decision rules provides limited precision regarding the
% percentage of time spent in each regime. Indeed, approximations by spline or
% Chebyshev polynomials fluctuate around the exact value between collocation
% nodes. So between two nodes at which the model should be at a bound, the
% approximated decision rule can yield a result very close to the bound without
% actually satisfy it. Below are the moments obtained from simulating the model
% in <sto2.html STO2> with the both methods and with decision rules approximated
% on 30 nodes. The moments are quite similar, but the approximated decision
% rules underestimate widely the time spent at the bounds:
sto2simu(2)

%% References
% <http://www.jstor.org/stable/1885726 Wright, B. D. and Williams, J. C. (1984).
% The Welfare Effects of the Introduction of Storage. _The Quarterly Journal of
% Economics_, 99(1), 169-192.>