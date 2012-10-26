%% Simulate the model

%% Running a simulation
% Once the rational expectations equilibrium has been found, it is possible to
% simulate the model using the function |recsSimul| with the following call:
%
%  [ssim,xsim] = recsSimul(model,interp,s0,nper);
%
% where |model| and |interp| have been defined previously as the model and the
% interpolation structure. |s0| is the initial state and |nper| is the number of
% simulation periods.
%
% |s0| is not necessarily a unique vector of state variables. It is possible to
% provide a matrix of initial states from which simulations will be carried for
% |nper| periods. Even when starting from a unique state, it is useful to use
% this feature to speed-up simulations.
%
% *Shocks* 
%
% The shocks used in the simulation can be provided as the fifth input to the function:
%
%  [ssim,xsim] = recsSimul(model,interp,s0,nper,shocks);
%
% However, by default, |recsSimul| uses the function |funrand| defined in the
% model structure (see <ug_model_struct.html Define the model structure>) to
% draw random shocks. If this function is not provided random draws are made
% from the shock discretization, using the associated probabilities.
%
% If you want to reproduce previously run results, remember to reset the random
% number generator using the Matlab function |reset|.
%
% *Asymptotic statistics*
%
% If in the options the field |stat| is set to 1 or if 5 arguments are required
% as output of |recsSimul|, then some statistics over the asymptotic
% distribution are calculated. For this, the first 20 observations are
% discarded. The statistics calculated are the mean, standard deviation,
% skweness, kurtosis, minimum, maximum, percentage of time spent at the lower
% and upper bounds, correlation matrix, and the 5 first-order autocorrelation
% coefficients. In addition, histograms of the variables distribution are drawn.
%
% These statistics are available as a structure in the fifth output of
% |recsSimul|:
%
%  [ssim,xsim,esim,fsim,stat] = recsSimul(model,interp,s0,nper);
% 
% Since RECS does not retain the variables names, when displaying the statistics
% variables are organized as follows: first state variables, then response
% variables, and for both they follow the order of their definition in the Yaml
% file.

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
% This is the default simulation method. It is chosen in the options structure
% by setting the field |simulmethod| to |interpolation|.
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
% Draw a shock realization $e_{t+1}$ from its distribution and update the state
% variable:
%
% $$s_{t+1}=g\left(s_t,x_t,e_{t+1}\right)$$
%
% This simulation method is chosen in the options structure by setting the field
% |simulmethod| to |solve|.
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