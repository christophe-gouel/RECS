%% Simulate the model

%% Running a simulation
% Once the rational expectations equilibrium has been found, it is possible to
% simulate the model using the function <matlab:doc('recsSimul') |recsSimul|>
% with the following call:
%
%  [ssim,xsim] = recsSimul(model,interp,s0,nper);
%
% where |model| and |interp| have been defined previously as the model and the
% interpolation structure. |s0| is the initial state and |nper| is the number of
% simulation periods.
%
% |s0| is not necessarily a unique vector of the state variables. It is possible
% to provide a matrix of initial states on the basis of which simulations can be
% run for |nper| periods. Even starting from a unique state, using this feature
% will speed up the simulations.
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
% To reproduce previously run results, it is necessarily to reset the random
% number generator using the MATLAB function |reset|.
%
% *Asymptotic statistics*
%
% If in the options the field |stat| is set to 1, or if five arguments are
% required as the output of |recsSimul|, then some statistics over the
% asymptotic distribution are calculated. The first 20 observations are
% discarded. The statistics calculated are the mean, standard deviation,
% skewness, kurtosis, minimum, maximum, percentage of time spent at the lower
% and upper bounds, correlation matrix, and the five first-order autocorrelation
% coefficients. In addition, |recsSimul| draws the histograms of the variables
% distribution.
%
% These statistics are available as a structure in the fifth output of
% |recsSimul|:
%
%  [ssim,xsim,esim,fsim,stat] = recsSimul(model,interp,s0,nper);
%
% Since RECS does not retain the variable names, when displaying the statistics,
% the variables are organized as follows: first state variables, followed by
% response variables, both of which follow the order of their definition in the
% Yaml file.

%% Choice of simulation techniques
% There are two main approaches to simulate the model once an approximated
% rational expectations equilibrium has been found.
%
% *Using approximated decision rule*
%
% The first, and most common, method consists of simulating the model by
% applying the approximated decision rules recursively.
%
% Starting from a known $s_1$, for $t=1:T$
%
% $$x_t=\mathcal{X}\left(s_t,c_{x}\right)$$
%
% Draw a shock realization $e_{t+1}$ from its distribution and update the state
% variable:
%
% $$s_{t+1}=g\left(s_t,x_t,e_{t+1}\right)$$
%
% This is the default simulation method. It is chosen in the options structure
% by setting the field |simulmethod| to |interpolation|.
%
% *Solve equilibrium equations using approximated expectations*
%
% The second method uses the function being approximated (decisions rules,
% conditional expectations, or expectations function as shown in
% <ug_methods.html Solution methods>) to approximate next-period expectations
% and solve the equilibrium equations to find the current decisions (in a
% time-iteration approach). For example, if approximated expectations
% ($z\approx\mathcal{Z}\left(s,c_z\right)$) are used to replace conditional
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
% The two approaches differ in their speed and in precision. Simulating a model
% through the recursive application of approximated decision rules is very fast
% because it does not require solving nonlinear equations. The second approach
% is slower because each period requires that a system of complementarity
% equations be solved. However, its level of precision is much higher, because
% the approximation is used only for the expectations of next-period
% conditions. The precision gain is even more important when decision rules have
% kinks, which makes them difficult to approximate. Wright and Williams (1984),
% who first proposed the parameterized expectations algorithm, suggest using the
% second method to simulate the storage model.
%
% *When does this distinction matter?*
%
% Most of the time, the default approach of recursively applying the
% approximated decision rules should be used. Simulating the model by solving
% the equations should be considered in the following two cases:
%
% * If the model has been solved with low precision, using the approximated
% decision rules can lead to large errors, while simulating the model by solving
% the equilibrium equations using only the approximation in expectations will
% yield a more precise solution. This applies to the following example inspired
% from <sto2.html STO2>: a price-floor policy backed by public storage. If the
% problem is approximated by a spline with a small grid of 6 points, it is not
% possible to expect a good approximation of the highly nonlinear storage
% rules. However, if the approximation is used only in expectations, it yields a
% very precise solution:
sto2simu(1)
%%
% * If you are interested in the percentage of time spent in various regimes,
% solving the equilibrium equations is the only way to achieve a precise
% estimate. Even when solved with high precision, a simulation using the
% approximated decision rules provides limited precision regarding the
% percentage of time spent in each regime. Indeed, approximations using spline
% or Chebyshev polynomials fluctuate around the exact value between collocation
% nodes. So between two nodes at which the model should be at a bound, the
% approximated decision rule can yield a result very close to the bound without
% actually satisfying it. Below are the moments obtained from simulating the
% model in <sto2.html STO2> with both methods and with decision rules
% approximated on 30 nodes. The moments are quite similar, but the approximated
% decision rules widely underestimate the time spent at the bounds:
sto2simu(2)

%% References
% <http://www.jstor.org/stable/1885726 Wright, B. D. and Williams, J. C. (1984).
% The Welfare Effects of the Introduction of Storage. _The Quarterly Journal of
% Economics_, 99(1), 169-192.>