%% Define the interpolation structure
% Following the previous two steps (<ug_model_files.html Writing RECS model
% files> and <ug_model_struct.html Defining the model object>) to completely
% define the problem, it remains only to define the domain over which it will be
% approximated and the precision of the approximation.

%% Create the interpolation structure
% This task is done by the function |recsinterpinit|, which requires at least
% three inputs: the number of points on the grid of approximation, the lower
% bounds and the upper bounds of the grid.
%
% The structure of the call to |recsinterpinit| is then
%
%  [interp,s] = recsinterpinit(n,smin,smax,method);
%
% The inputs are as follows: |n| designates the order of approximation (if it is
% a scalar, the same order is applied for all dimensions), |smin| and |smax| are
% size-d vectors of left and right endpoints of the state space, and the
% (optional) string |method| defines the interpolation method, either spline
% (|'spli'|, default), or Chebyshev polynomials (|'cheb'|).
%
% This function call returns two variables: the structure |interp|, which
% defines the interpolation structure, and the matrix |s|, which represents the
% state variables on the grid.

%% An example
% We now define the interpolation structure for the stochastic growth model
% example (<gro1.html GRO1>). Using the model object defined in
% <ug_model_struct.html the preceding step>, we choose bounds for capital 15%
% below and above the steady-state value (|model.sss(1)|), and for productivity,
% which follows an AR(1), we choose 4 times below minimum and 4 times above
% maximum discretized shocks (|model.shocks.e|):
smin          = [0.85*model.sss(1) min(model.shocks.e)*4];
smax          = [1.15*model.sss(1) max(model.shocks.e)*4];
%%
% Using 10 nodes for each dimension and Chebyshev polynomials, the function call
% is:
[interp,s] = recsinterpinit(10,smin,smax,'cheb');
%%
% The interpolation structure has the following fields:
disp(interp)

%% Choice Spline/Chebyshev polynomials
% To be completed

%% Notes on functions used
% For interpolation, RECS relies entirely on the tools programmed in the
% CompEcon toolbox. So, to know more about the underlying functions, see
% CompEcon documentation and Miranda and Fackler (2002).

%% References
% Miranda, M. J. and Fackler, P. L. (2002). _Applied Computational Economics and
% Finance_. Cambridge: MIT Press.

