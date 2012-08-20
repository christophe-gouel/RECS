function [interp,X,exitflag] = recsSolveREEFiniteHorizon(interp,model,s,x,xT,T,options)
% RECSSOLVEREE finds the rational expectations equilibrium (REE) of a model
%
% RECSSOLVEREE implementes various approximation schemes, and equation solvers to
% find the REE of a model.
%
% INTERP = RECSSOLVEREE(INTERP,MODEL,S,X) tries to find the rational expectations
% equilibrium of the model defined in the structure MODEL, by using the
% interpolation structure defined in the structure INTERP. The problem is solved
% on the grid of state variables provided in matrix S. Matrix X is used as a
% first guess of response variables on the grid. RECSSOLVEREE returns the
% interpolation structure containing the coefficient matrices cx
% and cz, and ch if this field was initially included in INTERP.
% INTERP is a structure, which has to include the following field:
%    fspace       : a definition structure for the interpolation family (created
%                   by the function fundef)
% Optionally INTERP can also include first guess for the coefficients of
% approximation. If absent, an approximation is made from X.
%    ch, cx or cz : a coefficient matrix providing a first guess of the
%                   approximation of the expectations function for ch, of the
%                   response variables for cx, or of the expectations for cz
% MODEL is a structure, which has to include the following fields:
%    [e,w]  : discrete distribution with finite support with e the values and w the
%             probabilities (it could be also the discretisation of a continuous
%             distribution through quadrature or Monte Carlo drawings)
%    func   : function name or anonymous function that defines the model's equations
%    params : model's parameters, it is preferable to pass them as a cell array
%             (compulsory with the functional option) but other formats are
%             acceptable
%
% INTERP = RECSSOLVEREE(INTERP,MODEL,S,X,OPTIONS) solves the problem with the
% parameters defined by the structure OPTIONS. The fields of the structure are
%    display          : 1 to show iterations (default: 1)
%    eqsolver         : 'fsolve', 'lmmcp' (default), 'ncpsolve' or 'path'
%    eqsolveroptions  : options structure to be passed to eqsolver
%    extrapolate      : 1 or 2 if extrapolation is allowed outside the
%                       interpolation space, 0 or -1 to forbid it (default: 1).
%                       For -1 and 2, RECSSOLVEREE displays a warning if state
%                       variables exceed the interpolation space.
%    funapprox        : 'expapprox', 'expfunapprox', 'resapprox-simple'
%                       or 'resapprox-complete' (default)
%    functional       : 1 if the equilibrium equations are a functional equation
%                       problem (default: 0)
%    loop_over_s      : 0 (default) to solve all grid points at once or 1 to loop
%                       over each grid points
%    reemethod        : 'iter' (default), 'iter-newton' or '1-step'
%    reesolver        : 'krylov', 'mixed', 'SA' (default) or 'fsolve' (in test)
%    reesolveroptions : options structure to be passed to reesolver
%    useapprox        : (default: 1) behaviour dependent of the chosen function to
%                       approximate. If 0 and funapprox is 'expapprox' then
%                       next-period responses are calculated by equations solve and
%                       not just interpolated. If 1 and funapprox is 'resapprox', the
%                       guess of response variables is found with the new
%                       approximation structure
%
% [INTERP,X] = RECSSOLVEREE(INTERP,MODEL,S,X,...) returns the value of the response
% variables on the grid.
%
% [INTERP,X,Z] = RECSSOLVEREE(INTERP,MODEL,S,X,...) returns the value of the
% expectations variables on the grid.
%
% [INTERP,X,Z,F] = RECSSOLVEREE(INTERP,MODEL,S,X,...) returns the value of the
% equilibrium equations on the grid.
%
% [INTERP,X,Z,F,EXITFLAG] = RECSSOLVEREE(INTERP,MODEL,S,X,...) returns EXITFLAG,
% which describes the exit conditions. Possible values are
%    1 : RECSSOLVEREE converges to the REE
%    0 : Failure to converge
%
% See also RECSCHECK, RECSSIMUL, RECSSS.

% Copyright (C) 2011-2012 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
if nargin <=6, options = struct([]); end

defaultopt = struct(                           ...
    'display'           , 1                   ,...
    'eqsolver'          , 'lmmcp'             ,...
    'eqsolveroptions'   , struct([])          ,...
    'extrapolate'       , 1                   ,...
    'functional'        , 0                   ,...
    'loop_over_s'       , 0                   ,...
    'funapprox'         , 'resapprox-complete',...
    'reemethod'         , 'iter'              ,...
    'reesolver'         , 'sa'                ,...
    'reesolveroptions'  , struct([])          ,...
    'useapprox'         , 1);
warning('off','catstruct:DuplicatesFound')

options = catstruct(defaultopt,options);

extrapolate        = options.extrapolate;
funapprox          = lower(options.funapprox);
functional         = options.functional;
reemethod          = lower(options.reemethod);

% Extract fields of model
e      = model.e;
params = model.params;
w      = model.w;
if isa(model.func,'char')
  model.func = str2func(model.func);
elseif ~isa(model.func,'function_handle')
  error('model.func must be either a string or a function handle')
end
func = model.func;

if strcmp(reemethod,'1-step') && ...
  (strcmp(funapprox,'expapprox') || strcmp(funapprox,'resapprox'))
  warning('RECS:Switching2Iterative',...
          ['Solving the rational expectations problem is not implemented when ' ...
           'approximating this funtion. Switching to the iterative scheme.'])
end

% Identify variables dimensions
[n,m]  = size(x);
[~,d]  = size(s);
output = struct('F',1,'Js',0,'Jx',0,'Jsn',0,'Jxn',0,'hmult',0);
k      = length(w);               % number of shock values
z      = zeros(n,0);
p      = size(func('h',s(1,:),x(1,:),[],e(1,:),s(1,:),x(1,:),params,output),2);

X      = zeros(n,m,T);
S      = zeros(n,d,T);

X(:,:,T) = xT;

% Extract fields of interp
fspace     = interp.fspace;
interp.Phi = funbasx(fspace);
Phi        = interp.Phi;

% Check last period situation
if any(reshape(isnan(X(:,:,T)),n*m,1))
  [LB,UB] = func('b',s,[],[],[],[],[],params);
  LB(~isnan(X(:,:,T))) = X(~isnan(X(:,:,T)));
  UB(~isnan(X(:,:,T))) = X(~isnan(X(:,:,T)));
  z = zeros(n,p);
  c = x;
  optionsT = catstruct(options,struct('funapprox','expapprox'));
  X(:,:,T) = recsSolveEquilibriumFH(s,x,z,func,params,c,e,w,fspace,LB,UB,optionsT);
end

[interp.cX,X,exitflag] = recsSolveREEIterFiniteHorizon(interp,model,s,x,X,options);

if exitflag~=1
  warning('RECS:FailureREE','Failure to find a rational expectations equilibrium');
end

