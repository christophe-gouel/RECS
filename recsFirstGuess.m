function [interp,x,z,exitflag,output] = recsFirstGuess(interp,model,s,sss,xss,options)
% RECSFIRSTGUESS finds a first guess using the perfect foresight solution
%
% RECSFIRSTGUESS tries to find a first-guess for the rational expectations model
% by solving the corresponding perfect foresight solution on all the grid points
% of the state variables. By default, it considers that the model goes back to
% its deterministic steady state in 50 periods.
%
% INTERP = RECSFIRSTGUESS(INTERP,MODEL,S,SSS,XSS) uses the interpolation structure
% defined in the structure INTERP to fit the perfect foresight solution of the
% model defined in the structure MODEL. The grid of state variables is provided in
% the n-by-d matrix S. A first guess for the steady state of the model is provided
% in SSS and XSS for state and response variables. RECSFIRSTGUESS returns an
% interpolation structure, INTERP, containing the first guess.
% INTERP is a structure, which has to include the following fields:
%    fspace       : a definition structure for the interpolation family (created
%                   by the function fundef)
% MODEL is an object created by recsmodel.
%
% INTERP = RECSFIRSTGUESS(INTERP,MODEL,S,SSS,XSS,OPTIONS) solves the problem
% with the parameters defined by the structure OPTIONS. The fields of the
% structure are
%    eqsolver         : 'fsolve', 'lmmcp' (default), 'ncpsolve' or 'path'
%    eqsolveroptions  : options structure to be passed to eqsolver (default:
%                       empty structure)
%    fgmethod         : 'auto' (default), 'perturbation', 'perfect-foresight' or
%                       'steady-state'
%    T                : integer defining the time horizon at which the model is
%                       supposed to converge to its steady state (default: 50)
%
% [INTERP,X] = RECSFIRSTGUESS(INTERP,MODEL,S,SSS,XSS,...) returns X, n-by-m matrix,
% containing the value of the response variables in the first period.
%
% [INTERP,X,Z] = RECSFIRSTGUESS(INTERP,MODEL,S,SSS,XSS,...) returns Z, n-by-p matrix,
% containing the value of the expectations variables in the first period.
%
% [INTERP,X,Z,EXITFLAG] = RECSFIRSTGUESS(INTERP,MODEL,S,SSS,XSS,...)
%
% [INTERP,X,Z,EXITFLAG,OUTPUT] = RECSFIRSTGUESS(INTERP,MODEL,S,SSS,XSS,...)
%
% See also RECSSOLVEDETERMINISTICPB, RECSSOLVELOCAL, RECSSOLVEREE, RECSSS.

% Copyright (C) 2011-2013 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
% Get s from interp structure
if nargin<=2 || isempty(s), s = interp.s; end

if nargin <=3 || isempty(sss)
  if isfield(model,'sss') || (isobject(model) && ~isempty(model.sss))
    sss = model.sss;
  else
  sss = [];
  end
end
if nargin <=4 || isempty(xss)
  if isfield(model,'xss') || (isobject(model) && ~isempty(model.xss))
    xss = model.xss;
  else
    xss = [];
  end
end
defaultopt = struct('fgmethod','auto',...
                    'T'       ,50);
if nargin <=5
  options = defaultopt;
else
  warning('off','catstruct:DuplicatesFound')
  options = catstruct(defaultopt,options);
end

fgmethod = options.fgmethod;
T        = options.T;

n        = size(s,1);

p        = model.dim{3};
params   = model.params;

%% Solve for the deterministic steady state
[sss,xss,zss] = recsSS(model,sss,xss,catstruct(options,struct('display',0)));

%% Find first-guess
[LB,UB]    = model.functions.b(s,params);
if strcmpi(fgmethod,'auto') && all(isinf([LB(:); UB(:)]))
  fgmethod = 'perturbation';
elseif strcmpi(fgmethod,'auto')
  fgmethod = 'perfect-foresight';
end

switch fgmethod
  case 'perturbation'
    %% First-order perturbation
    model    = recsSolveLocal(model);
    x        = min(max(model.LinearSolution.X(s),LB),UB);
    z        = model.LinearSolution.Z(s);
    exitflag = 1;
    output   = [];

  case 'perfect-foresight'
    %% Solve the perfect foresight problem on each point of the grid
    x        = zeros(n,size(xss,2));
    z        = zeros(n,size(zss,2));
    exitflag = zeros(n,1);
    N        = zeros(n,1);

    parfor i=1:n
      [X,~,Z,~,exitflag(i),N(i)] = recsSolveDeterministicPb(model,s(i,:),...
                                                        T,xss,zss,sss,options);
      x(i,:) = X(1,:);
      z(i,:) = Z(1,:);
    end
    output = struct('exitflag',exitflag,...
                    'N'       ,N);
    exitflag = ~any(exitflag~=1);

  case 'steady-state'
    %% Response variables equal to their steady-state values
    x        = min(max(repmat(xss,n,1),LB),UB);
    z        = NaN(n,p);
    exitflag = 1;
    output   = [];

end % switch fgmethod

%% Prepare output
interp.cx = funfitxy(interp.fspace,interp.Phi,x);
interp.cz = funfitxy(interp.fspace,interp.Phi,z);
interp.x  = x;
interp.z  = z;

% Check if it is possible to approximate the expectations function
if strcmp(model.infos.model_type,'fgh1')
  hv        = model.functions.h(zeros(n,0),[],[],s,x,params);
  interp.ch = funfitxy(interp.fspace,interp.Phi,hv);
end

