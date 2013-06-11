function [interp,X,exitflag] = recsSolveREEFiniteHorizon(interp,model,s,x,xT,T,options)
% RECSSOLVEREEFINITEHORIZON finds the finite-horizon rational expectations equilibrium of a model

% Copyright (C) 2011-2013 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
defaultopt = struct(                                        ...
    'eqsolver'          , 'lmmcp'                          ,...
    'eqsolveroptions'   , struct('DerivativeCheck', 'off' ,...
                                 'Jacobian'       , 'on')  ,...
    'extrapolate'       , 1                                ,...
    'loop_over_s'       , 0                                ,...
    'funapprox'         , 'resapprox-complete');
if nargin <=6
  options = defaultopt;
else
  warning('off','catstruct:DuplicatesFound')
  if isfield(options,'eqsolveroptions')
    options.eqsolveroptions = catstruct(defaultopt.eqsolveroptions,options.eqsolveroptions);
  end
  options = catstruct(options,struct('funapprox','resapprox-complete'));
  options = catstruct(defaultopt,options);
end

% Extract fields of model
b         = model.b;
e         = model.e;
f         = model.f;
g         = model.g;
h         = model.h;
ixforward = model.ixforward;
params    = model.params;
w         = model.w;

% Identify variables dimensions
n      = size(x,1);
z      = zeros(n,0);
[m,p]  = model.dim{2:3};

X      = zeros(n,m,T);
cX     = zeros(n,m,T);

% Extract fields of interp
fspace     = interp.fspace;
interp.Phi = funbasx(fspace);
Phi        = interp.Phi;

%% Check last period situation
if any(isnan(xT(:)))
  [LB,UB] = b(s,params);
  LB(~isnan(xT)) = xT(~isnan(xT));
  UB(~isnan(xT)) = xT(~isnan(xT));
  z = zeros(n,p);
  c = x;
  xinit = x;
  xinit(~isnan(xT)) = xT(~isnan(xT));
  optionsT = catstruct(options,struct('funapprox','expapprox'));
  X(:,:,T) = recsSolveEquilibrium(s,xinit,z,b,f,g,h,params,c,e,w,fspace,...
                                  ixforward,optionsT,LB,UB);
else
  X(:,:,T) = xT;
end

%% Solve for the rational expectations equilibrium
for t=T-1:-1:1
  cX(:,:,t+1)           = funfitxy(fspace,Phi,X(:,:,t+1));
  [X(:,:,t),~,exitflag] = recsSolveEquilibrium(s,X(:,:,t+1),z,b,f,g,h,params,...
                                               cX(:,ixforward,t+1),e,w,fspace,...
                                               ixforward,options);
  if exitflag~=1
    warning('RECS:FailureREE','Failure to find a rational expectations equilibrium');
  end
end
cX(:,:,1)               = funfitxy(fspace,Phi,X(:,:,1));

interp.cX = cX;


