function [interp,X,exitflag] = recsSolveREEFiniteHorizon(interp,model,s,x,xT,T,options)
% RECSSOLVEREEFINITEHORIZON finds the finite-horizon rational expectations equilibrium of a model

% Copyright (C) 2011-2012 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
defaultopt = struct(                                      ...
    'eqsolver'          , 'lmmcp'                        ,...
    'eqsolveroptions'   , struct('DerivativeCheck','off'),...
    'extrapolate'       , 1                              ,...
    'loop_over_s'       , 0                              ,...
    'funapprox'         , 'resapprox-complete');
if nargin <=6
  options = defaultopt;
else
  warning('off','catstruct:DuplicatesFound')
  options = catstruct(options,struct('funapprox','resapprox-complete'));
  options = catstruct(defaultopt,options);
end

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

% Identify variables dimensions
[n,m]  = size(x);
output = struct('F',1,'Js',0,'Jx',0,'Jsn',0,'Jxn',0,'hmult',0);
z      = zeros(n,0);
p      = size(func('h',s(1,:),x(1,:),[],e(1,:),s(1,:),x(1,:),params,output),2);

X      = zeros(n,m,T);
cX     = zeros(n,m,T);

% Extract fields of interp
fspace     = interp.fspace;
interp.Phi = funbasx(fspace);
Phi        = interp.Phi;

%% Check last period situation
if any(isnan(xT(:)))
  [LB,UB] = func('b',s,[],[],[],[],[],params);
  LB(~isnan(xT)) = xT(~isnan(xT));
  UB(~isnan(xT)) = xT(~isnan(xT));
  z = zeros(n,p);
  c = x;
  xinit = x;
  xinit(~isnan(xT)) = xT(~isnan(xT));
  optionsT = catstruct(options,struct('funapprox','expapprox'));
  X(:,:,T) = recsSolveEquilibriumFH(s,xinit,z,func,params,c,e,w,fspace,LB,UB,optionsT);
else
  X(:,:,T) = xT;
end

%% Solve for the rational expectations equilibrium
for t=T-1:-1:1
  cX(:,:,t+1)           = funfitxy(fspace,Phi,X(:,:,t+1));
  [X(:,:,t),~,exitflag] = recsSolveEquilibrium(s,X(:,:,t+1),z,func,params,...
                                               cX(:,:,t+1),e,w,fspace,options);
  if exitflag~=1
    warning('RECS:FailureREE','Failure to find a rational expectations equilibrium');
  end
end
cX(:,:,1)               = funfitxy(fspace,Phi,X(:,:,1));

interp.cX = cX;


