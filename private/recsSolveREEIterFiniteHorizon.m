function [cX,X,exitflag] = recsSolveREEIterFiniteHorizon(interp,model,s,x,X,options)
% RECSSOLVEREEITER Finds the REE of a model by iteration between equilibrium equations and rational expectations
%
% RECSSOLVEREEITER is called by recsSolveREE. It is not meant to be called directly
% by the user.
%
% See also RECSSOLVEEREE, RECSSOLVEREEITERNEWTON, RECSSOLVEEREEFULL.

% Copyright (C) 2011-2012 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
funapprox          = lower(options.funapprox);
useapprox          = options.useapprox;

e      = model.e;
func   = model.func;
params = model.params;
w      = model.w;

fspace = interp.fspace;
Phi    = interp.Phi;

[n,~,T] = size(X);
z       = zeros(n,0);

cX      = zeros(size(X));

%% Solve for the rational expectations equilibrium
for t=T-1:-1:1
  cX(:,:,t+1)         = funfitxy(fspace,Phi,X(:,:,t+1));
  [X(:,:,t),~,exitEQ] = recsSolveEquilibrium(s,X(:,:,t+1),z,func,params,cX(:,:,t+1),e,w,fspace,options);
end

if exitEQ==1
  exitflag = 1;
else
  exitflag = 0;
end

