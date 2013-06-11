function [c,x,fval,exitflag] = recsSolveREEFull(interp,model,s,x,c,options)
% RECSSOLVEREEFULL finds the REE of a model as one single problem and not by iterating on two subproblems
%
% RECSSOLVEREEFULL is called by RECSSOLVEREE. It is not meant to be called directly
% by the user.
%
% See also RECSSOLVEEREE, RECSSOLVEEREEITER.

% Copyright (C) 2011-2013 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
eqsolver         = lower(options.eqsolver);
eqsolveroptions  = options.eqsolveroptions;
extrapolate      = options.extrapolate;
funapprox        = lower(options.funapprox);
functional       = options.functional;

b         = model.b;
e         = model.e;
f         = model.f;
g         = model.g;
h         = model.h;
ixforward = model.ixforward;
params    = model.params;
w         = model.w;

fspace = interp.fspace;
Phi    = interp.Phi;

[n,m]  = size(x);

[~,grid] = spblkdiag(zeros(m,m,n),[],0);
[LB,UB]  = model.b(s,params);

if strcmp(funapprox,'resapprox-complete')
  c = c(:,ixforward);
end

if strcmp(funapprox,'resapprox-complete') && all(isinf(LB(:))) && all(isinf(UB(:)))
  %% Reshape inputs
  X              = x;
  X(:,ixforward) = c;
  X              = reshape(X',n*m,1);
  B              = inf(n*m,1);

  %% Solve for the rational expectations equilibrium
  [X,F,exitflag] = runeqsolver(@FullCompactPb,X,-B,B,eqsolver,eqsolveroptions,...
                               s,b,f,g,h,params,grid,e,w,fspace,funapprox,Phi,...
                               m,functional,extrapolate,ixforward);

  %% Reshape outputs
  x              = reshape(X,m,n)';
  c              = x(:,ixforward);
  x(:,ixforward) = funeval(c,fspace,Phi);
  fval  = reshape(F,m,n)';
else
  %% Reshape inputs
  X        = [reshape(x',[n*m 1]); reshape(c',[],1)];
  LB       = [reshape(LB',[n*m 1]); -inf(n*size(c,2),1)];
  UB       = [reshape(UB',[n*m 1]); +inf(n*size(c,2),1)];

  %% Solve for the rational expectations equilibrium
  [X,G,exitflag] = runeqsolver(@FullPb,X,LB,UB,eqsolver,eqsolveroptions,...
                               s,b,f,g,h,params,grid,e,w,fspace,funapprox,Phi,...
                               m,functional,extrapolate,ixforward);

  %% Reshape outputs
  x     = reshape(X(1:n*m),m,n)';
  c     = reshape(X(n*m+1:end),[],n)';
  fval  = reshape(G(1:n*m),m,n)';
end

if strcmp(funapprox,'resapprox-complete')
  c = funfitxy(fspace,Phi,x);
end

function [G,J] = FullPb(X,s,b,f,g,h,params,grid,e,w,fspace,funapprox,Phi,m,functional,extrapolate,ixforward)
% FULLPB evaluates the equations and Jacobian of the complete rational expectations problem

%% Initialization
n     = size(s,1);
x     = X(1:m*n);
x     = reshape(x,m,n)';
c     = reshape(X(m*n+1:end),[],n)';
if functional, params{end} = c; end

%% Evaluate equations and Jacobian
if nargout==2
  %% With Jacobian
  [F,Fx,Fc] = recsEquilibrium(reshape(x',[n*m 1]),s,zeros(n,0),b,f,g,h,params,...
                              grid,c,e,w,fspace,funapprox,extrapolate,ixforward);
  [R,Rx,Rc] = recsResidual(s,x,h,params,c,fspace,funapprox,Phi,ixforward);
  J = [Fx Fc;
       Rx Rc];
else
  %% Without Jacobian
  F = recsEquilibrium(reshape(x',[n*m 1]),s,zeros(n,0),b,f,g,h,params,...
                      grid,c,e,w,fspace,funapprox,extrapolate,ixforward);
  R = recsResidual(s,x,h,params,c,fspace,funapprox,Phi,ixforward);
end

% Concatenate equilibrium equations and rational expectations residual
G = [F; R];


function [F,J] = FullCompactPb(X,s,b,f,g,h,params,grid,e,w,fspace,funapprox,Phi,m,functional,extrapolate,ixforward)
% FULLCOMPACTPB evaluates the equations and Jacobian of the complete rational expectations problem as a single set of equations (equilibrium and residual equation are confounded)

%% Initialization
n              = size(s,1);
x              = reshape(X,m,n)';
c              = x(:,ixforward);
if functional, params{end} = c; end
x(:,ixforward) = funeval(c,fspace,Phi);
indx           = 1:n*m;
indxforward    = repmat(ixforward,1,n);
indxstatic     = ~indxforward;

%% Evaluate equations and Jacobian
if nargout==2
  %% With Jacobian
  [F,Fx,Fc] = recsEquilibrium(reshape(x',[n*m 1]),s,zeros(n,0),b,f,g,h,params,...
                              grid,c,e,w,fspace,funapprox,extrapolate,ixforward);
  [~,~,Rc]  = recsResidual(s,x,h,params,c,fspace,funapprox,Phi,ixforward);
  Fc = sparse(Fc);

  J(:,[indx(indxstatic) indx(indxforward)]) = [Fx(:,indxstatic) Fx(:,indxforward)*Rc+Fc];
  
else
  %% Without Jacobian
  F = recsEquilibrium(reshape(x',[n*m 1]),s,zeros(n,0),b,f,g,h,params,...
                      grid,c,e,w,fspace,funapprox,extrapolate,ixforward);
end


