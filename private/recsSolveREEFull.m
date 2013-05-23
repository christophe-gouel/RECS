function [c,x,f,exitflag] = recsSolveREEFull(interp,model,s,x,c,options)
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

e      = model.e;
h      = model.h;
func   = model.func;
params = model.params;
w      = model.w;

fspace = interp.fspace;
Phi    = interp.Phi;

[n,m]  = size(x);

[~,grid] = spblkdiag(zeros(m,m,n),[],0);
[LB,UB]  = model.b(s,params);
if strcmp(funapprox,'resapprox-complete') && all(isinf(LB(:))) && all(isinf(UB(:)))
  %% Reshape inputs
  C        = reshape(c',n*m,1);
  B        = inf(n*m,1);

  %% Solve for the rational expectations equilibrium
  [C,F,exitflag] = runeqsolver(@FullCompactPb,C,-B,B,eqsolver,eqsolveroptions,...
                               s,func,params,grid,e,w,fspace,funapprox,Phi,...
                               m,functional,extrapolate);

  %% Reshape outputs
  c     = reshape(C,m,n)';
  x     = funeval(c,fspace,Phi);
  f     = reshape(F,m,n)';
else
  %% Reshape inputs
  X        = [reshape(x',[n*m 1]); reshape(c',[],1)];
  LB       = [reshape(LB',[n*m 1]); -inf(n*size(c,2),1)];
  UB       = [reshape(UB',[n*m 1]); +inf(n*size(c,2),1)];

  %% Solve for the rational expectations equilibrium
  [X,G,exitflag] = runeqsolver(@FullPb,X,LB,UB,eqsolver,eqsolveroptions,...
                               s,func,params,grid,e,w,fspace,funapprox,Phi,...
                               m,functional,extrapolate);
  
  %% Reshape outputs
  x     = reshape(X(1:n*m),m,n)';
  c     = reshape(X(n*m+1:end),[],n)';
  f     = reshape(G(1:n*m),m,n)';
end


function [G,J] = FullPb(X,s,func,params,grid,e,w,fspace,funapprox,Phi,m,functional,extrapolate)
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
  [F,Fx,Fc] = recsEquilibrium(reshape(x',[n*m 1]),s,zeros(n,0),func,params,...
                              grid,c,e,w,fspace,funapprox,extrapolate);
  [R,Rx,Rc] = recsResidual(s,x,h,params,c,fspace,funapprox,Phi);
  J = [Fx Fc;
       Rx Rc];
else
  %% Without Jacobian
  F = recsEquilibrium(reshape(x',[n*m 1]),s,zeros(n,0),...
                       func,params,grid,c,e,w,fspace,funapprox,extrapolate);
  R = recsResidual(s,x,h,params,c,fspace,funapprox,Phi);
end

% Concatenate equilibrium equations and rational expectations residual
G = [F; R];


function [F,J] = FullCompactPb(C,s,func,params,grid,e,w,fspace,funapprox,Phi,m,functional,extrapolate)
% FULLCOMPACTPB evaluates the equations and Jacobian of the complete rational expectations problem as a single set of equations (equilibrium and residual equation are confounded)

%% Initialization
n     = size(s,1);
c     = reshape(C,m,n)';
if functional, params{end} = c; end
x     = funeval(c,fspace,Phi);

%% Evaluate equations and Jacobian
if nargout==2
  %% With Jacobian
  [F,Fx,Fc] = recsEquilibrium(reshape(x',[n*m 1]),s,zeros(n,0),func,params,...
                              grid,c,e,w,fspace,funapprox,extrapolate);
  [~,~,Rc]  = recsResidual(s,x,h,params,c,fspace,funapprox,Phi);
  
  Fc = sparse(Fc);
  J  = Fx*Rc+Fc;

else
  %% Without Jacobian
  F = recsEquilibrium(reshape(x',[n*m 1]),s,zeros(n,0),...
                       func,params,grid,c,e,w,fspace,funapprox,extrapolate);
end


