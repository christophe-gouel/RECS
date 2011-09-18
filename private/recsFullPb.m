function [G,J] = recsFullPb(X,s,func,params,grid,e,w,fspace,funapprox,Phi,m,functional,extrapolate)
% RECSFULLPB evaluates the equations and Jacobian of the complete rational expectations problem
%
% RECSFULLPB is called by recsSolveREEFull. It is not meant to be called directly
% by the user.
%
% See also RECSSOLVEREEFULL.

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

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
  [R,Rx,Rc] = recsResidual(s,x,func,params,c,fspace,funapprox,Phi);

  J = [Fx Fc;
       Rx Rc];
else
  %% Without Jacobian
  F = recsEquilibrium(reshape(x',[n*m 1]),s,zeros(n,0),...
                       func,params,grid,c,e,w,fspace,funapprox,extrapolate);
  R = recsResidual(s,x,func,params,c,fspace,funapprox,Phi);
end

% Concatenate equilibrium equations and rational expectations residual
G = [F; R];
