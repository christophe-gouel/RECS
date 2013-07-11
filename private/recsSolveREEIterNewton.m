function [c,x,fval,exitflag] = recsSolveREEIterNewton(interp,model,s,x,c,options)
% RECSSOLVEREEITERNEWTON finds the REE of a model by doing a Newton iteration for the rational expectations step
%
% RECSSOLVEREEITERNEWTON is called by recsSolveREE. It is not meant to be called
% directly by the user.
%
% This function is not yet adapted to MCP problems. For these problems, we have
% to rely on numerical derivatives.
%
% See also RECSSOLVEREE, RECSSOLVEREEITER, RECSSOLVEREEFULL.

% Copyright (C) 2011-2013 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
extrapolate        = options.extrapolate;
functional         = options.functional;
funapprox          = lower(options.funapprox);
reesolver          = lower(options.reesolver);
reesolveroptions   = catstruct(struct('showiters'      ,options.display,...
                                      'DerivativeCheck','off'          ,...
                                      'Display'        ,'iter'         ,...
                                      'Jacobian'       , 'on'          ,...
                                      'lmeth'          , 3)         ,...
                               options.reesolveroptions);

NewtonMethod       = ~any(strcmpi(reesolver,{'krylov','sa'}));

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

[n,m]    = size(x);
z        = zeros(n,0);
[~,grid] = spblkdiag(zeros(m,m,n),[],0);

if strcmp(funapprox,'resapprox-complete'),  c = c(:,ixforward); end

%% Solve for the rational expectations equilibrium
[c,~,exitREE] = runeqsolver(@ResidualFunction,c(:),...
                            -inf(numel(c),1),inf(numel(c),1),...
                            reesolver,...
                            reesolveroptions);

exitflag = and(exitREE,exitEQ);

if strcmp(funapprox,'resapprox-complete'), c = funfitxy(fspace,Phi,x); end


%% Nested function
function [R,dRdc] = ResidualFunction(cc)
% RESIDUALFUNCTION Calculates the residual of the model with regards to rational expectations

  cc    = reshape(cc,n,[]);
  if functional, params{end} = cc; end

  [x,fval,exitEQ]  = recsSolveEquilibrium(s,x,z,b,f,g,h,params,cc,e,w,fspace,...
                                          ixforward,options);
  if nargout==1
    %% Without Jacobian
    R = recsResidual(s,x,h,params,cc,fspace,funapprox,Phi,ixforward,NewtonMethod);
  else
    %% With Jacobian
    [R,Rx,Rc] = recsResidual(s,x,h,params,cc,fspace,funapprox,Phi,ixforward,true);
    [~,Fx,Fc] = recsEquilibrium(x,s,z,b,f,g,h,params,grid,cc,e,w,fspace,...
                                funapprox,extrapolate,ixforward);
    Fc        = sparse(Fc);
    dRdc      = -Rx*(Fx\Fc)+Rc;
  end

end

end
