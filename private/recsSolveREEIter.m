function [c,x,z,fval,exitflag] = recsSolveREEIter(interp,model,s,x,c,options)
% RECSSOLVEREEITER the REE of a model by iteration between equilibrium equations and rational expectations
%
% RECSSOLVEREEITER is called by recsSolveREE. It is not meant to be called
% directly by the user.
%
% See also RECSSOLVEREE, RECSSOLVEREEFULL.

% Copyright (C) 2011-2013 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
explicit           = options.explicit;
extrapolate        = options.extrapolate;
funapprox          = lower(options.funapprox);
functional         = options.functional;
reesolver          = lower(options.reesolver);
reesolveroptions   = catstruct(struct('showiters'      ,options.display,...
                                      'Diagnostics'    , 'off'         ,...
                                      'DerivativeCheck','off'          ,...
                                      'Display'        ,'iter'         ,...
                                      'Jacobian'       , 'on'          ,...
                                      'lmeth'          , 3)         ,...
                               options.reesolveroptions);
useapprox          = options.useapprox;

NewtonMethod       = ~any(strcmpi(reesolver,{'sa','krylov','mixed'}));

b         = model.b;
e         = model.e;
f         = model.f;
g         = model.g;
h         = model.h;
ixforward = model.ixforward;
params    = model.params;
w         = model.w;
[m,p]     = model.dim{2:3};

fspace = interp.fspace;
Phi    = interp.Phi;

n   = size(x,1);
vec = @(X) X(:);

%% Precalculations
if ~(explicit || strcmp(funapprox,'expapprox'))
  z        = zeros(n,0);
  [~,grid] = spblkdiag(zeros(m,m,n),[],0);
else
  k       = length(w);
  xnext   = zeros(n*k,m);
  [LB,UB] = b(s,params);
  ind     = (1:n);
  ind     = ind(ones(1,k),:);
  ss      = s(ind,:);
  ee      = e(repmat(1:k,1,n),:);
  exitEQ  = 1;
  fval    = zeros(size(x));
end

if strcmpi(funapprox,'resapprox'),  c = c(:,ixforward); end

%% Solve for the rational expectations equilibrium
[c,~,exitREE] = runeqsolver(@ResidualFunction,vec(c'),...
                            -inf(numel(c),1),inf(numel(c),1),...
                            reesolver,...
                            reesolveroptions);

exitflag = and(exitREE,exitEQ);

if strcmpi(funapprox,'resapprox'), c = funfitxy(fspace,Phi,x);
else                               c = reshape(c,[],n)';
end


%% Nested function
function [R,dRdc] = ResidualFunction(cc)
% RESIDUALFUNCTION Calculates the residual of the model with regards to rational expectations

  cc    = reshape(cc,[],n)';
  if functional, params{end} = cc; end

  if ~(explicit || strcmpi(funapprox,'expapprox'))
    %% Non-explicit models with response variable or expectations function approximation
    [x,fval,exitEQ]  = recsSolveEquilibrium(s,x,z,b,f,g,h,params,cc,e,w,fspace,...
                                            ixforward,options);
    if nargout==1
      %% Without Jacobian
      R = recsResidual(s,x,h,params,cc,fspace,funapprox,Phi,ixforward,NewtonMethod);
    else
      %% With Jacobian
      [R,Rx,Rc] = recsResidual(s,x,h,params,cc,fspace,funapprox,Phi,ixforward,true);
      [~,Fx,Fc] = recsEquilibrium(vec(x'),s,z,b,f,g,h,params,grid,cc,e,w,...
                                  fspace,funapprox,extrapolate,ixforward);
      Fc        = sparse(Fc);
      dRdc      = -Rx*(Fx\Fc)+Rc;
    end

  elseif strcmpi(funapprox,'expapprox')
    %% Expectations approximation (PEA)

    % Calculation of z by interpolation
    z     = funeval(cc,fspace,Phi);

    % Calculation of x
    [x,fval,exitEQ] = recsSolveEquilibrium(s,x,z,b,f,g,h,params,cc,e,w,...
                                           fspace,ixforward,options);

    % Calculation of snext
    xx      = x(ind,:);
    snext   = g(ss,xx,ee,params);

    % Calculation of xnext
    if extrapolate>=1, snextinterp = snext;
    else
      snextinterp = max(min(snext,fspace.b(ones(n*k,1),:)),...
                        fspace.a(ones(n*k,1),:));
    end % extrapolate
    [LBnext,UBnext] = b(snext,params);
    if useapprox % xnext calculated by interpolation
      xnext(:,ixforward) = min(max(funeval(funfitxy(fspace,Phi,x(:,ixforward)),...
                                           fspace,snextinterp),...
                                   LBnext(:,ixforward)),UBnext(:,ixforward));
    else  % xnext calculated by equation solve
      xnext = min(max(funeval(funfitxy(fspace,Phi,x),fspace,snextinterp),...
                      LBnext),UBnext);
      xnext = recsSolveEquilibrium(snext,xnext,...
                                   funeval(cc,fspace,snextinterp),...
                                   b,f,g,h,params,cc,e,w,fspace,...
                                   ixforward,options);
    end

    % Calculation of z
    hv    = h(ss,xx,ee,snext,xnext,params);
    z     = reshape(w'*reshape(hv,k,n*p),n,p);

    % Prepare output
    R     = vec((funfitxy(fspace,Phi,z)-cc)');

  else
    %% Explicit models

    % Calculation of snext
    xx     = x(ind,:);
    snext  = g(ss,xx,ee,params);

    % Calculation of xnext
    if extrapolate>=1, snextinterp = snext;
    else
      snextinterp = max(min(snext,fspace.b(ones(n*k,1),:)),...
                        fspace.a(ones(n*k,1),:));
    end % extrapolate
    [LBnext,UBnext]    = b(snext,params);
    xnext(:,ixforward) = min(max(funeval(cc,fspace,snextinterp),...
                                 LBnext(:,ixforward)),UBnext(:,ixforward));

    % Calculation of z
    hv     = h(ss,xx,ee,snext,xnext,params);
    z      = reshape(w'*reshape(hv,k,n*p),n,p);

    % Calculation of x by explicit formula
    x      = min(max(f(s,[],z,params),LB),UB);

    % Prepare output
    R      = vec((funfitxy(fspace,Phi,x(:,ixforward))-cc)');

  end % if ~(explicit || strcmpi(funapprox,'expapprox'))

end

end
