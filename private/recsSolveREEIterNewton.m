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
explicit           = options.explicit;
extrapolate        = options.extrapolate;
funapprox          = lower(options.funapprox);
functional         = options.functional;
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
[m,p]     = model.dim{2:3};

fspace = interp.fspace;
Phi    = interp.Phi;

n = size(x,1);

%% Precalculations
if ~explicit
  z        = zeros(n,0);
  [~,grid] = spblkdiag(zeros(m,m,n),[],0);
else
  k       = length(w);
  xnext   = zeros(n*k,m);
  output  = struct('F',1,'Js',0,'Jx',0,'Jz',0,'Jsn',0,'Jxn',0,'hmult',1);
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
[c,~,exitREE] = runeqsolver(@ResidualFunction,c(:),...
                            -inf(numel(c),1),inf(numel(c),1),...
                            reesolver,...
                            reesolveroptions);

exitflag = and(exitREE,exitEQ);

if strcmpi(funapprox,'resapprox'), c = funfitxy(fspace,Phi,x); end


%% Nested function
function [R,dRdc] = ResidualFunction(cc)
% RESIDUALFUNCTION Calculates the residual of the model with regards to rational expectations

  cc    = reshape(cc,n,[]);
  if functional, params{end} = cc; end

  if ~explicit
    %% Non-explicit models
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
  
  else
    %% Explicit models
    
    % Calculation of snext
    xx     = x(ind,:);
    snext  = g(ss,xx,ee,params,output);

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
    if nargout(h)<6
      hv                 = h(ss,xx,ee,snext,xnext,params,output);
    else
      [hv,~,~,~,~,hmult] = h(ss,xx,ee,snext,xnext,params,output);
      hv                 = hv.*hmult;
    end
    z      = reshape(w'*reshape(hv,k,n*p),n,p);

    % Calculation of x by explicit formula
    x      = min(max(f(s,[],z,params,output),LB),UB);

    % Prepare output
    R      = funfitxy(fspace,Phi,x(:,ixforward))-cc;
    R      = R(:);

  end % if ~explicit

end

end
