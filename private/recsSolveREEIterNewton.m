function [c,x,f,exitflag] = recsSolveREEIterNewton(interp,model,s,x,c,options)
% RECSSOLVEREEITERFULL 
  
% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
extrapolate        = options.extrapolate;
functional         = options.functional;
funapprox          = lower(options.funapprox);
reesolveroptions   = catstruct(struct('showiters' , options.display,...
                                      'atol'      , sqrt(eps),...
                                      'lmeth'     , 3,...
                                      'rtol'      , eps),...
                               options.reesolveroptions);
useapprox          = options.useapprox;

e      = model.e;
func   = model.func;
params = model.params;
w      = model.w;

fspace = interp.fspace;
Phi    = interp.Phi;

[n,m]    = size(x);
z        = zeros(n,0);
[~,grid] = spblkdiag(zeros(m,m,n),[],0);

%% Solve for the rational expectations equilibrium
if options.display==1
  reesolveroptions = optimset('display','iter-detailed',...
                              'Diagnostics','on',...
                              'Jacobian','on',...
                              'DerivativeCheck','on');
end
[c,~,exitflag] = fsolve(@ResidualFunction, reshape(c',[],1), reesolveroptions);

c = reshape(c,[],n)';
c = c(:);

%%
function [R,dRdc] = ResidualFunction(cc)
% RESIDUALFUNCTION Calculates the residual of the model with regards to rational expectations

  cc    = reshape(cc,[],n)';
  if functional, params{end} = cc; end
    
  if useapprox && strcmp(funapprox,'resapprox-complete') % x calculated by interpolation
    [LB,UB] = func('b',s,[],[],[],[],[],params);
    x       = min(max(funeval(cc,fspace,Phi),LB),UB);
  end % if not previous x is used

  [x,f]  = recsSolveEquilibrium(s,x,z,func,params,cc,e,w,fspace,options);
  if nargout==1
    R = recsResidual(s,x,func,params,cc,fspace,funapprox,Phi);
  else
    [R,Rx,Rc] = recsResidual(s,x,func,params,cc,fspace,funapprox,Phi);
    [~,Fx,Fc] = recsEquilibrium(x,s,z,func,params,grid,cc,e,w,fspace,funapprox,extrapolate);
    dRdc = -Rx*(Fx\Fc)+Rc;
  end

end

end
