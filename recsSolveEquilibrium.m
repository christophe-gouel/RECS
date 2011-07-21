function [x,f] = recsSolveEquilibrium(s,x,z,func,params,c,e,w,fspace,options)
% RECSSOLVEEQUILIBRIUM Solves the system of equilibrium equations using x as starting values
%
% RECSSOLVEEQUILIBRIUM is called by RECSSIMUL and RECSSOLVEREE. It is
% not meant to be called directly by the user.
%
% See also RECSSIMUL, RECSSOLVEREE.

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

eqsolver         = lower(options.eqsolver);
eqsolveroptions  = options.eqsolveroptions;
extrapolate      = options.extrapolate;
method           = lower(options.method);

[n,m]   = size(x);
[LB,UB] = func('b',s,[],[],[],[],[],params);

if options.loop_over_s % Solve separately for each point on the grid
  f                 = zeros(size(x));
  [~,grid]          = spblkdiag(zeros(m,m,1),[],0);
  for i=1:n
    [x(i,:),f(i,:)] = eqsolve(x(i,:),s(i,:),z(i,:),func,params,eqsolver,grid,c,e,w,...
                              fspace,method,eqsolveroptions,LB(i,:),UB(i,:),extrapolate);
  end
else % Solve all the grid in one step
  [~,grid] = spblkdiag(zeros(m,m,n),[],0);
  [x,f]    = eqsolve(x,s,z,func,params,eqsolver,grid,c,e,w,fspace,method, ...
                     eqsolveroptions,LB,UB,extrapolate);
end

function [x,f] = eqsolve(x,s,z,func,params,eqsolver,grid,c,e,w,fspace,method, ...
                         eqsolveroptions,LB,UB,extrapolate)

[n,m]   = size(x);
x       = reshape(x',[n*m 1]);
LB      = reshape(LB',[n*m 1]);
UB      = reshape(UB',[n*m 1]);

switch eqsolver
 case 'fsolve'
  options = optimset('Display','off',...
                     'Jacobian','on');
  options = optimset(options,eqsolveroptions);
  [x,f,exitflag] = fsolve(@(X) recsEquilibrium(X,s,z,func,params,grid,c,e,w,...
                                               fspace,method,extrapolate),...
                          x,options);
  if exitflag~=1, disp('No convergence'); end
 case 'lmmcp'
  [x,f,exitflag] = lmmcp(@(X) recsEquilibrium(X,s,z,func,params,grid,c,e,w,...
                                              fspace,method,extrapolate),...
                         x,LB,UB,eqsolveroptions);
  if exitflag~=1, disp('No convergence'); end
 case 'ncpsolve'
  [x,f] = ncpsolve(@(X) ncpsolvetransform(X,@recsEquilibrium,s,z,func,params,...
                                          grid,c,e,w,fspace,method,extrapolate),...
                   LB,UB,x);
  f     = -f;
 case 'path'
  global par
  par   = {@recsEquilibrium,s,z,func,params,grid,c,e,w,fspace,method,extrapolate};
  [x,f] = pathmcp(x,LB,UB,'pathtransform');
  clear global par
end
x    = reshape(x,m,n)';
f    = reshape(f,m,n)';

return