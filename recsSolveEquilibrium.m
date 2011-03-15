function [x,f] = recsSolveEquilibrium(s,x,z,func,params,eqsolver,c,e,w,fspace, ...
                                      method,eqsolveroptions,loop_over_s)
% RECSSOLVEEQUILIBRIUM Solves the system of equilibrium equations using x as starting values

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

if nargin < 13
  loop_over_s = 0;
  if nargin < 12
    eqsolveroptions = struct([]);
  end
end

[n,m]   = size(x);
[LB,UB] = func('b',s,[],[],[],[],[],params);

if loop_over_s % Solve separately for each point on the grid
  f                 = zeros(size(x));
  [~,grid]          = spblkdiag(zeros(m,m,1),[],0);
  for i=1:n
    [x(i,:),f(i,:)] = eqsolve(x(i,:),s(i,:),z(i,:),func,params,eqsolver,grid,c, ...
                              e,w,fspace,method,eqsolveroptions,LB(i,:),UB(i,:));
  end
else % Solve all the grid in one step
  [~,grid] = spblkdiag(zeros(m,m,n),[],0);
  [x,f]    = eqsolve(x,s,z,func,params,eqsolver,grid,c,e,w,fspace,method, ...
                     eqsolveroptions,LB,UB);
end

function [x,f] = eqsolve(x,s,z,func,params,eqsolver,grid,c,e,w,fspace,method, ...
                         eqsolveroptions,LB,UB)

[n,m]   = size(x);
x       = reshape(x',[n*m 1]);
LB      = reshape(LB',[n*m 1]);
UB      = reshape(UB',[n*m 1]);

switch lower(eqsolver)
 case 'fsolve'
  options = optimset('Display','off',...
                     'Jacobian','on');
  options = optimset(options,eqsolveroptions);
  [x,f,exitflag] = fsolve(@(X) recsEquilibrium(X,s,z,func,params,grid,c,e,w,...
                                               fspace,method),...
                          x,options);
  if exitflag~=1, disp('No convergence'); end
 case 'lmmcp'
  [x,f,exitflag] = lmmcp(@(X) recsEquilibrium(X,s,z,func,params,grid,c,e,w,...
                                              fspace,method),...
                         x,LB,UB,eqsolveroptions);
  if exitflag~=1, disp('No convergence'); end
 case 'ncpsolve'
  [x,f] = ncpsolve(@(X) ncpsolvetransform(X,@recsEquilibrium,s,z,func,...
                                          params,grid,c,e,w,fspace,method),...
                   LB,UB,x);
  f     = -f;
 case 'path'
  global par
  par   = {@recsEquilibrium,s,z,func,params,grid,c,e,w,fspace,method};
  [x,f] = pathmcp(x,LB,UB,'pathtransform');
  clear global par
end
x    = reshape(x,m,n)';
f    = reshape(f,m,n)';

return