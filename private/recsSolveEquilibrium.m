function [x,f] = recsSolveEquilibrium(s,x,z,func,params,c,e,w,fspace,options)
% RECSSOLVEEQUILIBRIUM Solves the system of equilibrium equations using x as starting values
%
% RECSSOLVEEQUILIBRIUM is called by recsSimul and recsSolveREE. It is
% not meant to be called directly by the user.
%
% See also RECSSIMUL, RECSSOLVEREE.

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
eqsolver         = lower(options.eqsolver);
eqsolveroptions  = options.eqsolveroptions;
extrapolate      = options.extrapolate;
funapprox        = lower(options.funapprox);

[n,m]   = size(x);
[LB,UB] = func('b',s,[],[],[],[],[],params);

%% Solve equilibrium equations on grid points
if options.loop_over_s
  %% Solve separately for each point on the grid
  f                 = zeros(size(x));
  [~,grid]          = spblkdiag(zeros(m,m,1),[],0);
  for i=1:n
    [x(i,:),f(i,:)] = eqsolve(x(i,:),s(i,:),z(i,:),func,params,eqsolver,grid,c,e,w,...
                              fspace,funapprox,eqsolveroptions,LB(i,:),UB(i,:),extrapolate);
  end
else
  %% Solve all the grid in one step
  [~,grid] = spblkdiag(zeros(m,m,n),[],0);
  [x,f]    = eqsolve(x,s,z,func,params,eqsolver,grid,c,e,w,fspace,funapprox, ...
                     eqsolveroptions,LB,UB,extrapolate);
end

function [x,f] = eqsolve(x,s,z,func,params,eqsolver,grid,c,e,w,fspace,funapprox, ...
                         eqsolveroptions,LB,UB,extrapolate)
%% Eqsolve

[n,m]   = size(x);
x       = reshape(x',[n*m 1]);
LB      = reshape(LB',[n*m 1]);
UB      = reshape(UB',[n*m 1]);

[x,f,exitflag] = runeqsolver(@recsEquilibrium,x,LB,UB,eqsolver,eqsolveroptions,s,...
                             z,func,params,grid,c,e,w,fspace,funapprox,extrapolate);

if exitflag~=1, disp('No convergence'); end

x    = reshape(x,m,n)';
if ~isempty(f), f = reshape(f,m,n)'; end

return