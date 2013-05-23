function [x,f,exitflag] = recsSolveEquilibrium(s,x,z,func,params,c,e,w,fspace,options,LB,UB)
% RECSSOLVEEQUILIBRIUM Solves the system of equilibrium equations using x as starting values
%
% RECSSOLVEEQUILIBRIUM is called by recsSimul, recsSolveREEFiniteHorizon,
% recsSolveREEIter and recsSolveREEIterNewton. It is not meant to be called
% directly by the user.
%
% See also RECSSIMUL, RECSSOLVEREE.

% Copyright (C) 2011-2013 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
eqsolver         = lower(options.eqsolver);
eqsolveroptions  = options.eqsolveroptions;
extrapolate      = options.extrapolate;
funapprox        = lower(options.funapprox);
loop_over_s      = options.loop_over_s;

[n,m]            = size(x);
if nargin<=10, [LB,UB] = func('b',s,[],[],[],[],[],params); end

%% Solve equilibrium equations on grid points
if loop_over_s
  if loop_over_s==1
    %% Solve separately for each point on the grid
    f            = zeros(size(x));
    [~,grid]     = spblkdiag(zeros(m,m,1),[],0);
    exitflag     = zeros(n,1);
    parfor i=1:n
      [x(i,:),f(i,:),exitflag(i)] = eqsolve(x(i,:),s(i,:),z(i,:),...
                                            func,params,eqsolver,...
                                            grid,c,e,w,fspace,funapprox,...
                                            eqsolveroptions,...
                                            LB(i,:),UB(i,:),extrapolate);
    end
  else
    %% Solve separately for blocks of grid points
    sizeBlocks   = fix(n/loop_over_s)*ones(loop_over_s,1);
    sizeBlocks(1:rem(n,loop_over_s)) = sizeBlocks(1:rem(n,loop_over_s))+1;
    exitflag     = zeros(loop_over_s,1);
    s            = mat2cell(s ,sizeBlocks);
    x            = mat2cell(x ,sizeBlocks);
    z            = mat2cell(z ,sizeBlocks);
    LB           = mat2cell(LB,sizeBlocks);
    UB           = mat2cell(UB,sizeBlocks);
    f            = cell(loop_over_s,1);
    parfor i=1:loop_over_s
      [~,grid]   = spblkdiag(zeros(m,m,sizeBlocks(i)),[],0);
      [x{i},f{i},exitflag(i)] = eqsolve(x{i},s{i},z{i},...
                                        func,params,eqsolver,...
                                        grid,c,e,w,fspace,funapprox,...
                                        eqsolveroptions,...
                                        LB{i},UB{i},extrapolate);
    end
    x            = cell2mat(x);
    f            = cell2mat(f);
  end % if loop_over_s==1
  exitflag       = all(exitflag);    
else
  %% Solve all the grid in one step
  [~,grid]       = spblkdiag(zeros(m,m,n),[],0);
  [x,f,exitflag] = eqsolve(x,s,z,func,params,eqsolver,grid,c,e,w,fspace,...
                           funapprox,eqsolveroptions,LB,UB,extrapolate);
end % if loop_over_s

function [x,f,exitflag] = eqsolve(x,s,z,func,params,eqsolver,grid,c,e,w,fspace,...
                                  funapprox,eqsolveroptions,LB,UB,extrapolate)
%% Eqsolve

[n,m]          = size(x);
x              = reshape(x',[n*m 1]);
LB             = reshape(LB',[n*m 1]);
UB             = reshape(UB',[n*m 1]);

[x,f,exitflag] = runeqsolver(@recsEquilibrium,x,LB,UB,eqsolver,eqsolveroptions,s,...
                             z,func,params,grid,c,e,w,fspace,funapprox,extrapolate);

if exitflag~=1, disp('No convergence'); end

x              = reshape(x,m,n)';
if ~isempty(f), f = reshape(f,m,n)'; end

return