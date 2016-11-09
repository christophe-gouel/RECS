function [x,fval,exitflag] = recsSolveEquilibrium(s,x,z,b,f,g,h,params,c,e,w,fspace,ixforward,options,LB,UB)
% RECSSOLVEEQUILIBRIUM Solves the system of equilibrium equations using x as starting values
%
% RECSSOLVEEQUILIBRIUM is called by recsSimul, recsSolveREEFiniteHorizon,
% and recsSolveREEIter. It is not meant to be called directly by the user.
%
% See also RECSSIMUL, RECSSOLVEREE.

% Copyright (C) 2011-2016 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
[n,m]            = size(x);
if nargin<=14, [LB,UB] = b(s,params); end

ArrayProblem     = options.ArrayProblem;
eqsolver         = lower(options.eqsolver);
eqsolveroptions  = options.eqsolveroptions;
eqsolveroptions.ArrayProblem = ArrayProblem;
extrapolate      = options.extrapolate;
funapprox        = lower(options.funapprox);
loop_over_s      = options.loop_over_s;
switch lower(options.UseParallel)
  case 'never'
    UseParallel = 0;
  case 'always'
    UseParallel = n;
end

%% Solve equilibrium equations on grid points
if loop_over_s
  if loop_over_s==1
    %% Solve separately for each point on the grid
    fval         = zeros(size(x));
    [~,grid]     = spblkdiag(zeros(m,m,1),[],0);
    exitflag     = zeros(n,1);
    parfor (i=1:n, UseParallel)
      [x(i,:),fval(i,:),exitflag(i)] = eqsolve(x(i,:),s(i,:),z(i,:),...
                                               b,f,g,h,params,eqsolver,...
                                               grid,c,e,w,fspace,funapprox,...
                                               eqsolveroptions,...
                                               LB(i,:),UB(i,:),extrapolate,...
                                               ixforward);
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
    fval         = cell(loop_over_s,1);
    parfor (i=1:loop_over_s, UseParallel)
      [~,grid]   = spblkdiag(zeros(m,m,sizeBlocks(i)),[],0);
      [x{i},fval{i},exitflag(i)] = eqsolve(x{i},s{i},z{i},...
                                           b,f,g,h,params,eqsolver,...
                                           grid,c,e,w,fspace,funapprox,...
                                           eqsolveroptions,...
                                           LB{i},UB{i},extrapolate,...
                                           ixforward);
    end
    x            = cell2mat(x);
    fval         = cell2mat(fval);
  end % if loop_over_s==1
  exitflag       = all(exitflag);
else
  %% Solve all the grid in one step
  if ~ArrayProblem
    [~,grid] = spblkdiag(zeros(m,m,n),[],0); 
  else
    grid     = [];
  end
  [x,fval,exitflag] = eqsolve(x,s,z,b,f,g,h,params,eqsolver,grid,c,e,w,fspace,...
                              funapprox,eqsolveroptions,LB,UB,extrapolate,...
                              ixforward,ArrayProblem);
end % if loop_over_s

function [x,fval,exitflag] = eqsolve(x,s,z,b,f,g,h,params,eqsolver,grid,c,e,w,...
                                     fspace,funapprox,eqsolveroptions,LB,UB,...
                                     extrapolate,ixforward,ArrayProblem)
%% Eqsolve

if ~ArrayProblem
  [n,m]          = size(x);
  x              = reshape(x',[n*m 1]);
  LB             = reshape(LB',[n*m 1]);
  UB             = reshape(UB',[n*m 1]);
end

[x,fval,exitflag] = runeqsolver(@recsEquilibrium,x,LB,UB,eqsolver,eqsolveroptions,...
                                s,z,b,f,g,h,params,grid,c,e,w,fspace,funapprox,...
                                extrapolate,ixforward,ArrayProblem);

if exitflag~=1, disp('No convergence'); end

if ~ArrayProblem
  x              = reshape(x,m,n)';
  if ~isempty(fval), fval = reshape(fval,m,n)'; end
end

return