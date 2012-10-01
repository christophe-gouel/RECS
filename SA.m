function [x,fval,exitflag] = SA(f,x,options,varargin)
% SA Solves a system of equations by successive approximation
%
% SA solves a system of equations by iterating on the evaluations of
% the equations: x[i+1] = x[i]+lambda*f(x[i]) until x has converged.
%
% X = SA(F,X0) tries to solve the system of equations F(X)=0 and
% start at the vector X0. F accepts a vector X and return a vector
% of equation values F evaluated at X. SA returns the vector X the
% root of F.
%
% X = SA(F,X0,OPTIONS) solves the problem using the options defined
% in the structure OPTIONS. Fields can be
%      maxit     : maximum number of iterations (default: 1000)
%      atol      : absolute convergence tolerance (default: sqrt(eps))
%      rtol      : relative convergence tolerance (default: sqrt(eps))
%      showiters : 1 to display results of each iteration, 0 (default) if not
%      lambda    : adjustment parameter (default: 1)
%
% X = SA(F,X0,OPTIONS,VARARGIN) provides additional arguments for
% F, which, in this case, takes the following form: F(X,VARARGIN).
%
% [X,FVAL] = SA(F,X0,...) returns the value of the equations F at X.
%
% [X,FVAL,EXITFLAG] = SA(F,X0,...) returns EXITFLAG that describes
% the exit conditions. Possible values listed below.
%       1 : SA converged to a root
%       0 : Too many iterations

% Copyright (C) 2011-2012 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
defaultopt = struct(...
    'maxit',1000,...
    'atol',sqrt(eps),...
    'rtol',sqrt(eps),...
    'showiters',0,...
    'lambda',1);
if nargin < 3
  options = defaultopt; 
else
  warning('off','catstruct:DuplicatesFound')
  options = catstruct(defaultopt,options);
end

maxit     = options.maxit;
atol      = options.atol;
rtol      = options.rtol;
showiters = options.showiters;
lambda    = options.lambda;

fval     = feval(f,x,varargin{:});
fnrm     = norm(fval);
stop_tol = atol + rtol*fnrm;
it       = 0;

%% Iterations
while(fnrm > stop_tol && it < maxit)
  x_old = x;
  it    = it+1;
  if it==1
    x   = x+lambda*fval;
  else
    x   = x+lambda*feval(f,x,varargin{:});
  end
  fnrm  = norm(x-x_old);
  if showiters
    if it==1
      fprintf(1,'Successive approximation\n');
      fprintf(1,'  Iteration\tResidual\n');
    end
    fprintf(1,'%8i\t%5.1E\n',it,fnrm);
  end
end

if it==maxit
  exitflag = 0;
else
  exitflag = 1;
end
