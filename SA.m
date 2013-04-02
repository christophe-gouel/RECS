function [x,fval,exitflag] = SA(f,x,options,varargin)
% SA Solves a system of equations by successive approximation
%
% SA solves a system of equations by iterating on the evaluations of the
% equations: x[i+1] = x[i]+dx[i] until norm(f(x)) has converged to 0, where
% dx = lambda*f(x).
%
% To enhance convergence SA uses a backtracking line search that divides dx by 2
% until norm(f(x+dx)) decreases with respect to norm(f(x)).
%
% X = SA(F,X0) tries to solve the system of equations F(X)=0 and
% start at the vector X0. F accepts a vector X and return a vector
% of equation values F evaluated at X. SA returns the vector X the
% root of F.
%
% X = SA(F,X0,OPTIONS) solves the problem using the options defined
% in the structure OPTIONS. Fields can be
%      maxit     : maximum number of iterations (default: 1000)
%      maxsteps  : maximum number of backstepping iterations (default: 3)
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
%       0 : Failure to converge because of too many iterations or equations not
%           defined at starting point

% Copyright (C) 2011-2013 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
defaultopt = struct(      ...
    'atol'     ,sqrt(eps),...
    'lambda'   ,1        ,...
    'maxit'    ,1000     ,...
    'maxsteps' ,3        ,...
    'rtol'     ,sqrt(eps),...
    'showiters',0);
if nargin < 3 || isempty(options)
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
maxsteps  = max(1,options.maxsteps);

fval      = feval(f,x,varargin{:});
isdomerr  = any(reshape(isinf(fval) | (imag(fval)~=0) | isnan(fval),[],1));
if isdomerr
  exitflag = 0;
  if showiters
    fprintf(1,'Equations not defined at starting point\n');
  end
  return
end
fnrm        = norm(fval(:));
stop_tol    = atol + rtol*fnrm;
it          = 0;
if showiters
  fprintf(1,'Successive approximation\n');
  fprintf(1,'  Major\t Minor\tResidual\n');
  fprintf(1,'%7i\t%6i\t%8.1E (Input point)\n',0,0,fnrm);
end

%% Iterations
while(fnrm > stop_tol && it < maxit)
  it      = it+1;
  dx      = lambda*fval;
  %% Backstepping
  for itback=1:maxsteps 
    fvalnew     = feval(f,x+dx,varargin{:});
    fnrmnew     = norm(fvalnew(:));
    isnotdomerr = ~any(reshape(isinf(fvalnew) | (imag(fvalnew)~=0) | isnan(fvalnew),...
                               [],1));
    if (fnrmnew<fnrm && isnotdomerr) || itback==maxsteps
      %% End of backstepping because of success or maximum iteration
      fval  = fvalnew;
      fnrm  = fnrmnew;
      x     = x+dx;
      if showiters, fprintf(1,'%7i\t%6i\t%8.1E\n',it,itback,fnrm); end
      break
    end
    if itback>1 && fnrmold<fnrmnew
      %% End of backstepping because it does not improve residual
      fval = fvalold;
      fnrm = fnrmold;
      x    = x+dx*2;
      if showiters, fprintf(1,'%7i\t%6i\t%8.1E\n',it,itback-1,fnrm); end
      break
    end
    fvalold = fvalnew;
    fnrmold = fnrmnew;
    dx      = dx/2;
  end % itback
end %it

%% Output treatment
if it==maxit
  exitflag = 0;
  if showiters
    fprintf(1,'Too many iterations\n');
  end
else
  exitflag = 1;
  if showiters
    fprintf(1,'Solution found');
    if fnrm < atol
      fprintf(1,' - Residual lower than absolute tolerance\n');
    else
      fprintf(1,' - Residual lower than relative tolerance\n');
    end
  end
end
