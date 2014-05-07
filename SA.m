function [x,fval,exitflag] = SA(f,x,options,varargin)
% SA Solves a system of equations by successive approximation
%
% SA solves a system of equations by iterating on the evaluations of the
% equations: x[i] = x[i-1]+dx[i-1] until norm(f(x[i]))<=TolFun+RelTolFun*norm(f(x[0])),
% where dx = lambda*f(x).
%
% To enhance convergence SA uses a backtracking line search that divides dx by 2
% until norm(f(x+dx)) decreases with respect to norm(f(x)) or until MaxSteps
% iterations have been done.
%
% X = SA(F,X0) tries to solve the system of equations F(X)=0 and
% start at the vector X0. F accepts a vector X and return a vector
% of equation values F evaluated at X. SA returns the vector X the
% root of F.
%
% X = SA(F,X0,OPTIONS) solves the problem using the options defined
% in the structure OPTIONS. Fields can be
%      Display   : 1 to display results of each iteration, 0 (default) if not
%      lambda    : adjustment parameter (default: 1)
%      MaxIter   : maximum number of iterations (default: 1000)
%      MaxSteps  : maximum number of backstepping iterations (default: 3)
%      RelTolFun : relative convergence tolerance (default: sqrt(eps))
%      TolFun    : absolute convergence tolerance (default: sqrt(eps))
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

% Copyright (C) 2011-2014 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
defaultopt = struct(      ...
    'Display'  , 0        ,...
    'lambda'   , 1        ,...
    'MaxIter'  , 1000     ,...
    'MaxSteps' , 3        ,...
    'RelTolFun', sqrt(eps),...
    'TolFun'   , sqrt(eps));
if nargin < 3 || isempty(options)
  options = defaultopt;
else
  warning('off','catstruct:DuplicatesFound')
  options = catstruct(defaultopt,options);
end

lambda    = options.lambda;
MaxIter   = options.MaxIter;
MaxSteps  = max(1,options.MaxSteps);
RelTolFun = options.RelTolFun;
Display   = options.Display;
TolFun    = options.TolFun;

fval      = feval(f,x,varargin{:});
isdomerr  = any(reshape(isinf(fval) | (imag(fval)~=0) | isnan(fval),[],1));
if isdomerr
  exitflag = 0;
  if Display
    fprintf(1,'Equations not defined at starting point\n');
  end
  return
end
fnrm        = norm(fval(:));
stop_tol    = TolFun + RelTolFun*fnrm;
it          = 0;
if Display
  fprintf(1,'Successive approximation\n');
  fprintf(1,'  Major\t Minor\tLipschitz\t Residual\n');
  fprintf(1,'%7i\t%6i\t         \t%9.2E (Input point)\n',0,0,fnrm);
end

%% Iterations
while(fnrm > stop_tol && it < MaxIter)
  it      = it+1;
  dx      = lambda*fval;
  %% Backstepping
  for itback=1:MaxSteps
    fvalnew     = feval(f,x+dx,varargin{:});
    fnrmnew     = norm(fvalnew(:));
    isnotdomerr = ~any(reshape(isinf(fvalnew) | (imag(fvalnew)~=0) | isnan(fvalnew),...
                               [],1));
    if (fnrmnew<fnrm && isnotdomerr) || itback==MaxSteps
      %% End of backstepping because of success or maximum iteration
      k     = norm(fvalnew-fval)/norm(dx);
      fval  = fvalnew;
      fnrm  = fnrmnew;
      x     = x+dx;
      if Display, fprintf(1,'%7i\t%6i\t%8.4f\t%9.2E\n',it,itback,k,fnrm); end
      break
    end
    if itback>1 && fnrmold<fnrmnew
      %% End of backstepping because it does not improve residual
      k    = norm(fvalold-fval)/norm(dx*2);
      fval = fvalold;
      fnrm = fnrmold;
      x    = x+dx*2;
      if Display, fprintf(1,'%7i\t%6i\t%9.4f\t%8.2E\n',it,itback-1,k,fnrm); end
      break
    end
    fvalold = fvalnew;
    fnrmold = fnrmnew;
    dx      = dx/2;
  end % itback
end %it

%% Output treatment
if it==MaxIter
  exitflag = 0;
  if Display
    fprintf(1,'Too many iterations\n');
  end
else
  exitflag = 1;
  if Display
    fprintf(1,'Solution found');
    if fnrm < TolFun
      fprintf(1,' - Residual lower than absolute tolerance\n');
    else
      fprintf(1,' - Residual lower than relative tolerance\n');
    end
  end
end
