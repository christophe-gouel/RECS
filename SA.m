function [x,fval,exitflag] = SA(f,x,options,varargin)
% SA Solves a system of equations by successive approximation
% USAGE
%   [x,fval,exitflag] = SA(f,x,options,varargin)
% INPUTS
%   f        : name of function of form:
%               fval=f(x,optional additional parameters)
%   x        : initial guess for root (d by 1)
%   varargin : additional arguments for f [optional]
% OUTPUTS
%   x        : root of f (d by 1)
%   fval     : function value estimate (d by 1)
%   options  : structure defining options. Fields can be
%               maxit     : maximum number of iterations (default: 1000)
%               atol      : absolute convergence tolerance
%               rtol      : relative convergence tolerance
%               showiters : display results of each iteration
%               lambda    : adjustment parameter
%
%   exitflag : integer describing exit condition. Possible values listed below.
%                1  SA converged to a root
%                0  Too many iterations

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

% ------------Initialization----------------
defaultopt = struct(...
    'maxit',1000,...
    'atol',sqrt(eps),...
    'rtol',sqrt(eps),...
    'showiters',0,...
    'lambda',1);

if nargin < 3, options = []; end

options = catstruct(defaultopt,options);

maxit     = options.maxit;
atol      = options.atol;
rtol      = options.rtol;
showiters = options.showiters;
lambda    = options.lambda;

fval     = feval(f,x,varargin{:});
fnrm     = norm(fval);
stop_tol = atol + rtol*fnrm;
it       = 0;

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
