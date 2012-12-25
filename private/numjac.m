function J = numjac(FUN,x,epsilon,SerialProblem,varargin)
% NUMJAC Computes two-sided finite difference Jacobian
%
% J = NUMJAC(FUN,X) computes the Jacobian J by two-sided finite difference
% applied to function FUN (function name or function handle) that accepts a vector
% X as input and return a vector F. The Jacobian is evaluated at X.
%
% J = NUMJAC(FUN,X,epsilon) uses epsilon as step size (default value=srqt(eps))
%
% J = NUMJAC(FUN,X,epsilon,SerialProblem) considers the problem as a serial
% problem: X is a n-by-d matrix and FUN(X) outputs a n-by-d matrix, such that
% each line is independent from the other. In this case, NUMJAC returns a sparse
% block diagonal Jacobian.

% Copyright (C) 2011-2012 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
if isa(FUN,'char')
  FUN = str2func(FUN);
elseif ~isa(FUN,'function_handle')
  error('FUN must be either a string or a function handle')
end

if nargin<=2 || isempty(epsilon), epsilon = sqrt(eps); end

if nargin<4, SerialProblem = 0; end

scale = max(abs(x),1);
step  = epsilon.*scale;
x1    = x+step;
x0    = x-step;
  
%% Computation of finite difference Jacobian
if SerialProblem
  %% Serial problem
  [n,d]      = size(x);
  J          = zeros(n,d,d);
  for j=1:d
    xx       = x;
    xx(:,j)  = x1(:,j);
    f1       = FUN(xx,varargin{:});
    xx(:,j)  = x0(:,j);
    f0       = FUN(xx,varargin{:});
    J(:,:,j) = (f1-f0)./(2*step(:,j));
  end
  J = spblkdiag(permute(J,[2 3 1]));
else
  %% Traditional problem
  d     = length(x);
  for j=1:d
    xx     = x;
    xx(j)  = x1(j);
    f1     = FUN(xx,varargin{:});
    xx(j)  = x0(j);
    f0     = FUN(xx,varargin{:});
    if j==1, J = zeros(length(f0),d); end
    J(:,j) = (f1-f0)/(2*step(j));
  end
end

