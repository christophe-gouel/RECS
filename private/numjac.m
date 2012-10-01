function J = numjac(FUN,x,epsilon,varargin)
% NUMJAC Computes two-sided finite difference Jacobian
%
% J = NUMJAC(FUN,x) computes the Jacobian J by two-sided finite difference
% applied to function FUN (function name or function handle) that accepts a vector
% x as input and return a vector F. The Jacobian is evaluated at x.
%
% J = NUMJAC(FUN,x,epsilon) uses epsilon as step size (default value=srqt(eps))

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
if isa(FUN,'char')
  FUN = str2func(FUN);
elseif ~isa(FUN,'function_handle')
  error('FUN must be either a string or a function handle')
end

if nargin<=2 || isempty(epsilon), epsilon = sqrt(eps); end

scale = max(abs(x),1);
step  = epsilon.*scale;
x1    = x+step;
x0    = x-step;
n     = length(x);

%% Computation of finite difference Jacobian
for j=1:n
   xx     = x;
   xx(j)  = x1(j);
   f1     = FUN(xx,varargin{:});
   xx(j)  = x0(j);
   f0     = FUN(xx,varargin{:});
   if j==1, J = zeros(length(f0),n); end
   J(:,j) = (f1-f0)/(2*step(j));
end