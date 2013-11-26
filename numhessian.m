function H = numhessian(FUN,x,options,varargin)
% NUMHESSIAN Computes a 4-point central finite differences Hessian
%
% H = NUMHESSIAN(FUN,X) computes the Hessian H by 4-point central finite
% difference applied to function FUN that accepts a vector X as input and return
% a scalar. The Hessian is evaluated at X. The elements of the Hessian are
% evaluated using the following formula:
%
% d2f/dxdy = [f(.,x+h,y+k,.)+f(.,x-h,y-k,.)-f(.,x+h,y-k,.)-f(.,x-h,y+k,.)]/4*h*k
%
% The order of approximation is O(hk).
%
% H = NUMHESSIAN(FUN,X,OPTIONS) uses the structure OPTIONS to calculate the
% Hessian. The fields of the structure are
%  DiagOnly       : true or false (default), if true calculates only the diagonal
%                   elements of the Hessian.
%  FinDiffRelStep : scalar or vector step size factor, the default value is
%                   FinDiffRelStep = 1E-4.
%  UseParallel    : 'always' to use parallel calculation (require Parallel
%                    Computing Toolbox)' or never' (default)

% Copyright (C) 2011-2013 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
% FUN
if isa(FUN,'char')
  FUN = str2func(FUN);
elseif ~isa(FUN,'function_handle')
  error('FUN must be either a string or a function handle')
end

% x
x  = x(:);
n  = length(x);

% options
defaultopt = struct(             ...
    'DiagOnly'       , false    ,...
    'FinDiffRelStep' , 1E-4     ,...
    'UseParallel'    , 'never');
if nargin <= 2
  options = defaultopt;
else
  warning('off','catstruct:DuplicatesFound')
  options = catstruct(defaultopt,options);
end
DiagOnly = options.DiagOnly;
switch lower(options.UseParallel)
  case 'never'
    UseParallel = 0;
  case 'always'
    UseParallel = n;
end

%% Compute the stepsize
step = options.FinDiffRelStep.*max(abs(x),1);
ee   = sparse(1:n,1:n,step,n,n);

%% Calculation of the hessian
f = FUN(x,varargin{:});

if ~DiagOnly
  %% All hessian

  H = 4*(step*step');
  parfor (i=1:n, UseParallel)
    for j=1:n
      fpp    = FUN(x+ee(:,i)+ee(:,j),varargin{:});
      fmm    = FUN(x-ee(:,i)-ee(:,j),varargin{:});
      if i~=j
        fpm    = FUN(x+ee(:,i)-ee(:,j),varargin{:});
        fmp    = FUN(x-ee(:,i)+ee(:,j),varargin{:});
        H(i,j) = (fmm+fpp-fpm-fmp)/H(i,j);
      else
        H(i,j) = (fmm+fpp-2*f)/H(i,j);
      end
    end
  end
  % Impose symmetry (an alternative and faster approach would be to do above a
  % loop for j=i:n)
  H = (H+H')/2;

else
  %% Diagonal only

  H = 4*step.^2;
  parfor (i=1:n, UseParallel)
    fpp    = FUN(x+2*ee(:,i),varargin{:}); %#ok<*PFBNS>
    fmm    = FUN(x-2*ee(:,i),varargin{:});
    H(i)   = (fmm+fpp-2*f)/H(i);
  end

end
