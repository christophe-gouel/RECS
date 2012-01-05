function [z,f,exitflag,J,mu] = pathmcp(z,l,u,cpfj,A,b,t,mu)
% PATHMCP solves a polyhedrally constrained variational inequality using PATH
%
% Z = PATHMCP(Z,L,U,CPFJ) tries to solve, using Z as a starting point, the mixed
% complementarity problem of the form:
% L =Z    =>   CPFJ(Z)>0,
% L<=Z<=U =>   CPFJ(Z)=0,
%    Z =U =>   CPFJ(Z)<0.
% L and U are the lower and upper bounds on Z. PATHMCP returns Z the solution.
% CPFJ is the name of the m-file for evaluating the function F and its Jacobian
% J (without .m-extension). The m-file must be supplied (where default name is
% 'mcp_funjac.m' unless stated otherwise in the variable CPFJ). 'mcp_funjac.m'
% contains function [F,J,DOMERR]=CPFJ(Z,JACFLAG) that computes the function F
% and if JACFLAG=1 the sparse Jacobian J at the point Z. DOMERR returns the
% number of domain violations.
% Solver options can be defined through an option file present in the working
% directory and named 'path.opt'. Many options are described in the following file:
% http://www.cs.wisc.edu/~ferris/path/options.pdf
%
% Z = PATHMCP(Z,L,U,CPFJ,A,B,T,MU)
%  A  - constraint matrix
%  B  - right hand side of the constraints
%  T  - types of the constraints
%      <0 : less than or equal
%      =0 : equal to
%      >0 : greater than or equal
%      We have AZ ? B, ? is the type of constraint
%  MU - multipliers on the constraints
%
% [Z,F] = PATHMCP(Z,L,U,CPFJ,...) returns F the function evaluation at the
% solution.
%
% [Z,F,EXITFLAG] = PATHMCP(Z,L,U,CPFJ,...) returns EXITFLAG that describes
% the exit conditions. Possible values listed below.
%       1 : solved
%       0 : failed to solve
%
% [Z,F,EXITFLAG,J] = PATHMCP(Z,L,U,CPFJ,...) returns J the Jacobian evaluation
% at the solution.
%
% [Z,F,EXITFLAG,J,MU] = PATHMCP(Z,L,U,CPFJ,A,B,T,MU) returns MU the multipliers
% on the constraints at the solution.
%
% For more information, see the following references 
% Dirkse, S. P. and Ferris, M. C. (1995), The PATH solver: A non-monotone
%   stabilization scheme for mixed complementarity problems, Optimization Methods
%   and Software 5, 123-156. DOI: <a href="http://dx.doi.org/10.1080/10556789508805606">10.1080/10556789508805606</a>
% Ferris, M. C. and Munson, T. S. (1999), Interfaces to PATH 3.0: Design,
%   Implementation and Usage, Computational Optimization and Applications 12,
%   207-227. DOI: <a href="http://dx.doi.org/10.1023/A:1008636318275">10.1023/A:1008636318275</a>
% Munson, T. S. (2002), Algorithms and Environments for Complementarity, PhD
%   thesis, University of Wisonsin-Madison.

% Copyright (C) Michael C. Ferris and Todd S. Munson

%% Initialization
Big = 1e20;

if (nargin < 1)
  error('one input arguments required for mcp(z)');
end

z = full(z(:));
n = length(z);

if (n == 0)
  error('empty model');
end

if (nargin < 2 || isempty(l))
  l = zeros(n,1);
end

if (nargin < 3 || isempty(u))
  u = Big*ones(n,1);
end

l = full(l(:)); 
u = full(u(:));
if (length(l) ~= n || length(u) ~= n)
  error('Input arguments are of incompatible sizes');
end

if any(u-l<0)
  error('Some upper bounds are inferior to their corresponding lower bounds')
end

l = max(l,-Big*ones(n,1));
u = min(u,Big*ones(n,1));
z = min(max(z,l),u);

if (nargin < 4 || isempty(cpfj))
  cpfj = 'mcp_funjac';
end

m   = 0;
mu  = [];
l_p = [];
u_p = [];

%%
if (nargin > 4)
  %% With Polyhedral constraints
  if (nargin < 6)
    error('Polyhedral constraints require A and b');
  end

  if (~issparse(A))
    A = sparse(A);
  end
  b = full(b(:));

  m = length(b);

  if (m > 0)

    [am, an] = size(A);

    if (am ~= m || an ~= n)
      error('Polyhedral constraints of incompatible sizes');
    end

    if (nargin < 7 || isempty(t))
      t = ones(m,1);
    end

    if (nargin < 8 || isempty(mu))
      mu = zeros(m,1);
    end

    t = full(t(:)); mu = full(mu(:));
    if (length(t) ~= m || length(mu) ~= m)
      error('Polyhedral input arguments are of incompatible sizes');
    end

    l_p = -Big*ones(m,1);
    u_p =  Big*ones(m,1);

    idx = find(t > 0);
    if (length(idx) > 0)
      l_p(idx) = zeros(length(idx),1);
    end

    idx = find(t < 0);
    if (length(idx) > 0)
      u_p(idx) = zeros(length(idx),1);
    end

    mu = min(max(mu,l_p),u_p);
  else
    if (nargin >= 8 && ~isempty(mu))
      error('No polyhedral constraints -- multipliers set.');
    end

    if (nargin >= 7 && ~isempty(t))
      error('No polyhedral constraints -- equation types set.');
    end
  end
else
  A = [];
end

%% Check for domain error at starting point and Jacobian
% this is a fix, nnz may be bigger than this
[~,J,domerr] = feval(cpfj,z+1e-5*ones(size(z))+1e-5*abs(z),1);

if (domerr > 0)
  [~,J,domerr] = feval(cpfj,z,1);
end

if (domerr > 0)
  error([cpfj ' not defined at starting point']);
end

if ~issparse(J)
  error([cpfj ' must return a sparse Jacobian']);
end

nnzJ = nzmax(J);

row = n + m;
% Change CG
%ele = nnzJ + 2*nzmax(A);
ele = nnzJ*1.05 + 2*nzmax(A);

%% Solve the problem
init = [z; mu];
low  = [l; l_p];
upp  = [u; u_p];

if m > 0
  global mcp_vifunc;
  global mcp_viconn;
  global mcp_viconm;
  global mcp_viconA;
  global mcp_viconb;

  mcp_vifunc = cpfj;
  mcp_viconn = n;
  mcp_viconm = m;
  mcp_viconA = A;
  mcp_viconb = b;

  [exitflag, ~, f, J] = mcppath(row, ele, init, low, upp, 'mcp_vifunjac');
else
  [exitflag, ~, f, J] = mcppath(row, ele, init, low, upp, cpfj);
end

mu = [];
z  = init;

%%
if m > 0
  mu = init(n+1:n+m);
  z  = init(1:n);

  J  = J(1:n,1:n);
  f  = f(1:n) + A'*mu;
end
