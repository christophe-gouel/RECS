% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

clear interp model options

% ENTER MODEL PARAMETERS
alpha = -0.4;
k     = 0.02;
delta = 0.02;
r     = 0.05;
mu    = 10;

% COMPUTE SHOCK DISTRIBUTION
Mu                = 1;
sigma             = 0.1;
[model.e,model.w] = qnwnorm(5,Mu,sigma^2);
model.funrand     = @(nrep) normrnd(Mu(ones(nrep,1),:),sigma(ones(nrep,1),:));

% PACK MODEL STRUCTURE
model.func   = 'mtest';                               % model functions
model.params = {alpha,k,delta,r,mu};               % other parameters

% DEFINE APPROXIMATION SPACE
order         = 30;                                          % degree of approximation
smin          = 0.5;
smax          = 2;
interp.fspace = fundefn('spli',order,smin,smax);                 % function space
snodes        = funnode(interp.fspace);                             % state collocaton nodes
s             = gridmake(snodes);
interp.Phi    = funbasx(interp.fspace);
n             = order;

xinit  = [zeros(n,1) ones(n,1) s.^(1/alpha)];
interp.cz = ones(n,2);
interp.cx = xinit;
interp.ch = [s.^(1/alpha) s.^(1/alpha)];
options.eqsolver = 'path';

[sss,xss,zss] = SteadyState(model,1,[0 1 1],options)

s0 = 2;
p  = 2;
T  = 10;
e = model.w'*model.e;
X = [xss zss sss];
X = X(ones(T,1),:)';
X = reshape(X,[],1);

[F,J] = DeterministicProblem(X,'mtest',s0,xss,p,e,model.params);

[x,s,z,f] = SolveDeterministicProblem(model,s0,T,xss,zss,sss,struct('eqsolver','path'));
%[s x z]

A = linspace(0.7,2,100)';
S = zeros(size(A));
H = zeros(size(A));

tic
  for i=1:length(A)
  x = SolveDeterministicProblem(model,A(i),T,xss,zss,sss,struct('eqsolver','path'));
  S(i) = x(1,1);
  H(i) = x(1,2);
  end
toc

plot(A,H);