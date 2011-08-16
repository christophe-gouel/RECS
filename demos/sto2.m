% STO2 Competitive storage with floor-price backed by public storage

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

disp('STO2 Competitive storage with floor-price backed by public storage');
clear interp model options

% ENTER MODEL PARAMETERS
k     = 0.02;
delta = 0.02;
r     = 0.05;
mu    = 10;
alpha = -0.4;
PF    = 1.02;
Sgbar = 0.4;

% COMPUTE SHOCK DISTRIBUTION
Mu                = 1;
sigma             = 0.05;
[model.e,model.w] = qnwnorm(5,Mu,sigma^2);
model.funrand     = @(nrep) Mu+sigma*randn(nrep,1);

% PACK MODEL STRUCTURE
model.func   = @sto2model;                               % model functions
model.params = {k,delta,r,mu,alpha,PF,Sgbar};               % other parameters

% DEFINE APPROXIMATION SPACE
order         = 50;                                          % degree of approximation
smin          = min(model.e)*0.95;
smax          = 2;
interp.fspace = fundefn('spli',order,smin,smax);                 % function space
snodes        = funnode(interp.fspace);                             % state collocaton nodes
s             = gridmake(snodes);
interp.Phi    = funbasx(interp.fspace);

optset('ncpsolve','type','minmax'); % 'minmax' / 'smooth'
options.simulmethod = 'solve';

[interp,xinit,zinit] = recsFirstGuess(interp,model,s,1,[0 1 1 0],5);

tic
recsSolveREE(interp,model,s,xinit);
toc

options.method = 'resapprox-simple';
tic
recsSolveREE(interp,model,s,xinit,options);
toc

options.method = 'expapprox';
tic
interp = recsSolveREE(interp,model,s,xinit,options);
toc

options.eqsolver = 'ncpsolve';
interp.cx = funfitxy(interp.fspace,interp.Phi,xinit);
interp.cz = funfitxy(interp.fspace,interp.Phi,zinit);
tic
interp = recsSolveREE(interp,model,s,xinit,options);
toc

tic
[ssim,xsim,esim] = recsSimul(model,interp,1,200,[],options);
toc
