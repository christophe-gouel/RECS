% STO2 Competitive storage with floor-price backed by public storage

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

disp('STO2 Competitive storage with floor-price backed by public storage');
clear interp model options

%% Enter model parameters
k     = 0.02;
delta = 0.02;
r     = 0.05;
mu    = 10;
alpha = -0.4;
PF    = 1.02;
Sgbar = 0.4;

%% Compute shock distribution
Mu                = 1;
sigma             = 0.05;
[model.e,model.w] = qnwnorm(5,Mu,sigma^2);
model.funrand     = @(nrep) Mu+sigma*randn(nrep,1);

%% Pack model structure
model.func   = @sto2model;                               % model functions
model.params = {k,delta,r,mu,alpha,PF,Sgbar};               % other parameters

%% Define approximation space
order         = 50;                                          % degree of approximation
smin          = min(model.e)*0.95;
smax          = 2;
interp.fspace = fundefn('spli',order,smin,smax);                 % function space
snodes        = funnode(interp.fspace);                             % state collocaton nodes
s             = gridmake(snodes);
interp.Phi    = funbasx(interp.fspace);

%% Find a first guess through the perfect foresight solution
[interp,xinit,zinit] = recsFirstGuess(interp,model,s,1,[0 1 1 0],5);

%% Solve for rational expectations
tic
recsSolveREE(interp,model,s,xinit);
toc

options.funapprox = 'resapprox-simple';
tic
recsSolveREE(interp,model,s,xinit,options);
toc

options.funapprox = 'expapprox';
tic
interp = recsSolveREE(interp,model,s,xinit,options);
toc

options.eqsolver = 'ncpsolve';
optset('ncpsolve','type','minmax'); % 'minmax' / 'smooth'
interp.cx = funfitxy(interp.fspace,interp.Phi,xinit);
interp.cz = funfitxy(interp.fspace,interp.Phi,zinit);
tic
interp = recsSolveREE(interp,model,s,xinit,options);
toc

%% Simulate the model
options.simulmethod = 'solve';
tic
[ssim,xsim,esim] = recsSimul(model,interp,1,200,[],options);
toc
