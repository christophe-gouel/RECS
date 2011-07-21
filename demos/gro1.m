% GRO1 Stochastic growth model

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

disp('GRO1 Stochastic growth model');
clear interp model options

% COMPUTE SHOCK DISTRIBUTION
Mu                = 0;
sigma             = 0.007;
[model.e,model.w] = qnwnorm(5,Mu,sigma^2);
model.funrand     = @(nrep) Mu+sigma*randn(nrep,1);

% PACK MODEL STRUCTURE
model.func   = @gro1model;
model.params = gro1model('params');
a            = model.params(1);
delta        = model.params(3);

options = struct('reesolver','krylov');

disp('Deterministic steady-state')
[sss,xss,zss] = recsSS(model,[1 0],a-delta,options)

% DEFINE APPROXIMATION SPACE
order         = [10 10];                                 % degree of approximation
smin          = [0.85*sss(1) min(model.e)*4];
smax          = [1.15*sss(1) max(model.e)*4];
interp.fspace = fundefn('cheb',order,smin,smax);           % function space
snodes        = funnode(interp.fspace);                    % state collocaton nodes
s             = gridmake(snodes);
interp.Phi    = funbasx(interp.fspace);
n             = prod(order);

disp('Find the perfect-foresight solution as first guess')
tic
[interp,x,z] = recsFirstGuess(interp,model,s,sss,xss,50,options);
toc

disp('Solve the stochastic problem:')
tic
[interp.cx,x,z] = recsSolveREE(interp,model,s,x,options);
toc
interp.cz = funfitxy(interp.fspace,interp.Phi,z);

[~,~,~,~,stat] = recsSimul(model,interp,sss(ones(1000,1),:),200);


