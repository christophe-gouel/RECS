% GRO1 Stochastic growth model

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

disp('GRO1 Stochastic growth model');
clear interp model options

% ENTER MODEL PARAMETERS
tau   = 2;
beta  = 0.9896;
beta  = 0.95;
alpha = 0.4;
%delta = 0.0196;
delta = 0.1;
rho   = 0.95;
a     = (1/beta-1+delta)/alpha;

% COMPUTE SHOCK DISTRIBUTION
Mu                = 0;
sigma             = 0.007;
[model.e,model.w] = qnwnorm(5,Mu,sigma^2);
model.funrand     = @(nrep) Mu+sigma*randn(nrep,1);

% PACK MODEL STRUCTURE
model.func   = @gro1model;                               % model functions
model.params = {tau,beta,alpha,delta,rho,a};             % other parameters

options = struct(...
    'eqsolver','path',...
    'method','resapprox-complete',...
    'reesolver','SA',...
    'reesolveroptions',struct('lambda',0.5));

disp('Deterministic steady-state')
[out1,out2,out3] = model.func('ss',[],[],[],[],[],[],model.params);
[sss,xss,zss] = recsSS(model,out1,out2,options)

% Check derivatives
recsCheck(model,sss,xss,zss);
return
% DEFINE APPROXIMATION SPACE
order         = [10 10];                                 % degree of approximation
smin          = [0.5*sss(1) min(model.e)*0.95];
smax          = [2*sss(1)   max(model.e)*10.5];
smin          = [0.95*sss(1) -0.009];
smax          = [1.05*sss(1) 0.009];
interp.fspace = fundefn('spli',order,smin,smax);                 % function space
snodes        = funnode(interp.fspace);                             % state collocaton nodes
s             = gridmake(snodes);
interp.Phi    = funbasx(interp.fspace);
n             = prod(order);

[x,s,z,F] = recsSolveDeterministicPb(model,[sss(1)*1.01 sss(2)],10,xss,zss,sss, ...
                                     options)
return
[interp,x,z] = recsFirstGuess(interp,model,s,sss,xss,30,options)
return

xinit  = xss*ones(n,1);
interp.cz = funfitxy(interp.fspace,interp.Phi,zss*ones(n,1));
interp.cx = funfitxy(interp.fspace,interp.Phi,xinit);

tic
[interp.cz,x] = recsSolveREEFull(interp,model,s,xinit,options);
toc

tic
[interp.cz,x] = recsSolveREE(interp,model,s,xinit,options);
toc
interp.cx = funfitxy(interp.fspace,interp.Phi,x);

