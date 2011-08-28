% STO1 Competitive storage with supply reaction and explicit market equilibrium equation

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

disp('STO1 Competitive storage with supply reaction and explicit market equilibrium equation');
clear interp model options

%% Enter model parameters
alpha = -0.4;
k     = 0.02;
delta = 0.02;
r     = 0.05;
mu    = 10;

%% Compute shock distribution
Mu                = 1;
sigma             = 0.1;
[model.e,model.w] = qnwnorm(5,Mu,sigma^2);
model.funrand     = @(nrep) Mu+sigma*randn(nrep,1);

%% Pack model structure
model.func   = @sto1model;                               % model functions
model.params = {alpha,k,delta,r,mu};               % other parameters

%% Define approximation space
order         = 30;                                          % degree of approximation
smin          = 0.5;
smax          = 2;
interp.fspace = fundefn('spli',order,smin,smax);                 % function space
snodes        = funnode(interp.fspace);                             % state collocaton nodes
s             = gridmake(snodes);
interp.Phi    = funbasx(interp.fspace);
n             = order;

%% Provide a first guess
xinit  = [zeros(n,1) ones(n,1) s.^(1/alpha)];
interp.cz = funfitxy(interp.fspace,interp.Phi,ones(n,2));
interp.cx = funfitxy(interp.fspace,interp.Phi,xinit);
interp.ch = funfitxy(interp.fspace,interp.Phi,[s.^(1/alpha) s.^(1/alpha)]);

%% Find deterministic steady-state
disp('Deterministic steady-state')
[sss,xss,zss] = recsSS(model,1,[0 1 1])

%% Check derivatives
recsCheck(model,sss,xss,zss);

%% Solve for rational expectations
tic
recsSolveREE(interp,model,s,xinit);
toc

options = struct(...
    'method','expapprox',...
    'reesolver','krylov');
tic
recsSolveREE(interp,model,s,xinit,options);
toc

options.method = 'expfunapprox';
tic
recsSolveREE(interp,model,s,xinit,options);
toc

options.method = 'resapprox-simple';
tic
interp = recsSolveREE(interp,model,s,xinit,options);
toc

if exist('fsolve','file')
  options.reesolver = 'fsolve';
  interp.cz     = ones(n,2);
  interp.cx     = xinit;
  tic
    interp = recsSolveREE(interp,model,s,xinit,options);
  toc
end

if exist('kinsol','file')
  options.reesolver = 'kinsol';
  interp.cz     = ones(n,2);
  interp.cx     = xinit;
  tic
    interp = recsSolveREE(interp,model,s,xinit,options);
  toc
end

%% Simulate the model
reset(RandStream.getDefaultStream);
[ssim,xsim,esim,~,stat] = recsSimul(model,interp,ones(1000,1),200);

% Solution accuracy
recsAccuracy(model,interp,ssim)

figure
subplot(1,3,1)
plot(ssim(:),reshape(xsim(:,1,:),[],1),'.')
title('Storage')
subplot(1,3,2)
plot(ssim(:),reshape(xsim(:,2,:),[],1),'.')
title('Planned production')
subplot(1,3,3)
plot(ssim(:),reshape(xsim(:,3,:),[],1),'.')
title('Price')

% Simulation using the shocks discretisation and not the random number generator
% function provided by the user
model = rmfield(model,'funrand');
[ssim,xsim,esim] = recsSimul(model,interp,ones(1000,1),200);
