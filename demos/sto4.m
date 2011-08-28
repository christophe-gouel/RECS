% STO4 One small-country storage-trade model

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

disp('STO4 One small-country storage-trade model');
clear interp model options

%% Enter model parameters
delta = 0;
r     = 0.05;
k     = 0.02;
alpha = -0.4;
tau   = 0.2;
rho   = 0.6;
sigma = 0.16;

%% Compute shock distribution
Mu                = [1 0];
Sigma             = [0.05 0;0 sigma];
[model.e,model.w] = qnwnorm([5 5],Mu,Sigma^2);
model.funrand     = @(nrep) Mu(ones(nrep,1),:)+randn(nrep,2)*Sigma;

%% Pack model structure
model.func   = @sto4model;                                     % model functions
model.params = {delta,r,k,alpha,tau,rho,sigma};               % other parameters

%% Define approximation space
order         = [15; 15];                                      % degree of approximation
smin          = [min(model.e(:,1))*0.95; 0.4];
smax          = [1.6; 2.05];
interp.fspace = fundefn('spli',order,smin,smax);               % function space
snodes        = funnode(interp.fspace);                        % state collocaton nodes
s             = gridmake(snodes);
interp.Phi    = funbasx(interp.fspace);
n             = prod(order);

%% Provide a first guess
xinit         = [zeros(n,1) max(min(s(:,1).^(1/alpha),s(:,2)+tau),s(:,2)-tau) zeros(n,2)];
interp.cz     = max(min(ones(n,1),s(:,2).^rho+tau),s(:,2).^rho-tau);
interp.cx     = xinit;
interp.ch     = max(min(s(:,1).^(1/alpha),s(:,2)+tau),s(:,2)-tau);

%% Solve for rational expectations
options = struct('simulmethod','solve',...
                 'method','expfunapprox',...
                 'stat',1);

% Solve by Full Newton
if exist('mcppath','file')
  options.eqsolver = 'path';
  tic
    recsSolveREEFull(interp,model,s,xinit,options);
  toc
  options.eqsolver = 'lmmcp';
end

tic
interp = recsSolveREE(interp,model,s,xinit,options);
toc

if exist('KINSol','file')
  options.reesolver = 'kinsol';
  interp.ch     = max(min(s(:,1).^(1/alpha),s(:,2)+tau),s(:,2)-tau);
  tic
    interp = recsSolveREE(interp,model,s,xinit,options);
  toc
end

%% Simulate the model
reset(RandStream.getDefaultStream);
[ssim,xsim,esim] = recsSimul(model,interp,[1 1],200,[],options);

