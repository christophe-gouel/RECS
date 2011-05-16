% STO4 One small-country storage-trade model

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

disp('STO4 One small-country storage-trade model');
clear interp model options

% ENTER MODEL PARAMETERS
delta = 0;
r     = 0.05;
k     = 0.02;
alpha = -0.4;
tau   = 0.2;
rho   = 0.6;
sigma = 0.16;

% COMPUTE SHOCK DISTRIBUTION
Mu                = [1 0];
Sigma             = [0.05 0;0 sigma];
[model.e,model.w] = qnwnorm([5 5],Mu,Sigma^2);
model.funrand     = @(nrep) Mu(ones(nrep,1),:)+randn(nrep,2)*Sigma;

% PACK MODEL STRUCTURE
model.func   = @sto4model;                                     % model functions
model.params = {delta,r,k,alpha,tau,rho,sigma};               % other parameters

% DEFINE APPROXIMATION SPACE
order         = [15; 15];                                      % degree of approximation
smin          = [min(model.e(:,1))*0.95; 0.4];
smax          = [1.6; 2];
interp.fspace = fundefn('spli',order,smin,smax);               % function space
snodes        = funnode(interp.fspace);                        % state collocaton nodes
s             = gridmake(snodes);
interp.Phi    = funbasx(interp.fspace);
n             = prod(order);

xinit         = [zeros(n,1) max(min(s(:,1).^(1/alpha),s(:,2)+tau),s(:,2)-tau) zeros(n,2)];
interp.cz     = max(min(ones(n,1),s(:,2).^rho+tau),s(:,2).^rho-tau);
interp.cx     = xinit;
interp.ch     = max(min(s(:,1).^(1/alpha),s(:,2)+tau),s(:,2)-tau);

options = struct('eqsolver','lmmcp',...
                 'reesolver','SA',...
                 'simulmethod','solve',...
                 'method','expfunapprox',...
                 'stat',1);

tic
[interp.ch,x,z] = recsSolveREE(interp,model,s,xinit,options);
toc
interp.cx = funfitxy(interp.fspace,interp.Phi,x);
interp.cz = funfitxy(interp.fspace,interp.Phi,z);

if exist('KINSol')
  options.reesolver = 'kinsol';
  interp.ch     = max(min(s(:,1).^(1/alpha),s(:,2)+tau),s(:,2)-tau);
  tic
    [interp.ch,x,z] = recsSolveREE(interp,model,s,xinit,options);
  toc
  interp.cx = funfitxy(interp.fspace,interp.Phi,x);
  interp.cz = funfitxy(interp.fspace,interp.Phi,z);
end

reset(RandStream.getDefaultStream);
[ssim,xsim,esim] = recsSimul(model,interp,[1 1],200,[],options);

