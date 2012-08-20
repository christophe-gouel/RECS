%% STO4 Storage-trade model of a small country

%% Model's structure
%
% *Response variables* Storage ($S$), Price ($P$), Import ($M$) and Export
% ($X$).
%
% *State variable* Availability ($A$) and World price ($P^{w}$).
%
% *Shock* Production shocks ($\epsilon$) and Innovation on world price ($\nu$).
%
% *Parameters* Unit physical storage cost ($k$), Depreciation share ($\delta$),
% Interest rate ($r$), Demand price elasticity ($\alpha$), World price
% autocorrelation ($\rho$) and Trade cost ($\theta$).
%
% *Equilibrium equations*
%
% $$S_{t}: S_{t}\ge 0 \quad \perp \quad P_{t}+k-\frac{1-\delta}{1+r}\mathrm{E}_{t}\left(P_{t+1}\right)\ge 0,$$
%
% $$P_{t}: A_{t}+M_t={P_{t}}^{\alpha}+S_{t}+X_t,$$
%
% $$M_{t}: M_{t}\ge 0 \quad \perp \quad P^w_{t}+\theta\ge P_{t},$$
%
% $$X_{t}: X_{t}\ge 0 \quad \perp \quad P_{t}\ge P^w_{t}-\theta.$$
%
% *Transition equation*
%
% $$A_{t}: A_{t}=\left(1-\delta\right)S_{t-1}+\epsilon_{t},$$
%
% $$P^w_{t}: \ln P^w_{t} = \rho\ln P^w_{t-1}+\nu_{t}.$$

%% Writing the model
% The model is defined in a Matlab file: <sto4model.html sto4model.m>.

%% Enter model parameters
delta = 0;
r     = 0.05;
k     = 0.02;
alpha = -0.4;
tau   = 0.2;
rho   = 0.6;
sigma = 0.16;

%% Define shock distribution
Mu                = [1 0];
Sigma             = [0.05 0;0 sigma];
[model.e,model.w] = qnwnorm([5 5],Mu,Sigma^2);
model.funrand     = @(nrep) Mu(ones(nrep,1),:)+randn(nrep,2)*Sigma;

%% Pack model structure
% Model functions
model.func   = @sto4model;
%%
% Model parameters
model.params = {delta,r,k,alpha,tau,rho,sigma};

%% Define approximation space
% Approximation order
order         = [15; 15];
n             = prod(order);
%%
% Minimum and maximum values of the state variable grid
smin          = [min(model.e(:,1))*0.95; 0.4];
smax          = [1.6; 2.05];
%%
% Interpolation structure
interp.fspace = fundefn('spli',order,smin,smax);
%%
% State collocation nodes
s             = gridmake(funnode(interp.fspace));

%% Provide a first guess
xinit         = [zeros(n,1) max(min(s(:,1).^(1/alpha),s(:,2)+tau),s(:,2)-tau) zeros(n,2)];
interp.cz     = funfitxy(interp.fspace,s,max(min(ones(n,1),s(:,2).^rho+tau),s(:,2).^rho-tau));
interp.ch     = funfitxy(interp.fspace,s,max(min(s(:,1).^(1/alpha),s(:,2)+tau),s(:,2)-tau));

%% Solve for rational expectations
% Define options
options = struct('funapprox','expfunapprox',...
                 'simulmethod','solve',...
                 'stat',1);
%%
% If available, we will use successively different methods to find the solution.
%
% * Solve by Full Newton
if exist('mcppath','file')
  options.eqsolver  = 'path';
  options.reemethod = '1-step';
  tic
  recsSolveREE(interp,model,s,xinit,options);
  toc
  options.eqsolver = 'lmmcp';
  options.reemethod = 'iter';
end

%%
% * Solve by Newton-Krylov
if exist('KINSol','file')
  options.reesolver = 'kinsol';
  tic
  recsSolveREE(interp,model,s,xinit,options);
  toc
  options.reesolver = 'SA';
end

%%
% * Solve by successive approximation
tic
interp = recsSolveREE(interp,model,s,xinit,options);
toc

%% Simulate the model
recsSimul(model,interp,repmat([1 1],10,1),200,[],options);
subplot(2,3,1)
xlabel('Availability')
ylabel('Frequency')
subplot(2,3,2)
xlabel('World price')
ylabel('Frequency')
subplot(2,3,3)
xlabel('Storage')
ylabel('Frequency')
subplot(2,3,4)
xlabel('Price')
ylabel('Frequency')
subplot(2,3,5)
xlabel('Import')
ylabel('Frequency')
subplot(2,3,6)
xlabel('Export')
ylabel('Frequency')

%%
% Copyright (C) 2011-2012 Christophe Gouel
%
% Licensed under the Expat license, see <LICENSE.txt>
