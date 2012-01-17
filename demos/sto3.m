%% STO3 Competitive storage with price-band backed by public storage
% This model is an extension of Miranda and Helmberger [1] to include a capacity
% constraint on the public stock level.

%% Model's structure
%
% *Response variables* Speculative storage ($S^{\mathrm{S}}$), Increase in
% public stock level ($\Delta S^{\mathrm{G}+}$), Decrease in public stock level
% ($\Delta S^{\mathrm{G}-}$), Planned production ($H$) and Price ($P$).
%
% *State variable* Availability ($A$) and Public stock ($S^{\mathrm{G}}$).
%
% *Shock* Productivity shocks ($\epsilon$).
%
% *Parameters* Unit physical storage cost ($k$), Depreciation share ($\delta$),
% Interest rate ($r$), Scale parameter for the production cost function ($h$),
% Inverse of supply elasticity ($\mu$), Demand price elasticity ($\alpha$),
% Floor price ($P^{\mathrm{F}}$), Ceiling price ($P^{\mathrm{c}}$) and Capacity
% constraint on public stock ($\bar{S}^{\mathrm{G}}$).
%
% *Equilibrium equations*
%
% $$S^{\mathrm{S}}_{t}: S^{\mathrm{S}}_{t}\ge 0 \quad \perp \quad P_{t}+k-\frac{1-\delta}{1+r}\mathrm{E}_{t}\left(P_{t+1}\right)\ge 0,$$
%
% $$H_{t}: \frac{1}{1+r}\mathrm{E}_{t}\left(P_{t+1}\epsilon_{t+1}\right)=h {H_{t}}^{\mu},$$
%
% $$P_{t}: A_{t}+\Delta S^{\mathrm{G}-}_{t}={P_{t}}^{\alpha}+S^{\mathrm{S}}_{t}+\Delta S^{\mathrm{G}+}_{t},$$
%
% $$\Delta S^{\mathrm{G}+}_{t}:0\le\Delta S^{\mathrm{G}+}_{t}\le\bar{S}^{\mathrm{G}}-\left(1-\delta\right)S^{\mathrm{G}}_{t-1} \quad \perp \quad P_{t}-P^{\mathrm{F}},$$
%
% $$\Delta S^{\mathrm{G}-}_{t}:0\le\Delta S^{\mathrm{G}-}_{t}\le\left(1-\delta\right)S^{\mathrm{G}}_{t-1} \quad \perp \quad P^{\mathrm{c}}-P_{t}.$$
%
% *Transition equation*
%
% $$A_{t}: A_{t}=\left(1-\delta\right)\left(S^{\mathrm{S}}_{t-1}+\right)+H_{t-1}\epsilon_{t},$$
%
% $$S^{\mathrm{G}}_{t}: S^{\mathrm{G}}_{t}= \left(1-\delta\right)S^{\mathrm{G}}_{t}+\Delta S^{\mathrm{G}+}_{t}-\Delta S^{\mathrm{G}-}_{t}.$$

%% Writing the model
% The model is defined in a Matlab file: <matlab:filetohelp('sto3model.txt')
% |sto3model.m|>.

%% Enter model parameters
delta = 0;
r     = 0.05;
alpha = -0.4;
k     = 0.02;
mu    = 10;
PF    = 0.9;
PC    = 1.1;
Sgbar = 0.4;

%% Define shock distribution
Mu                = 1;
sigma             = 0.05;
[model.e,model.w] = qnwnorm(5,Mu,sigma^2);
model.funrand     = @(nrep) Mu+sigma*randn(nrep,1);

%% Pack model structure
% Model functions
model.func   = @sto3model;
%%
% Model parameters
model.params = {delta,r,alpha,k,mu,PF,PC,Sgbar};

%% Define approximation space
% Approximation order
order         = [20; 20];
n             = prod(order);
%%
% Minimum and maximum values of the state variable grid
smin          = [0.74 0    ];
smax          = [1.4  Sgbar];
%%
% Interpolation structure
interp.fspace = fundefn('spli',order,smin,smax);
%%
% State collocation nodes
s             = gridmake(funnode(interp.fspace));

%% Provide a very simple first guess
xinit         = [zeros(n,1) ones(n,2) zeros(n,2)];
interp.cz     = funfitxy(interp.fspace,s,ones(n,2));
interp.cx     = funfitxy(interp.fspace,s,xinit);

%% Solve for rational expectations
% If |PATH| is installed, this will solved the model using |PATH| to solve the
% equilibrium equations:
if exist('mcppath','file')
  options.eqsolver = 'path';
  tic
  recsSolveREE(interp,model,s,xinit,options);
  toc
end

%%
% Solve the model with default options:
tic
interp = recsSolveREE(interp,model,s,xinit);
toc

%% Simulate the model
options.stat = 1;
recsSimul(model,interp,repmat([1 0],1000,1),200,[],options);
subplot(3,3,1)
xlabel('Availability')
ylabel('Frequency')
subplot(3,3,2)
xlabel('Public stock')
ylabel('Frequency')
subplot(3,3,3)
xlabel('Speculative storage')
ylabel('Frequency')
subplot(3,3,4)
xlabel('Planned production')
ylabel('Frequency')
subplot(3,3,5)
xlabel('Price')
ylabel('Frequency')
subplot(3,3,6)
xlabel('Increase in public stock level')
ylabel('Frequency')
subplot(3,3,7)
xlabel('Decrease in public stock level')
ylabel('Frequency')

%% References
%
% [1] Miranda, M. J. and Helmberger, P. G. (1988). The effects of commodity price
% stabilization programs. _The American Economic Review_, 78(1), 46-58.

%%
% Copyright (C) 2011-2012 Christophe Gouel
%
% Licensed under the Expat license, see <matlab:filetohelp('LICENSE.txt') LICENSE.txt>

