%% STO2 Competitive storage with floor-price backed by public storage
% This model is close to one of the models presented in Wright and Williams [1].

%% Model's structure
%
% *Response variables* Speculative storage ($S^{\mathrm{S}}$), Public storage
% ($S^{\mathrm{G}}$), Planned production ($H$) and Price ($P$).
%
% *State variable* Availability ($A$).
%
% *Shock* Productivity shocks ($\epsilon$).
%
% *Parameters* Unit physical storage cost ($k$), Depreciation share ($\delta$),
% Interest rate ($r$), Scale parameter for the production cost function ($h$),
% Inverse of supply elasticity ($\mu$), Demand price elasticity ($\alpha$),
% Floor price ($P^{\mathrm{F}}$) and Capacity constraint on public stock
% ($\bar{S}^{\mathrm{G}}$).
%
% *Equilibrium equations* 
%
% $$S^{\mathrm{S}}_{t}: S^{\mathrm{S}}_{t}\ge 0 \quad \perp \quad P_{t}+k-\frac{1-\delta}{1+r}\mathrm{E}_{t}\left(P_{t+1}\right)\ge 0,$$
%
% $$H_{t}: \frac{1}{1+r}\mathrm{E}_{t}\left(P_{t+1}\epsilon_{t+1}\right)=h {H_{t}}^{\mu},$$
%
% $$P_{t}: A_{t}={P_{t}}^{\alpha}+S^{\mathrm{S}}_{t}+S^{\mathrm{G}}_{t}.$$
%
% $$S^{\mathrm{G}}_{t}: 0\le S^{\mathrm{G}}_{t}\le \bar{S}^{\mathrm{G}} \quad \perp \quad P_t - P^{\mathrm{F}}$$
%
% *Transition equation*
%
% $$A_{t}: A_{t}=\left(1-\delta\right)\left(S^{\mathrm{S}}_{t-1}+S^{\mathrm{G}}_{t-1}\right)+H_{t-1}\epsilon_{t}.$$

%% Writing the model
% The model is defined in a Matlab file: <matlab:filetohelp('sto2model.txt')
% |sto2model.m|>.

%% Enter model parameters
k     = 0.02;
delta = 0.02;
r     = 0.05;
mu    = 10;
alpha = -0.4;
PF    = 1.02;
Sgbar = 0.4;

%% Define shock distribution
Mu                = 1;
sigma             = 0.05;
[model.e,model.w] = qnwnorm(5,Mu,sigma^2);
model.funrand     = @(nrep) Mu+sigma*randn(nrep,1);

%% Pack model structure
% Model functions
model.func   = @sto2model;                               
%%
% Model parameters
model.params = {k,delta,r,mu,alpha,PF,Sgbar};

%% Define approximation space
% Approximation order
order         = 50;                                        
n             = 30;  
%%
% Minimum and maximum values of the state variable grid
smin          = min(model.e)*0.95;
smax          = 2;
%%
% Interpolation structure
interp.fspace = fundefn('spli',order,smin,smax); 
%%
% State collocation nodes
s             = gridmake(funnode(interp.fspace));          

%% Find a first guess through the perfect foresight solution
[interp,xinit] = recsFirstGuess(interp,model,s,1,[0 1 1 0],5);

%% Solve for rational expectations
% We will use successively different methods to find the solution.
%
% * Default options (successive approximations of response variables):
tic
recsSolveREE(interp,model,s,xinit);
toc

%%
% * Response variables approximation applied to all response variables in expectations
options.funapprox = 'resapprox-simple';
tic
recsSolveREE(interp,model,s,xinit,options);
toc

%%
% * Expectations approximation
options.funapprox = 'expapprox';
tic
recsSolveREE(interp,model,s,xinit,options);
toc

%%
% * Solve equilibrium equations with |ncpsolve|
options.eqsolver = 'ncpsolve';
optset('ncpsolve','type','minmax'); % 'minmax' / 'smooth'
tic
interp = recsSolveREE(interp,model,s,xinit,options);
toc

%% Simulate the model
% Use high precision to be able to draw precisely decision rules 
options.simulmethod = 'solve';
[ssim,xsim] = recsSimul(model,interp,ones(20,1),200,[],options);

%% Plot storage rules
subplot(2,1,1)
plot(ssim(:),reshape(xsim(:,1,:),[],1),'.',...
     ssim(:),reshape(xsim(:,4,:),[],1),'.')
leg = legend('Private stock','Public stock');
set(leg,'Location','NorthWest')
set(leg,'Box','off')
xlabel('Availability')
ylabel('Stock')
subplot(2,1,2)
plot(ssim(:),reshape(xsim(:,3,:),[],1),'.')
xlabel('Availability')
ylabel('Price')

%% References
%
% [1] Wright, B. D. and Williams, J. C. (1988). The incidence of market-stabilising
% price support schemes. _The Economic Journal_, 98(393), 1183-1198.

%%
% Copyright (C) 2011-2012 Christophe Gouel
%
% Licensed under the Expat license, see <matlab:filetohelp('LICENSE.txt') LICENSE.txt>

