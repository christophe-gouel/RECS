%% STO5 Two-country storage-trade model
% Model close to Miranda and Glauber [1]. Countries $a$ and $b$ are indicated by
% the subscript $i \in \left\{a,b\right\}$. When variables corresponding to both
% countries appear in the same equation, the foreign country is indicated by the
% subscript $j$.

%% Model's structure
%
% *Response variables* Storage ($S_{i}$), Planned production ($H_{i}$), Price
% ($P_{i}$) and Export ($X_{i}$).
%
% *State variable* Availability ($A_{i}$).
%
% *Shock* Productivity shocks ($\epsilon_{i}$).
%
% *Parameters* Unit physical storage cost ($k$), Depreciation share ($\delta$),
% Interest rate ($r$), Scale parameter for the production cost function ($h$),
% Inverse of supply elasticity ($\mu$), Demand price elasticity ($\alpha$),
% Trade cost ($\theta$) and Tariff ($\tau_{i}$).
%
% *Equilibrium equations* 
%
% For $i \in \left\{a,b\right\}$ and $j \ne i$:
%
% $$S_{it}: S_{it}\ge 0 \quad \perp \quad P_{it}+k-\frac{1-\delta}{1+r}\mathrm{E}_{t}\left(P_{it+1}\right)\ge 0,$$
%
% $$H_{it}: \frac{1}{1+r}\mathrm{E}_{t}\left(P_{it+1}\epsilon_{it+1}\right)=h {H_{it}}^{\mu},$$
%
% $$P_{it}: A_{it}+X_{jt}={P_{it}}^{\alpha}+S_{it}+X_{it},$$
%
% $$X_{it}: X_{it}\ge 0 \quad \perp \quad P_{it}+\theta+\tau_{j}\ge P_{jt}.$$
%
% *Transition equation*
%
% For $i \in \left\{a,b\right\}$:
%
% $$A_{it}: A_{it}=\left(1-\delta\right)S_{it-1}+H_{it-1}\epsilon_{it}.$$

%% Writing the model
% The model is defined in a Matlab file: <matlab:filetohelp('sto5model.txt')
% |sto5model.m|>.

%% Enter model parameters
delta  = 0;
r      = 0.05;
k      = 0.02;
alpha  = -0.4;
theta  = 0.1;
%%
% Relative demand size, equals $D(b)/D(a)$ for the same price
lambda = 1;
%%
% Relative supply size, such that the unit production cost of eta $Q$ in $b$ is
% the same as $Q$ in $a$
eta    = 1;  
mu     = 10;
%%
% Tariffs in countries $a$ and $b$:
taua   = 0.06;
taub   = 0;

%% Define shock distribution
% Mean of production shocks in $a$ and $b$
Mu                = [1 1];
%%
% Std-deviation of production shocks in $a$ and $b$
sigma             = [0.05 0;
                     0    0.05];   
[model.e,model.w] = qnwnorm([5 5],Mu,sigma^2);
model.funrand     = @(nrep) Mu(ones(nrep,1),:)+randn(nrep,2)*sigma;

%% Pack model structure
% Model function
model.func   = @sto5model;                                    
%%
% Model parameters
model.params = {delta,r,k,alpha,theta,lambda,eta,mu,taua,taub};       

%% Find deterministic steady state
[sss,xss,zss] = recsSS(model,[1 1],[0 0 1 1 1 1 0 0])

%%
% *Check derivatives on steady state*
%
% Since model derivatives have been calculated by hand, it is important to check
% that there is no error.
recsCheck(model,sss,xss,zss);

%% Define approximation space
% Approximation order
order         = [15 15];
n             = prod(order);
%%
% Minimum and maximum values of the state variable grid (use the steady state to
% find good bounds)
smin          = sss*0.73;
smax          = sss*2;
%%
% Interpolation structure
interp.fspace = fundefn('spli',order,smin,smax);
%%
% State collocation nodes
s             = gridmake(funnode(interp.fspace));                      

%% Provide a first guess
xinit         = [0.1*zeros(n,2) max(s(:,1).^(1/alpha),0.8) ...
                 max((s(:,2)/lambda).^(1/alpha),0.8) zeros(n,2) ones(n,2)];
interp.cz     = funfitxy(interp.fspace,s,ones(n,4));
interp.cx     = funfitxy(interp.fspace,s,xinit);
interp.ch     = funfitxy(interp.fspace,s,ones(n,4));

%% Solve for rational expectations
% Define options:
%
% * Use |ncpsolve| to solve equilibrium equations.
%
% * Reduce the step of successive iterations to 0.5.
%
% * Do not use the approximation structure to find the new guess of response variables.
options =struct('eqsolver','ncpsolve',...
                'reesolveroptions',struct('lambda',0.5),...
                'useapprox',0);
%%
% * Approach the MCP problem by smooth transformation
optset('ncpsolve','type','smooth')

%%
% *Solve by Full Newton*
tic
options.reemethod = '1-step';
recsSolveREE(interp,model,s,xinit,options);
toc

%%
% *Solve by successive approximations* (this is the default option)
tic
options.reemethod = 'iter';
interp = recsSolveREE(interp,model,s,xinit,options);
toc

%% Simulate the model
[~,~,~,~,stat] = recsSimul(model,interp,sss(ones(100,1),:),1000,[],options);

%% References
%
% [1] Miranda, M. J. and Glauber, J. W. (1995). Solving stochastic models of
% competitive storage and trade by Chebychev collocation methods. _Agricultural
% and Resource Economics Review_, 24(1), 70-77.

%%
% Copyright (C) 2011-2012 Christophe Gouel
%
% Licensed under the Expat license, see <matlab:filetohelp('LICENSE.txt') LICENSE.txt>

