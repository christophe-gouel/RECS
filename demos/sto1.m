%% STO1 Competitive storage model
% This model is a standard competitive storage model with supply reaction. Its
% setting is close to Wright and Williams [1].

%% Model's structure
%
% *Response variables* Storage ($S$), Planned production ($H$) and Price ($P$).
%
% *State variable* Availability ($A$).
%
% *Shock* Productivity shocks ($\epsilon$).
%
% *Parameters* Unit physical storage cost ($k$), Depreciation share ($\delta$),
% Interest rate ($r$), Scale parameter for the production cost function ($h$),
% Inverse of supply elasticity ($\mu$) and Demand price elasticity ($\alpha$).
%
% *Equilibrium equations*
%
% $$S_{t}: S_{t}\ge 0 \quad \perp \quad P_{t}+k-\frac{1-\delta}{1+r}\mathrm{E}_{t}\left(P_{t+1}\right)\ge 0,$$
%
% $$H_{t}: \frac{1}{1+r}\mathrm{E}_{t}\left(P_{t+1}\epsilon_{t+1}\right)=h {H_{t}}^{\mu},$$
%
% $$P_{t}: A_{t}={P_{t}}^{\alpha}+S_{t}.$$
%
% *Transition equation*
%
% $$A_{t}: A_{t}=\left(1-\delta\right)S_{t-1}+H_{t-1}\epsilon_{t}.$$

%% Writing the model
% The model is defined in a Matlab file: <matlab:filetohelp('sto1model.txt') |sto1model.m|>.

%% Enter model parameters
alpha = -0.4;
k     = 0.02;
delta = 0.02;
r     = 0.05;
mu    = 10;

%% Define shock distribution
Mu                = 1;
sigma             = 0.1;
[model.e,model.w] = qnwnorm(5,Mu,sigma^2);
model.funrand     = @(nrep) Mu+sigma*randn(nrep,1);

%% Pack model structure
% Model function
model.func   = @sto1model;
%%
% Model parameters
model.params = {alpha,k,delta,r,mu};

%% Define approximation space
% Approximation order
n             = 30;
%%
% Minimum and maximum values of the state variable grid
smin          = 0.5;
smax          = 2;
%%
% Interpolation structure
interp.fspace = fundefn('spli',n,smin,smax);
%%
% State collocation nodes
s             = gridmake(funnode(interp.fspace));

%% Provide a first guess
xinit     = [zeros(n,1) ones(n,1) s.^(1/alpha)];
interp.cz = funfitxy(interp.fspace,s,ones(n,2));
interp.ch = funfitxy(interp.fspace,s,[s.^(1/alpha) s.^(1/alpha)]);

%% Find deterministic steady state
[sss,xss,zss] = recsSS(model,1,[0 1 1])

%%
% *Check derivatives on steady state*
%
% Since model derivatives have been calculated by hand, it is important to check
% that there is no error.
recsCheck(model,sss,xss,zss);

%% Solve for rational expectations
% We will use successively different methods to find the solution.
%
% * Default options (successive approximations of response variables):
tic
[~,x] = recsSolveREE(interp,model,s,xinit);
toc

%%
% * Expectations approximation - Newton-Krylov iterations
options = struct(...
    'funapprox','expapprox',...
    'reesolver','krylov');
tic
recsSolveREE(interp,model,s,xinit,options);
toc

%%
% * Expectations function approximation - Newton-Krylov iterations
options.funapprox = 'expfunapprox';
tic
recsSolveREE(interp,model,s,xinit,options);
toc

%%
% * Response variables approximation (used both for next- and current-period response) - Newton-Krylov iterations
options.funapprox = 'resapprox-simple';
tic
interp = recsSolveREE(interp,model,s,xinit,options);
toc

%%
% * Response variables approximation (used both for next- and current-period response) - Newton iterations (slow because it uses numerical derivatives)
if exist('fsolve','file')
  options.reesolver = 'fsolve';
  interp.cz     = ones(n,2);
  interp.cx     = xinit;
  tic
    interp = recsSolveREE(interp,model,s,xinit,options);
  toc
end

%%
% * Response variables approximation (used both for next- and current-period response) - Newton-Krylov iterations
if exist('kinsol','file')
  options.reesolver = 'kinsol';
  interp.cz     = ones(n,2);
  interp.cx     = xinit;
  tic
    interp = recsSolveREE(interp,model,s,xinit,options);
  toc
end

%% Plot the decision rules
figure
subplot(1,3,1)
plot(s,x(:,1))
title('Storage')
subplot(1,3,2)
plot(s,x(:,2))
title('Planned production')
subplot(1,3,3)
plot(s,x(:,3))
title('Price')

%% Simulate the model
reset(RandStream.getDefaultStream);
[ssim,~,~,~,stat] = recsSimul(model,interp,ones(1000,1),200);
subplot(2,2,1)
xlabel('Availability')
ylabel('Frequency')
subplot(2,2,2)
xlabel('Storage')
ylabel('Frequency')
subplot(2,2,3)
xlabel('Planned production')
ylabel('Frequency')
subplot(2,2,4)
xlabel('Price')
ylabel('Frequency')

%% Assess solution accuracy
[se,lEE] = recsAccuracy(model,interp,ssim);

%%
% Plot the Euler equation error (see Gouel [2] for definition)
figure
[~,i] = sort(se);
plot(se(i),lEE(i,1))
xlim([0.6 1.5])
ylim([-10 -2])
title('Euler equation error on storage-arbitrage equation')
xlabel('Availability')
ylabel(texlabel('log_{10}|EE|'))

%% References
%
% [1] Wright, B. D. and Williams, J. C. (1982). The economic role of commodity
% storage. _The Economic Journal_, 92(367), 596-614.
%
% [2] Gouel, C. (forthcoming). Comparing numerical methods for solving the competitive
% storage model. _Computational Economics_.

%%
% Copyright (C) 2011-2012 Christophe Gouel
%
% Licensed under the Expat license, see <matlab:filetohelp('LICENSE.txt') LICENSE.txt>
