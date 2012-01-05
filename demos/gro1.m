%% GRO1 Stochastic growth model

%% Model's structure
%
% *Response variable* Consumption ($C$).
%
% *State variable* Capital stock ($K$), Log of productivity ($Z$).
%
% *Shock* Innovation to productivity ($\epsilon$).
%
% *Parameters* Capital depreciation rate ($\delta$), Discount factor ($\beta$),
% Elasticity of intertemporal substitution ($\tau$), Capital share ($\alpha$),
% Scale parameter ($a$).
%
% *Equilibrium equation*
%
% $$C_{t}: C_{t}^{-\tau}=\beta\mathrm{E}_{t}\left[C_{t+1}^{-\tau}\left(1-\delta+a \alpha e^{Z_{t+1}}K_{t+1}^{\alpha-1}\right)\right].$$
%
% *Transition equation*
%
% $$K_{t}: K_{t}=a e^{Z_{t-1}}K_{t-1}^{\alpha}+\left(1-\delta\right)K_{t-1}-C_{t-1},$$
%
% $$Z_{t}: Z_{t}=\rho Z_{t-1}+\epsilon_{t}.$$

%% Writing the model
% The model is defined in a yaml file: <matlab:filetohelp('gro1.txt') |gro1.yaml|>.

%% Pack model structure
% Mean and standard deviation of the shocks
Mu                = 0;
sigma             = 0.007;

%%
% If your RECS installation is complete, you can generate the Matlab model file
% and pack the model structure with the following command
model = recsmodelinit('gro1.yaml',...
                      struct('Mu',Mu,'Sigma',sigma^2,'order',5));

%%
% This command creates a Matlab file, <matlab:filetohelp('gro1model.txt')
% |gro1model.m|>, containing the definition the model and all its Jacobians from
% the human readable file <matlab:filetohelp('gro1.txt') |gro1.yaml|>.

%%
% If your installation is not complete, you have to pack yourself the different
% elements inside a structure
[model.e,model.w] = qnwnorm(5,Mu,sigma^2);
model.funrand     = @(nrep) Mu+sigma*randn(nrep,1);
model.func        = @gro1model;
model.params      = gro1model('params');

%% Find deterministic steady-state
a             = model.params(1);
delta         = model.params(3);
[sss,xss,zss] = recsSS(model,[1 0],a-delta)

%% Define approximation space
% Degree of approximation
order         = [10 10];                                 
%%
% Limits of the state space
smin          = [0.85*sss(1) min(model.e)*4];
smax          = [1.15*sss(1) max(model.e)*4];
%%
% Function space
interp.fspace = fundefn('cheb',order,smin,smax);           
%%
% State collocation nodes
s             = gridmake(funnode(interp.fspace));

%% Find a first guess through the perfect foresight solution
[interp,x] = recsFirstGuess(interp,model,s,sss,xss,50);

%% Define options
% With Chebyshev polynomials, extrapolation outside the state space should not
% be allowed.
options = struct('reesolver','mixed',...
                 'extrapolate',0);

%% Solve for rational expectations
interp = recsSolveREE(interp,model,s,x,options);

%% Simulate the model
[~,~,~,~,stat] = recsSimul(model,interp,sss(ones(1000,1),:),200);
subplot(2,2,1)
xlabel('Capital stock')
ylabel('Frequency')
subplot(2,2,2)
xlabel('Log of productivity')
ylabel('Frequency')
subplot(2,2,3)
xlabel('Consumption')
ylabel('Frequency')

%%
% Copyright (C) 2011-2012 Christophe Gouel
%
% Licensed under the Expat license, see <matlab:filetohelp('LICENSE.txt') LICENSE.txt>
