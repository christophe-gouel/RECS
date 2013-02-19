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
% *Transition equations*
%
% $$K_{t}: K_{t}=a e^{Z_{t-1}}K_{t-1}^{\alpha}+\left(1-\delta\right)K_{t-1}-C_{t-1},$$
%
% $$Z_{t}: Z_{t}=\rho Z_{t-1}+\epsilon_{t}.$$

%% Writing the model
% The model is defined in a Yaml file: <gro1.txt gro1.yaml>.

%% Pack model structure
% Mean and standard deviation of the shocks
Mu                = 0;
sigma             = 0.007;

%%
% You generate the MATLAB model file and pack the model structure with the
% following command
model = recsmodelinit('gro1.yaml',...
                      struct('Mu',Mu,'Sigma',sigma^2,'order',5));

%%
% This command creates a MATLAB file, <gro1model.html gro1model.m>, containing
% the definition the model and all its Jacobians from the human readable file
% <gro1.txt gro1.yaml>.

%%
% It is also possible, if you are ready to write all the model's equations and
% jacobians in the file |gro1model.m|, to pack yourself the different elements
% inside a structure
[model.e,model.w] = qnwnorm(5,Mu,sigma^2);
model.funrand     = @(nrep) Mu+sigma*randn(nrep,1);
model.func        = @gro1model;
model.params      = gro1model('params');

%% Find deterministic steady-state
a             = model.params(1);
delta         = model.params(3);
[sss,xss,zss] = recsSS(model,[1 0],a-delta);

%% Define approximation space using Chebyshev polynomials
% Degree of approximation
order         = [10 10];
%%
% Limits of the state space
smin          = [0.85*sss(1) min(model.e)*4];
smax          = [1.15*sss(1) max(model.e)*4];
%%
[interp,s] = recsinterpinit(order,smin,smax,'cheb');

%% Find a first guess through the perfect foresight solution
[interp,x] = recsFirstGuess(interp,model,s,sss,xss,50);

%% Define options
% With high order Chebyshev polynomials, extrapolation outside the state space
% should not be allowed to prevent explosive values.
options = struct('reesolver','krylov',...
                 'extrapolate',0    ,...
                 'accuracy'   ,1);

%% Solve for rational expectations
interp = recsSolveREE(interp,model,s,x,options);

%% Simulate the model
[~,~,~,~,stat] = recsSimul(model,interp,sss(ones(1000,1),:),200,[],options);
subplot(2,2,1)
xlabel('Capital stock')
ylabel('Frequency')
subplot(2,2,2)
xlabel('Log of productivity')
ylabel('Frequency')
subplot(2,2,3)
xlabel('Consumption')
ylabel('Frequency')
