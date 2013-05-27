%% GRO3 Stochastic growth model with recursive preferences and stochastic volatility (SV)
% This model is similar to Caldara et al. (2012).

%% Model's structure
% *Response variable* Consumption ($C$), Labor ($L$), Marginal utility with
% respect to consumption ($U_C$), Instantaneous utility ($U$), Production ($Y$),
% Utility ($V$).
%
% *State variable* Capital stock ($K$), Log of productivity ($Z$), Log of
% productivity volatility ($\sigma$).
%
% *Shock* Innovation to productivity ($\epsilon$), Innovation to productivity
% volatility ($\omega$).
%
% *Parameters* Capital depreciation rate ($\delta$), Discount factor ($\beta$),
% Index of deviation from CRRA utility ($\nu$), Risk-aversion parameter
% ($\tau$), Capital share ($\alpha$), Share of consumption in utility
% ($\theta$), Persistence of log-productivity ($\rho_Z$), Persistence of SV
% ($\rho_{\sigma}$), Unconditional mean of SV ($\bar{\sigma}$), Standard
% deviation of SV ($\eta$).
%
% *Equilibrium equation*
%
% $$C_{t}: U_{C,t}=\beta\left(\mathrm{E}_{t}\tilde{V}_{t+1}\right)^{\left(\frac{1}{\theta}-1\right)}\mathrm{E}_{t}\left[\tilde{V}_{t+1}^{\frac{\theta-1}{\theta}}U_{C,t+1}\left(1-\delta+ \alpha e^{Z_{t+1}}K_{t+1}^{\alpha-1}L_{t+1}^{\alpha-1}\right)\right].$$
%
% $$L_t: (1-\theta) C_t = \theta (1-\alpha) Y_t \left(\frac{1}{L_t}-1\right)$$
%
% $$U_{C,t}: U_{C,t} = \theta \frac{1-\tau}{\nu} \frac{U_t^{1/\nu}}{C_t}$$
%
% $$U_t: U_t = \left[C_t^{\theta}(1-L_t)^{1-\theta}\right]^{1-\tau}$$
%
% $$Y_t: Y_t = e^{Z_{t}}K_{t}^{\alpha}L_t^{1-\alpha}$$
%
% $$V_t: \tilde{V}_t^{\frac{1}{\nu}} = (1-\beta) U_t^{\frac{1}{\nu}} + \beta \mathrm{E}_t \tilde{V}_t^{\frac{1}{\nu}}$$
%
% $$\tilde{V}_t: \tilde{V}_t=V_t^{1-\tau}$$
%
% *Transition equations*
%
% $$K_{t}: K_{t}=+\left(1-\delta\right)K_{t-1}-C_{t-1},$$
%
% $$Z_{t}: Z_{t}=\rho_Z Z_{t-1}+e^{\sigma_t}\epsilon_{t}.$$
%
% $$\sigma_t: \sigma_t = (1-\rho_{\sigma})\bar{\sigma}+\rho_{\sigma}\sigma_{t-1}+\eta\omega_t$$

%% Writing the model
% The model is defined in a Yaml file: <gro3.txt gro3.yaml>.

%% Create the model object
model = recsmodel('gro3.yaml',struct('Mu',[0 0],'Sigma',eye(2),'order',5));


%% Define approximation space using Chebyshev polynomials
smin       = [0.85*model.sss(1) -0.11 log(0.007)*1.15];
smax       = [1.20*model.sss(1)  0.11 log(0.007)*0.85];
%%
[interp,s] = recsinterpinit(4,smin,smax,'cheb');

%% Find a first guess through the perfect foresight solution
[interp,x] = recsFirstGuess(interp,model,s);

%% Define options
options = struct('reemethod','1-step',...
                 'accuracy' ,1,...
                 'stat'     ,1);

%% Solve for rational expectations
[interp,x,z] = recsSolveREE(interp,model,s,x,options);

%% Use simple continuation method to solve for higher values of risk aversion
% The procedure to solve for different parameters values has to be packed in a
% function. This is done in gro3problem.m.
type('gro3problem.m')
%%
% This function requires as input the cell array X:
X = {model interp s x options};
%%
% The function SCP starts from the known solution with a low risk aversion to
% find in 2 steps the solution with an higher risk aversion:
X = SCP(X,[5 0.5],[0.5 1/0.5],@gro3problem,2);
[model,interp,s,x,options] = X{:};

%% Simulate the model
recsSimul(model,interp,model.sss(ones(100,1),:),200,[],options);
subplot(3,4,1)
xlabel('Capital stock')
ylabel('Frequency')
subplot(3,4,2)
xlabel('Log of productivity')
ylabel('Frequency')
subplot(3,4,3)
xlabel('Log of productivity volatility')
ylabel('Frequency')
subplot(3,4,4)
xlabel('Consumption')
ylabel('Frequency')
subplot(3,4,5)
xlabel('Labor')
ylabel('Frequency')
subplot(3,4,6)
xlabel('Marginal utility wrt consumption')
ylabel('Frequency')
subplot(3,4,7)
xlabel('Instantaneous utility')
ylabel('Frequency')
subplot(3,4,8)
xlabel('Production')
ylabel('Frequency')
subplot(3,4,9)
xlabel('Intertemporal utility')
ylabel('Frequency')
subplot(3,4,10)
xlabel('')
ylabel('Frequency')

%% References
% <http://dx.doi.org/10.1016/j.red.2011.10.001 Caldara D.;
% Fernandez-Villaverdes, J.; Rubio-Ramirez, J. F. & Yao, W. (2012). Computing
% DSGE models with recursive preferences and stochastic volatility. _Review of
% Economic Dynamics_ 15(2), 188-206.>