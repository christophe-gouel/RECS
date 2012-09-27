%% CS1 Consumption/saving model with borrowing constraint
% This is an implementation of the model in Deaton (1991).

%% Model's structure
%
% *Response variable* Consumption ($C$).
%
% *State variable* Cash on hand ($X$).
%
% *Shock* Labor income ($Y$).
%
% *Parameters* Interest rate ($r$), Rate of time preference ($\delta$), and
% Elasticity of intertemporal substitution ($\rho$).
%
% *Equilibrium equation*
%
% $$C_{t}: C_{t}\le X_{t} \quad \perp \quad \frac{1+r}{1+\delta}\mathrm{E}_{t}\left(C_{t+1}^{-\rho}\right)-C_{t}^{-\rho}\le 0.$$
%
% *Transition equation*
%
% $$X_{t}: X_{t}=\left(1+r\right)\left(X_{t-1}-C_{t-1}\right)+Y_{t}.$$

%% Writing the model
% The model is defined in a yaml file: <cs1.yaml>.

%% Pack model structure
% Mean and standard deviation of the shocks
Mu                = 100;
sigma             = 10;

%%
% If your RECS installation is complete, you can generate the Matlab model file
% and pack the model structure with the following command
model = recsmodelinit('cs1.yaml',...
                      struct('Mu',Mu,'Sigma',sigma^2,'order',5));

%%
% This command creates a Matlab file, <cs1model.html cs1model.m>, containing the
% definition the model and all its Jacobians from the human readable file
% <cs1.yaml>.

%% Define approximation space
[interp,s] = recsinterpinit(20,model.sss/2,model.sss*2);

%% First-guess: Consumption equal to cash on hand
x         = s;
%%
% To force the solver to compute the approximation of the expectations
% function, it is necessary to add at least an empty value for |interp.ch|
interp.ch = [];

%% Solve for rational expectations
[interp,x] = recsSolveREE(interp,model,s,x);

%% Plot the decision rule
figure
plot(s,x,s,s)
legend('Policy rule','45 degree line')
legend('Location','NorthWest')
legend('boxoff')
xlabel('Cash on hand')
ylabel('Consumption')

%% Simulate the model
[~,~,~,~,stat] = recsSimul(model,interp,model.sss(ones(1000,1),:),200);
subplot(1,2,1)
xlabel('Cash on hand')
ylabel('Frequency')
subplot(1,2,2)
xlabel('Consumption')
ylabel('Frequency')

%% References
%
% Deaton, A. (1991). Saving and liquidity constraints. _Econometrica_, 59(5),
% 1221-1248.

%%
% Copyright (C) 2011-2012 Christophe Gouel
%
% Licensed under the Expat license, see <LICENSE.txt>
