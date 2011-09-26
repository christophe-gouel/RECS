%% CS1 Consumption/Saving model with borrowing constraint
% This is an implementation of the model in Deaton [1].

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
% The model is defined in a yaml file: <matlab:filetohelp('cs1.yaml') |cs1.yaml|>.

%% Pack model structure
Mu                = 100;
sigma             = 10;

%%
% If your RECS installation is complete, you can generate the Matlab model file
% and pack the model structure with the following command
model = recsmodelinit('cs1.yaml',...
                      struct('Mu',Mu,'Sigma',sigma^2,'order',5));

%%
% This command creates a Matlab file, <matlab:filetohelp('cs1model.m') |cs1model.m|>,
% containing the definition the model and all its Jacobian from the human readable file
% <matlab:filetohelp('cs1.yaml') |cs1.yaml|>.

%%
% If not, you have to pack yourself the different elements inside a structure
[model.e,model.w] = qnwnorm(5,Mu,sigma^2);
model.funrand     = @(nrep) Mu+sigma*randn(nrep,1);
model.func        = @cs1model;
model.params      = cs1model('params');


%% Define approximation space
% Define approximation function
interp.fspace = fundefn('spli',20,50,200);
%%
% State collocaton nodes
s             = gridmake(funnode(interp.fspace));

%% First-guess: Consumption equal to cash on hand
x         = s;
interp.cx = funfitxy(interp.fspace,s,x);
%%
% To force the solver to compute the approximation of the expectations
% function, it is necessary to add at least an empty value for |interp.ch|
interp.ch = [];

%% Solve for rational expectations
[interp,x] = recsSolveREE(interp,model,s,x);

%% Plot the decision rule
plot(s,x)
xlabel('Cash on hand')
ylabel('Consumption')

%% Simulate the model
[~,~,~,~,stat] = recsSimul(model,interp,100*ones(1000,1),200);
subplot(1,2,1)
xlabel('Cash on hand')
ylabel('Frequency')
subplot(1,2,2)
xlabel('Consumption')
ylabel('Frequency')

%% References
%
% [1] Deaton, A. (1991). Saving and liquidity constraints. _Econometrica_,
% 59(5), 1221-1248.

%%
% Copyright (C) 2011 Christophe Gouel
%
% Licensed under the Expat license, see <matlab:filetohelp('LICENSE.txt') LICENSE.txt>
