%% STO1 Competitive storage model
% This model is a standard competitive storage model with supply reaction. Its
% setting is close to Wright and Williams (1982).

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
% The model is defined in a Yaml file: <sto1.txt sto1.yaml>.

%% Create the model object
Mu                = 1;
sigma             = 0.05;

%%
model = recsmodel('sto1.yaml',struct('Mu',Mu,'Sigma',sigma^2,'order',7));

%% Define approximation space
[interp,s] = recsinterpinit(40,model.sss*0.7,model.sss*1.5);

%% Find a first guess through the perfect foresight solution
interp = recsFirstGuess(interp,model,s,model.sss,model.xss,struct('T',5));

%% Solve for rational expectations
interp = recsSolveREE(interp,model);

%% Plot the decision rules
recsDecisionRules(model,interp,[],[],[],struct('simulmethod','solve'));
subplot(2,2,1)
xlabel('Availability')
ylabel('Storage')
subplot(2,2,2)
xlabel('Availability')
ylabel('Planned production')
subplot(2,2,3)
xlabel('Availability')
ylabel('Price')

%% Simulate the model
[ssim,~,~,stat] = recsSimul(model,interp,model.sss(ones(1000,1),:),200);
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
recsAccuracy(model,interp,ssim);

%% References
% <http://www.jstor.org/stable/2232552 Wright, B. D. and Williams,
% J. C. (1982). The economic role of commodity storage. _The Economic Journal_,
% 92(367), 596-614.>
