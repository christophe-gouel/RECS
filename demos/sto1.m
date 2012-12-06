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

%% Pack model structure
Mu                = 1;
sigma             = 0.05;

%%
model = recsmodelinit('sto1.yaml',struct('Mu',Mu,'Sigma',sigma^2,'order',7));

%% Define approximation space
[interp,s] = recsinterpinit(40,model.sss*0.7,model.sss*1.5);

%% Find a first guess through the perfect foresight solution
interp = recsFirstGuess(interp,model,s,model.sss,model.xss,5);

%% Solve for rational expectations
[interp,x] = recsSolveREE(interp,model,s);

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
recsAccuracy(model,interp,ssim);

%% References
%
% Wright, B. D. and Williams, J. C. (1982). The economic role of commodity
% storage. _The Economic Journal_, 92(367), 596-614.
