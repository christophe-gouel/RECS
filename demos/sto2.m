%% STO2 Competitive storage with floor-price backed by public storage
% This model is close to one of the models presented in Wright and Williams (1988).

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
% $$P_{t}: A_{t}={P_{t}}^{\alpha}+S^{\mathrm{S}}_{t}+S^{\mathrm{G}}_{t},$$
%
% $$S^{\mathrm{G}}_{t}: 0\le S^{\mathrm{G}}_{t}\le \bar{S}^{\mathrm{G}} \quad \perp \quad P_t - P^{\mathrm{F}}.$$
%
% *Transition equation*
%
% $$A_{t}: A_{t}=\left(1-\delta\right)\left(S^{\mathrm{S}}_{t-1}+S^{\mathrm{G}}_{t-1}\right)+H_{t-1}\epsilon_{t}.$$

%% Writing the model
% The model is defined in a Matlab file: <sto2model.html sto2model.m>.

%% Pack model structure
Mu                = 1;
sigma             = 0.05;

%%
model = recsmodelinit('sto2.yaml',struct('Mu',Mu,'Sigma',sigma^2,'order',7));

%% Define approximation space
[interp,s] = recsinterpinit(50,model.sss*0.7,model.sss*1.6);

%% Find a first guess through the perfect foresight solution
interp = recsFirstGuess(interp,model,s,model.sss,model.xss,5);

%% Solve for rational expectations
[interp,x] = recsSolveREE(interp,model,s);

%% Plot storage rules
subplot(2,1,1)
plot(s,x(:,[1 4]))
leg = legend('Private stock','Public stock');
set(leg,'Location','NorthWest')
set(leg,'Box','off')
xlabel('Availability')
ylabel('Stock')
subplot(2,1,2)
plot(s,x(:,3))
xlabel('Availability')
ylabel('Price')

%% References
%
% Wright, B. D. and Williams, J. C. (1988). The incidence of market-stabilising
% price support schemes. _The Economic Journal_, 98(393), 1183-1198.
