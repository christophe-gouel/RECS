%% STO4 Competitive storage with price-band backed by public storage
% This model is an extension of Miranda and Helmberger (1988) to include a
% capacity constraint on the public stock level.

%% Model's structure
%
% *Response variables* Speculative storage ($S^{\mathrm{S}}$), Increase in
% public stock level ($\Delta S^{\mathrm{G}+}$), Decrease in public stock level
% ($\Delta S^{\mathrm{G}-}$), Planned production ($H$) and Price ($P$).
%
% *State variable* Availability ($A$) and Public stock ($S^{\mathrm{G}}$).
%
% *Shock* Productivity shocks ($\epsilon$).
%
% *Parameters* Unit physical storage cost ($k$), Interest rate ($r$), Scale
% parameter for the production cost function ($h$), Inverse of supply elasticity
% ($\mu$), Demand price elasticity ($\alpha$), Floor price ($P^{\mathrm{F}}$),
% Ceiling price ($P^{\mathrm{c}}$) and Capacity constraint on public stock
% ($\bar{S}^{\mathrm{G}}$).
%
% *Equilibrium equations*
%
% $$S^{\mathrm{S}}_{t}: S^{\mathrm{S}}_{t}\ge 0 \quad \perp \quad P_{t}+k-\frac{\mathrm{E}_{t}\left(P_{t+1}\right)}{1+r}\ge 0,$$
%
% $$H_{t}: \frac{1}{1+r}\mathrm{E}_{t}\left(P_{t+1}\epsilon_{t+1}\right)=h {H_{t}}^{\mu},$$
%
% $$P_{t}: A_{t}+\Delta S^{\mathrm{G}-}_{t}={P_{t}}^{\alpha}+S^{\mathrm{S}}_{t}+\Delta S^{\mathrm{G}+}_{t},$$
%
% $$\Delta S^{\mathrm{G}+}_{t}:0\le\Delta S^{\mathrm{G}+}_{t}\le\bar{S}^{\mathrm{G}}-S^{\mathrm{G}}_{t-1} \quad \perp \quad P_{t}-P^{\mathrm{F}},$$
%
% $$\Delta S^{\mathrm{G}-}_{t}:0\le\Delta S^{\mathrm{G}-}_{t}\le S^{\mathrm{G}}_{t-1} \quad \perp \quad P^{\mathrm{c}}-P_{t}.$$
%
% *Transition equation*
%
% $$A_{t}: A_{t}=S^{\mathrm{S}}_{t-1}+H_{t-1}\epsilon_{t},$$
%
% $$S^{\mathrm{G}}_{t}: S^{\mathrm{G}}_{t}= S^{\mathrm{G}}_{t}+\Delta S^{\mathrm{G}+}_{t}-\Delta S^{\mathrm{G}-}_{t}.$$

%% Writing the model
% The model is defined in a Yaml file: <sto4.txt sto4.yaml>.

%% Pack model structure
Mu                = 1;
sigma             = 0.05;
model = recsmodelinit('sto4.yaml',struct('Mu',Mu,'Sigma',sigma^2,'order',7));
%%
Sgbar = model.params(end-1);

%% Multiple steady states
% In this model, there is no stock depreciation. This assumption implies that
% there are multiple steady states: as long as the steady-state price is between
% the floor and ceiling prices, any public stock level between 0 and
% $\bar{S}^{\mathrm{G}}$ is a steady state. Actually, the steady-state response
% variables are unique, only the public stock level is indeterminate as we can
% see in the examples below:
[sss1,xss1] = recsSS(model,[1 0],model.xss)
%%
[sss2,xss2] = recsSS(model,[1 0.2],model.xss)

%% Define approximation space
[interp,s] = recsinterpinit(20,[0.74 0],[1.4 Sgbar]);

%% Find a first guess through the perfect foresight solution
interp = recsFirstGuess(interp,model,s,model.sss,model.xss,5);

%% Solve for rational expectations
[interp,x] = recsSolveREE(interp,model,s);

%% Simulate the model
options.stat = 1;
recsSimul(model,interp,repmat([1 0],1E3,1),200,[],options);
subplot(3,3,1)
xlabel('Availability')
ylabel('Frequency')
subplot(3,3,2)
xlabel('Public stock')
ylabel('Frequency')
subplot(3,3,3)
xlabel('Speculative storage')
ylabel('Frequency')
subplot(3,3,4)
xlabel('Planned production')
ylabel('Frequency')
subplot(3,3,5)
xlabel('Price')
ylabel('Frequency')
subplot(3,3,6)
xlabel('Increase in public stock level')
ylabel('Frequency')
subplot(3,3,7)
xlabel('Decrease in public stock level')
ylabel('Frequency')

%% References
%
% Miranda, M. J. and Helmberger, P. G. (1988). The effects of commodity price
% stabilization programs. _The American Economic Review_, 78(1), 46-58.
