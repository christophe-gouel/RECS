%% STO6 Two-country storage-trade model
% Model close to Miranda and Glauber (1995). Countries $a$ and $b$ are indicated
% by the subscript $i \in \left\{a,b\right\}$. When variables corresponding to
% both countries appear in the same equation, the foreign country is indicated
% by the subscript $j$.

%% Model's structure
%
% *Response variables* Storage ($S_{i}$), Planned production ($H_{i}$), Price
% ($P_{i}$) and Export ($X_{i}$).
%
% *State variable* Availability ($A_{i}$).
%
% *Shock* Productivity shocks ($\epsilon_{i}$).
%
% *Parameters* Unit physical storage cost ($k$), Depreciation share ($\delta$),
% Interest rate ($r$), Scale parameter for the production cost function ($h$),
% Inverse of supply elasticity ($\mu$), Demand price elasticity ($\alpha$),
% Scale parameter for demand function ($\gamma_i$), Trade cost ($\theta$) and
% Tariff ($\tau_{i}$).
%
% *Equilibrium equations*
%
% For $i \in \left\{a,b\right\}$ and $j \ne i$:
%
% $$S_{it}: S_{it}\ge 0 \quad \perp \quad P_{it}+k-\frac{1-\delta}{1+r}\mathrm{E}_{t}\left(P_{it+1}\right)\ge 0,$$
%
% $$H_{it}: \frac{1}{1+r}\mathrm{E}_{t}\left(P_{it+1}\epsilon_{it+1}\right)=h {H_{it}}^{\mu},$$
%
% $$P_{it}: A_{it}+X_{jt}=\gamma_i{P_{it}}^{\alpha}+S_{it}+X_{it},$$
%
% $$X_{it}: X_{it}\ge 0 \quad \perp \quad P_{it}+\theta+\tau_{j}\ge P_{jt}.$$
%
% *Transition equation*
%
% For $i \in \left\{a,b\right\}$:
%
% $$A_{it}: A_{it}=\left(1-\delta\right)S_{it-1}+H_{it-1}\epsilon_{it}.$$

%% Writing the model
% The model is defined in a Yaml file: <sto6.txt sto6.yaml>.

%% Create the model object
Mu                = [1 1];
sigma             = [0.05 0;
                     0    0.05];
model = recsmodel('sto6.yaml',struct('Mu',Mu,'Sigma',sigma^2,'order',7));

%% Define approximation space
[interp,s] = recsinterpinit(15,0.73*model.sss,2*model.sss);

%% Find a first guess through the perfect foresight solution
interp = recsFirstGuess(interp,model,s,model.sss,model.xss,struct('T',5));

%% Solve for rational expectations
[interp,x] = recsSolveREE(interp,model,s);

%% Simulate the model
[~,~,~,stat] = recsSimul(model,interp,model.sss(ones(1E4,1),:),100);

%% References
%
% Miranda, M. J. and Glauber, J. W. (1995). Solving stochastic models of
% competitive storage and trade by Chebychev collocation methods. _Agricultural
% and Resource Economics Review_, 24(1), 70-77.


