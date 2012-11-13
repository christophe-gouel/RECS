%% STO5 Storage-trade model of a small country

%% Model's structure
%
% *Response variables* Storage ($S$), Price ($P$), Import ($M$) and Export
% ($X$).
%
% *State variable* Availability ($A$) and World price ($P^{w}$).
%
% *Shock* Production shocks ($\epsilon$) and Innovation on world price ($\nu$).
%
% *Parameters* Unit physical storage cost ($k$), Interest rate ($r$), Demand
% price elasticity ($\alpha$), World price autocorrelation ($\rho$) and Trade
% cost ($\theta$).
%
% *Equilibrium equations*
%
% $$S_{t}: S_{t}\ge 0 \quad \perp \quad P_{t}+k-\frac{\mathrm{E}_{t}\left(P_{t+1}\right)}{1+r}\ge 0,$$
%
% $$P_{t}: A_{t}+M_t={P_{t}}^{\alpha}+S_{t}+X_t,$$
%
% $$M_{t}: M_{t}\ge 0 \quad \perp \quad P^w_{t}+\theta\ge P_{t},$$
%
% $$X_{t}: X_{t}\ge 0 \quad \perp \quad P_{t}\ge P^w_{t}-\theta.$$
%
% *Transition equation*
%
% $$A_{t}: A_{t}=S_{t-1}+\epsilon_{t},$$
%
% $$P^w_{t}: \ln P^w_{t} = \rho\ln P^w_{t-1}+\nu_{t}.$$

%% Writing the model
% The model is defined in a Yaml file: <sto5.txt sto5.yaml>.

%% Pack model structure
Mu                = [1 0];
sigma             = [0.05 0;
                     0    1];
model = recsmodelinit('sto5.yaml',struct('Mu',Mu,'Sigma',sigma^2,'order',7));

%% Define approximation space
% Minimum and maximum values of the state variable grid
smin          = [min(model.e(:,1))*0.95; 0.4 ];
smax          = [1.6;                    2.12];
%%
% Interpolation structure
[interp,s]    = recsinterpinit(15,smin,smax);

%% Find a first guess through the perfect foresight solution
% We use the solver |ncpsolve| because |lmmcp| does not work in this case:
options = struct('eqsolver','ncpsolve');
%%
interp = recsFirstGuess(interp,model,s,model.sss,model.xss,5,options);


%% Solve for rational expectations
% Define options
options = struct('funapprox'  ,'expfunapprox',...
                 'simulmethod','solve'       ,...
                 'stat'       ,1);
%%
[interp,x] = recsSolveREE(interp,model,s,[],options);

%% Simulate the model
recsSimul(model,interp,repmat([1 1],10,1),200,[],options);
subplot(2,3,1)
xlabel('Availability')
ylabel('Frequency')
subplot(2,3,2)
xlabel('World price')
ylabel('Frequency')
subplot(2,3,3)
xlabel('Storage')
ylabel('Frequency')
subplot(2,3,4)
xlabel('Price')
ylabel('Frequency')
subplot(2,3,5)
xlabel('Import')
ylabel('Frequency')
subplot(2,3,6)
xlabel('Export')
ylabel('Frequency')
