function sto2simu(situation)

%% Enter model parameters
k     = 0.02;
delta = 0.02;
r     = 0.05;
mu    = 10;
alpha = -0.4;
PF    = 1.02;
Sgbar = 0.4;

%% Define shock distribution
Mu                = 1;
sigma             = 0.05;
[model.e,model.w] = qnwnorm(5,Mu,sigma^2);
model.funrand     = @(nrep) Mu+sigma*randn(nrep,1);

%% Pack model structure
% Model functions
model.func   = @sto2model;
%%
% Model parameters
model.params = {k,delta,r,mu,alpha,PF,Sgbar};

%% Define approximation space
if situation==1
  n = 6;
else
  n = 30;
end
[interp,s] = recsinterpinit(n,min(model.e)*0.95,2);

%% Find a first guess through the perfect foresight solution
[interp,xinit] = recsFirstGuess(interp,model,s,1,[0 1 1 0],5);

%% Solve for rational expectations
[interp,x] = recsSolveREE(interp,model,s,xinit,struct('reemethod','1-step'));

if situation==1
  %% Plot storage rules
  figure
  st = (interp.fspace.a:0.01:2)';
  [~,xt] = recsSimul(model,interp,st,0,[],struct('simulmethod','solve'));
  [~,xs] = recsSimul(model,interp,st,0);
  plot(st,[xt(:,[1 4]) xs(:,[1 4])])
  ylim([0 0.55])
  leg = legend('Private stock (From equilibrium equations solve)',...
               'Public stock (From equilibrium equations solve)',...
               'Private stock (Approximated decision rule)',...
               'Public stock (Approximated decision rule)');
  set(leg,'Location','NorthWest')
  set(leg,'Box','off')                    
  xlabel('Availability')
  ylabel('Stock')
  title('Storage rules with a grid of 6 points')
else
  %% Simulate the model
  stream0 = RandStream('mt19937ar','Seed',0);
  RandStream.setGlobalStream(stream0);
  disp('Long-run statistics from equilibrium equations solve')
  disp('')
  stream0.reset
  [~,~,~,~,stat]  = recsSimul(model,interp,ones(1E3,1),118,[],...
                              struct('simulmethod','solve'));
  disp('Long-run statistics if simulated with approximated decision rules')
  disp('')
  stream0.reset
  [~,~,~,~,stat2] = recsSimul(model,interp,ones(1E3,1),118);
  close all
end