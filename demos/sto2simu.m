function sto2simu(situation)

%% Pack model structure
model = recsmodelinit('sto2.yaml',struct('Mu',1,'Sigma',0.05^2,'order',7));

%% Define approximation space
if situation==1
  n = 6;
else
  n = 30;
end
[interp,s] = recsinterpinit(n,min(model.e)*0.95,2);

%% Find a first guess through the perfect foresight solution
[interp,xinit] = recsFirstGuess(interp,model,s,model.sss,model.xss,5);

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