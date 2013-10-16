%% Quarterly storage model with annual inelastic supply
% This model represents the market of a storable commodity that is produced once
% a year and stored for a year-long consumption. Supply is stochastic and
% inelastic. Except for its distribution, no information about the coming
% harvest is known.

%% Model's structure
%

%% Writing the model
% The model is defined in a Yaml file: <sto7.txt sto7.yaml>.

%% Create the model object
model = recsmodel('sto7.yaml',struct('Mu',0,'Sigma',0.05^2,'order',5));

%% Define approximation space
interp = recsinterpinit(50,model.sss*0.7,model.sss*1.5);

%% Find a first guess through the perfect foresight solution
interp = recsFirstGuess(interp,model);

%% Solve for rational expectations
[interp,Xcat] = recsSolveREE(interp,model);

%% Plot the decision rules
recsDecisionRules(model,interp,[],[],[],struct('simulmethod','solve'));
for i=1:model.dim{2}
  subplot(3,4,i)
  xlabel(model.symbols.states{1});
  ylabel(model.symbols.controls{i});
end

%% Simulate the model
[ssim,~,~,stat] = recsSimul(model,interp,model.sss(ones(1000,1),:),200,[],...
                            struct('accuracy',1));
subplot(3,4,1)
xlabel(model.symbols.states{1});
for i=1:model.dim{2}
  subplot(3,4,i+1)
  xlabel(model.symbols.controls{i});
end

