%%
model = recsmodel('stospcat.yaml',struct('Mu',0,'Sigma',0.06^2,'order',5));

%% Define approximation space
[interp,s] = recsinterpinit(50,model.sss*0.75,model.sss*1.5);

%% Find a first guess through the perfect foresight solution
interp = recsFirstGuess(interp,model);

%% Solve for rational expectations
[interp,Xcat] = recsSolveREE(interp,model);

%% Plot the decision rules
recsDecisionRules(model,interp,[],[],[],struct('simulmethod','solve'));

%% Simulate the model
[ssim,~,~,stat] = recsSimul(model,interp,model.sss(ones(1000,1),:),200);

%% Assess solution accuracy
recsAccuracy(model,interp,ssim);
