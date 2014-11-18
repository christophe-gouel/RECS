%% STO6SP Quarterly storage model with informational subperiods and annual inelastic supply 
% This model represents the market of a storable commodity that is produced once
% a year and stored for a year-long consumption. Supply is stochastic and
% inelastic. Contrary to <sto6.html STO6>, there are informational shocks about
% the coming harvest that allows stocks to be adjusted before the full harvest
% is known.

%% Writing the model
% The model is defined in 4 Yaml files: sto6SP1.yaml, sto6SP2.yaml,
% sto6SP3.yaml, and sto6SP4.yaml.

%% Create the model object
model = recsmodelsp({'sto6SP1.yaml' 'sto6SP2.yaml' 'sto6SP3.yaml' 'sto6SP4.yaml'});
model.shocks = cell(model.nperiods,1);
model.bounds = cell(model.nperiods,2);

params = num2cell(model.params);
[k, delta, r, elastD, d] = params{:};

%% Define approximation space and shocks
clear('interp')
sigma = [eps eps eps 0.05];
n     = {50; 50; [50 5]; [50 5]};
smin  = {2.8; 2.05; [1.35 -0.15]; [0.67 -0.21]};
smax  = {6  ; 5.1 ; [3.9   0.15]; [2.72  0.21]}; 
for iperiod=1:model.nperiods
  % Shocks
  [model.shocks{iperiod}.e,model.shocks{iperiod}.w] = qnwnorm(5,0,sigma(iperiod)^2);
  model.shocks{iperiod}.funrand = @(nrep) randn(nrep,1)*sigma(iperiod);
end
% Interpolation structure
interp.fspace = cellfun(@(N,SMIN,SMAX) fundefn('spli',N,SMIN,SMAX),n,smin,smax,...
                        'UniformOutput', false);
interp.Phi    = cellfun(@(FSPACE) funbasx(FSPACE),interp.fspace,...
                        'UniformOutput', false);
interp.s      = cellfun(@(FSPACE) gridmake(funnode(FSPACE)),interp.fspace,...
                        'UniformOutput', false);

[model.ss.sss,model.ss.xss,model.ss.zss] = recsSSSP(model,{4; 3; [2 0]; [1 0]},...
                                                  {[3 1]; [2 1]; [1 1]; [0 1]});

[s1,s2,s3,s4] = interp.s{:};

%% Bounds
for iperiod=1:model.nperiods
  [LB,UB] = eval(['model.functions(iperiod).b(s' int2str(iperiod) '(1,:),params);']);
  model.bounds(iperiod,:) = {LB UB};
end

%% Provide a simple first guess
x4 = [zeros(size(s4,1),1)     (s4(:,1)/d).^(1/elastD)];
x3 = [s3(:,1)/2           ((s3(:,1)/2)/d).^(1/elastD)];
x2 = [s2(:,1)*2/3         ((s2(:,1)/3)/d).^(1/elastD)];
x1 = [s1(:,1)*3/4         ((s1(:,1)/4)/d).^(1/elastD)];

%% Solve for rational expectations
[interp,X] = recsSolveREESP(model,interp,{x1; x2; x3; x4});

%% Compare STO6 and STO6SP when informational shocks are removed
disp('Max absolute error in first subperiod storage and price (in log10)');
disp(log10(max(abs(Xcat(:,[1 5])-X{1}))));

%% Introduced information shocks
sigma = [eps 0.05/sqrt(3) 0.05/sqrt(3) 0.05/sqrt(3)];
for iperiod=1:model.nperiods
  % Shocks
  [model.shocks{iperiod}.e,model.shocks{iperiod}.w] = qnwnorm(5,0,sigma(iperiod)^2);
  model.shocks{iperiod}.funrand = @(nrep) randn(nrep,1)*sigma(iperiod);
end

[interp,X] = recsSolveREESP(model,interp,X);

[ssim,xsim,esim,stat,fsim] = recsSimulSP(model,interp,repmat(4,1000,1),200);