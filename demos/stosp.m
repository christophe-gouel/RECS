interpcat = interp;
clear model interp

model = recsmodelsp({'stosp1.yaml' 'stosp2.yaml' 'stosp3.yaml' 'stosp4.yaml'});
model.shocks = cell(model.nperiods,1);
model.bounds = cell(model.nperiods,2);

params = num2cell(model.params);
[k, delta, r, elastD, d] = params{:};

sigma = [0.02 0.03 0.06];
sigma = [eps eps eps 0.05];
n     = {50; 50; [50 5]; [50 5]};
% smin  = {0.75; [0.5 -0.05];   [0.25 -0.1]};
% smax  = {1.5;  [0.95 0.05];   [0.6   0.1]};
smin  = {2.8; 2.05; [1.35 -0.05]; [0.67 -0.1]};
smax  = {6  ; 5.1 ; [3.9   0.05]; [2.72  0.1]}; 
for iperiod=1:model.nperiods
  %% Shocks
  [model.shocks{iperiod}.e,model.shocks{iperiod}.w] = qnwnorm(5,0,sigma(iperiod)^2);
  model.shocks{iperiod}.funrand = @(nrep) randn(nrep,1)*sigma(1);

  %% Interpolation structure
  interp.fspace{iperiod} = fundefn('spli',n{iperiod},smin{iperiod},smax{iperiod});
  interp.Phi{iperiod}    = funbasx(interp.fspace{iperiod});
  interp.s{iperiod}      = gridmake(funnode(interp.fspace{iperiod}));
end

[s1,s2,s3,s4] = interp.s{:};

for iperiod=1:model.nperiods
  %% Bounds
  [LB,UB] = eval(['model.functions(iperiod).b(s' int2str(iperiod) '(1,:),params);']);
  model.bounds(iperiod,:) = {LB UB};
end

x4 = [zeros(size(s4,1),1)     (s4(:,1)/d).^(1/elastD)];
x3 = [s3(:,1)/2           ((s3(:,1)/2)/d).^(1/elastD)];
x2 = [s2(:,1)*2/3         ((s2(:,1)/3)/d).^(1/elastD)];
x1 = [s1(:,1)*3/4         ((s1(:,1)/4)/d).^(1/elastD)];

X = {x1 x2 x3 x4};

[interp,X] = recsSolveREESP(model,interp,X)