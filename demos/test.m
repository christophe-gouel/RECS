modelcat = model;
interpcat = interp;
clear model interp

model = recsmodelsp({'test1.yaml' 'test2.yaml' 'test3.yaml'});
model.bounds = cell(model.nperiods,2);
model.shocks = cell(model.nperiods,1);

params = num2cell(model.params);
[k, delta, r, elastD, d] = params{:};

n     = 50;
smin  = {0.75; 0.49;  0.25};
smax  = {1.5;  1;     0.5};
smin  = {0.75; 0.3;  0.1};
smax  = {1.5;  1.2;     0.7};
for iperiod=1:model.nperiods
  %% Shocks
  model.shocks{iperiod}.e = 0;
  model.shocks{iperiod}.w = 1;
  
  %% Interpolation structure
  interp.fspace{iperiod} = fundefn('spli',n,smin{iperiod},smax{iperiod});
  interp.Phi{iperiod}    = funbasx(interp.fspace{iperiod});
  interp.s{iperiod}      = gridmake(funnode(interp.fspace{iperiod}));

  stmp = interp.s{iperiod};
  %% Bounds
  [LB,UB] = model.functions(iperiod).b(stmp(1,:),params);
  model.bounds(iperiod,:) = {LB UB};
end

[model.shocks{model.nperiods}.e,model.shocks{model.nperiods}.w] = qnwnorm(5,0,0.06^2);

s3 = interp.s{3};
x3 = [zeros(size(s3,1),1) (s3(:,1)/d).^(1/elastD)];
s2 = interp.s{2};
x2 = [s2(:,1)/2 ((s2(:,1)/2)/d).^(1/elastD)];
s1 = interp.s{1};
x1 = [2*s1(:,1)/3 ((s1(:,1)/3)/d).^(1/elastD)];

X = {x1 x2 x3};

[interp,X] = recsSolveREESP(model,interp,X)