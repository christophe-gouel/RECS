model = recsmodelsp({'sto1sp.yaml'});

params = num2cell(model.params);
[k, delta, r, h, mu, elastD, elastS] = params{:};

sigma = {0.05};
n     = {40};
smin  = {0.7};
smax  = {1.5};
for iperiod=1:model.nperiods
  %% Shocks
  [model.shocks(iperiod).e,model.shocks(iperiod).w] = qnwnorm(5,1,sigma{iperiod});
  model.shocks(iperiod).funrand = @(nrep) randn(nrep,1)*sigma{1};

  %% Interpolation structure
  interp.fspace{iperiod} = fundefn('spli',n{iperiod},smin{iperiod},smax{iperiod});
  interp.Phi{iperiod}    = funbasx(interp.fspace{iperiod});
  interp.s{iperiod}      = gridmake(funnode(interp.fspace{iperiod}));
end

s = interp.s{1};
x = [zeros(size(s,1),1) ones(size(s,1),1) s.^(1/elastD)];

X = {x};

[interp,X] = recsSolveREESP(model,interp,X)