interpcat = interp;
clear model interp

model = recsmodelsp({'stosp1.yaml' 'stosp2.yaml' 'stosp3.yaml'});
model.shocks = cell(model.nperiods,1);

params = num2cell(model.params);
[k, delta, r, elastD, d] = params{:};

sigma = [0.02 0.03 0.06];
n     = {10; [10 5]; [10 5]};
smin  = {0.75; [0.5 -0.05];   [0.25 -0.1]};
smax  = {1.5;  [0.95 0.05];   [0.6   0.1]};
for iperiod=1:model.nperiods
  %% Shocks
  [model.shocks{iperiod}.e,model.shocks{iperiod}.w] = qnwnorm(5,0,sigma(iperiod)^2);
  model.shocks{iperiod}.funrand = @(nrep) randn(nrep,1)*sigma(1);

  %% Interpolation structure
  interp.fspace{iperiod} = fundefn('spli',n{iperiod},smin{iperiod},smax{iperiod});
  interp.Phi{iperiod}    = funbasx(interp.fspace{iperiod});
  interp.s{iperiod}      = gridmake(funnode(interp.fspace{iperiod}));
end

s3 = interp.s{3};
x3 = [zeros(size(s3,1),1) (s3(:,1)/d).^(1/elastD)];
s2 = interp.s{2};
x2 = [s2(:,1)/2 ((s2(:,1)/2)/d).^(1/elastD)];
s1 = interp.s{1};
x1 = [2*s1(:,1)/3 ((s1(:,1)/3)/d).^(1/elastD)];

fspace = fundefn('spli',10,0.25,0.6);
cx = funfitxy(fspace,Xcat(:,end),Xcat(:,[3 6]));
s = interp.s{3};
x3 = max(funeval(cx,fspace,s3(:,1)),0);

fspace = fundefn('spli',10,0.5,0.95);
cx = funfitxy(fspace,Xcat(:,end-1),Xcat(:,[2 5]));
x2 = max(funeval(cx,fspace,s2(:,1)),0);

x1 = max(funeval(interpcat.cx(:,[1 4]),interpcat.fspace,s1),0);

X = {x1 x2 x3};

[interp,X] = recsSolveREESP(model,interp,X)