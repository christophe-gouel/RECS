Mu                = 1;
sigma             = 0.05;

model = recsmodelinit('sto1FX.yaml');
[model.e,model.w] = qnwnorm(7,Mu,sigma^2);
model.funrand     = @(nrep) Mu(ones(nrep,1),:)+randn(nrep,1)*sigma;

[interp,s] = recsinterpinit(40,0.7,1.5);

x = [zeros(size(s)) max(0.7,s.^(1/model.params(4)))];

options = struct('explicit' , 1,...
                 'useapprox', 0);

[interp,x] = recsSolveREE(interp,model,s,x,options);

[~,~,~,~] = recsSimul(model,interp,ones(1000,1),200);