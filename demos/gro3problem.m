function [X,f,exitflag] = gro3problem(X,z)
% GRO3PROBLEM Solves the model GRO3 for different values of risk aversion and IES

[model,interp,s,x,options] = X{:};

tau = z(1);
Psi = z(2);
nu  = (1-tau)/(1-1/Psi);
model.params(1)         = tau;
model.params(end-1:end) = [nu Psi];

[interp,x,~,f,exitflag] = recsSolveREE(interp,model,s,x,options);

X = {model interp s x options};