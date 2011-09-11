% CS1 Consumption/Saving model with borrowing constraint

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

disp('CS1 Consumption/Saving model with borrowing constraint');
clear interp model options

%% Pack model structure
Mu                = 100;
sigma             = 10;
[model.e,model.w] = qnwnorm(5,Mu,sigma^2);
model.funrand     = @(nrep) Mu+sigma*randn(nrep,1);
model.func        = @cs1model;
model.params      = cs1model('params');

% If your RECS installation is complete, you can also generate the Matlab model
% file and pack the model structure with the following command
%{
model = recsmodelinit('cs1.yaml',...
                      struct('Mu',100,'Sigma',100,'order',5));
%}

%% Find deterministic steady-state
disp('Deterministic steady-state')
[sss,xss,zss] = recsSS(model,100,100)

%% Define approximation space
interp.fspace = fundefn('spli',20,sss*0.5,sss*2);           % function space
s             = gridmake(funnode(interp.fspace));                    % state collocaton nodes
interp.Phi    = funbasx(interp.fspace);

%% First-guess: Consumption equal to cash on hand
x         = s;
interp.cx = funfitxy(interp.fspace,interp.Phi,x);
interp.ch = [];  % To force the solve to compute the approximation of the
                 % expectations function

%% Solve for rational expectations
tic
[interp,x] = recsSolveREE(interp,model,s,x);
toc

figure
plot(s,x)
xlabel('Cash on hand')
ylabel('Consumption')

%% Simulate the model
[~,~,~,~,stat] = recsSimul(model,interp,sss(ones(1000,1),:),200);
