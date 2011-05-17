% STO5 Two-country storage-trade model with supply reaction, fixed tariffs

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

disp('STO5 Two-country storage-trade model with supply reaction, fixed tariffs');
clear interp model options

% ENTER MODEL PARAMETERS
delta = 0;
r     = 0.05;
k     = 0.02;
alpha = -0.4;
theta = 0.1; % per-unit shipping cost
lambda= 1;  % relative dmd size, equals D(b)/D(a) for the same price
eta   = 1;  % relative supply size, such that the unit production cost of eta Q in b is the same as Q in a
mu    = 10;
taua = 0.06; % tariff in country a
taub = 0; % tariff in country b

% COMPUTE SHOCK DISTRIBUTION
q                 = 2; % nb of shocks
Mu                = [1 1];
Sigma             = [0.05 0;0 0.05];   % std-deviation of production shocks in a and b
[model.e,model.w] = qnwnorm([5 5],Mu,Sigma^2);
model.funrand     = @(nrep) Mu(ones(nrep,1),:)+randn(nrep,2)*Sigma;

% PACK MODEL STRUCTURE
model.func   = @sto5model;                                     % model functions
model.params = {delta,r,k,alpha,theta,lambda,eta,mu,taua,taub};         % other parameters

% DEFINE APPROXIMATION SPACE
order         = [15; 15];                                      % degree of approximation
pstat         = ((1+lambda)/(1+eta))^(1/(1/mu-alpha));   % static equilibrium price
smin          = [pstat^(1/mu)*min(model.e(:,1))*0.9; eta*pstat^(1/mu)*min(model.e(:,2))*0.9]; % if the steady-state H differs from (one,eta), adjust the 0.9
smax          = [pstat^(1/mu)*2; pstat^(1/mu)*eta*2];
interp.fspace = fundefn('spli',order,smin,smax);               % function space
snodes        = funnode(interp.fspace);                        % state collocaton nodes
s             = gridmake(snodes);
interp.Phi    = funbasx(interp.fspace);
n             = prod(order);

xinit         = [0.1*zeros(n,2) max(s(:,1).^(1/alpha),0.8) max((s(:,2)/lambda).^(1/alpha),0.8) zeros(n,2) ones(n,2)];
interp.cz     = ones(n,4);
interp.cx     = xinit;
interp.ch     = ones(n,4);

options =struct('eqsolver','ncpsolve',...
                'reesolver','SA',...
                'method','resapprox-complete',...
                'stat',1,...
                'reesolveroptions',struct('lambda',0.5),...
                'useapprox',0);

tic
[interp.cx,x,z] = recsSolveREE(interp,model,s,xinit,options);
toc
interp.cz = funfitxy(interp.fspace,interp.Phi,z);

tic
s0 = [1 1] ;
nper = 1000 ; % nb simulated periods
[ssim,xsim,esim,fsim,stat] = recsSimul(model,interp,s0(ones(100,1),:),nper,[],options) ;
toc

