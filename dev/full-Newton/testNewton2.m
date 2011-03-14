% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

clear interp model options

% ENTER MODEL PARAMETERS
delta = 0;
r     = 0.05;
k     = 0.02;
alpha = -0.4;
tau   = 0.2;
rho   = 0.6;
sigma = 0.16;

% COMPUTE SHOCK DISTRIBUTION
Mu                = [1 0];
Sigma             = [0.05 0;0 sigma];
[model.e,model.w] = qnwnorm([5 5],Mu,Sigma^2);
model.funrand     = @(nrep) normrnd(Mu(ones(nrep,1),:),repmat(diag(Sigma)',nrep,1));
model.dim         = {2,4,1};

% PACK MODEL STRUCTURE
model.func   = 'mdsem06';                                     % model functions
model.params = {delta,r,k,alpha,tau,rho,sigma};               % other parameters

% DEFINE APPROXIMATION SPACE
order         = [35; 15];                                      % degree of approximation
smin          = [min(model.e(:,1))*0.95; 0.4];
smax          = [1.6; 2];
interp.fspace = fundefn('spli',order,smin,smax);               % function space
snodes        = funnode(interp.fspace);                        % state collocaton nodes
interp.s             = gridmake(snodes);
interp.Phi    = funbas(interp.fspace);
n             = prod(order);
s = interp.s;

xinit         = [zeros(n,1) max(min(s(:,1).^(1/alpha),s(:,2)+tau),s(:,2)-tau) zeros(n,2)];
interp.cz     = max(min(ones(n,1),s(:,2).^rho+tau),s(:,2).^rho-tau);
interp.cx     = xinit;
interp.ch     = max(min(s(:,1).^(1/alpha),s(:,2)+tau),s(:,2)-tau);

% This problem is difficult to solve. For some parameters values, it requires:
%  - Path as mcp-solver
%  - Approximation of the expectation function

options.eqsolver   = 'path';
options.reesolver   = 'SA';
options.simulmethod = 'solve';
optset('ncpsolve','type','spminmax'); % 'minmax' / 'smooth' / 'spminmax'

[d,m,p] = model.dim{:};
X = [reshape(xinit',[n*m 1])
     reshape(interp.cz',[n*p 1])
     reshape(interp.cz',[n*p 1])];
[LB,UB] = feval(model.func,'b',s,[],[],[],[],[],model.params{:});
LB = [reshape(LB',[n*m 1])
      -inf(2*n*p,1)];
UB = [reshape(UB',[n*m 1])
      +inf(2*n*p,1)];

global par
par = {'recsCompleteModel',interp,model};
tic
[X,f] = pathmcp(X,...
                LB, ...
                UB,...
                'pathtransform');
toc

x = reshape(X(1:n*m),m,n)';
z = reshape(X(n*m+1:n*(m+p)),p,n)';
interp.ch = reshape(X(n*(m+p)+1:n*(m+2*p)),p,n)';
interp.cx = interp.Phi\x;
interp.cz = interp.Phi\z;

order         = [15; 15];                                      % degree of approximation
interp1.fspace = fundefn('spli',order,smin,smax);               % function space
snodes        = funnode(interp1.fspace);                        % state collocaton nodes
interp1.s             = gridmake(snodes);
interp1.Phi    = funbas(interp1.fspace);
n             = prod(order);
s = interp1.s;

interp1.cx = funconv(interp.cx,interp.fspace,interp1.fspace);
interp1.cz = funconv(interp.cz,interp.fspace,interp1.fspace);
interp1.ch = funconv(interp.ch,interp.fspace,interp1.fspace);
X = [reshape(funeval(interp1.cx,interp1.fspace,s)',[n*m 1])
     reshape(funeval(interp1.cz,interp1.fspace,s)',[n*p 1])
     reshape(interp1.ch',[n*p 1])];
[LB,UB] = feval(model.func,'b',s,[],[],[],[],[],model.params{:});
LB = [reshape(LB',[n*m 1])
      -inf(2*n*p,1)];
UB = [reshape(UB',[n*m 1])
      +inf(2*n*p,1)];
par = {'recsCompleteModel',interp1,model};

% $$$ tic
% $$$ [X,f] = ncpsolve('ncpsolvetransform',...
% $$$                  LB, ...
% $$$                  UB,...
% $$$                  X,...
% $$$                  'recsCompleteModel',...
% $$$                  interp1,...
% $$$                  model);
% $$$ toc
% $$$ return
clear mcppath
return
tic
[X,f] = pathmcp(X,...
                LB, ...
                UB,...
                'pathtransform');
toc
return

options.method      = 'expfunapprox';
tic
[interp.ch,x,z] = recsSolveREE(interp,model,options,s,xinit);
toc
interp.cx = funfitxy(interp.fspace,interp.Phi,x);
interp.cz = funfitxy(interp.fspace,interp.Phi,z);

tic
[ssim,xsim,esim] = recsSimul(model,interp,[1 1],200,options);
toc
