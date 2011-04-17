% DEMSTO3 Competitive storage with price-band backed by public storage

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

disp('DEMSTO3 Competitive storage with price-band backed by public storage');
clear interp model options

% ENTER MODEL PARAMETERS
delta = 0;
r     = 0.05;
alpha = -0.4;
k     = 0.02;
mu    = 10;
PF    = 0.9;
PC    = 1.1;
Sgbar = 0.4;

% COMPUTE SHOCK DISTRIBUTION
Mu                = 1;
sigma             = 0.05;
[model.e,model.w] = qnwnorm(5,Mu,sigma^2);
model.funrand     = @(nrep) Mu+sigma*randn(nrep,1);

% PACK MODEL STRUCTURE
model.func   = @msto3;                                      % model functions
model.params = {delta,r,alpha,k,mu,PF,PC,Sgbar};               % other parameters

% DEFINE APPROXIMATION SPACE
order         = [20; 20];                                      % degree of approximation
smin          = [min(model.e(:,1))*0.95; 0];
smax          = [1.4; Sgbar];
interp.fspace = fundefn('spli',order,smin,smax);               % function space
snodes        = funnode(interp.fspace);                        % state collocaton nodes
s             = gridmake(snodes);
interp.Phi    = funbasx(interp.fspace);
n             = prod(order);

xinit         = [zeros(n,1) ones(n,2) zeros(n,2)];
interp.cz     = ones(n,2);
interp.cx     = xinit;

options = struct('eqsolver','lmmcp',...
                 'simulmethod','solve',...
                 'reesolver','SA',...
                 'method','resapprox-complete');

interp.cz = ones(n,2);
interp.cx = xinit;
tic
[interp.cx,x,z] = recsSolveREE(interp,model,s,xinit,options);
toc
interp.cz = funfitxy(interp.fspace,interp.Phi,z);
return
if exist('mcppath')
  options.eqsolver = 'path';
  interp.cz = ones(n,2);
  interp.cx = xinit;
  tic
    [interp.cx,x,z] = recsSolveREE(interp,model,s,xinit,options);
  toc
  interp.cz = funfitxy(interp.fspace,interp.Phi,z);
end

tic
[ssim,xsim,esim] = recsSimul(model,interp,[1 0],200,[],options);
toc
