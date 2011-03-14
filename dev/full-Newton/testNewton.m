% DEMrecs01 Competitive storage without supply reaction

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

clear interp model options

% ENTER MODEL PARAMETERS
alpha = -0.4;
k     = 0.02;
delta = 0;
r     = 0.05;

% COMPUTE SHOCK DISTRIBUTION
[model.e,model.w] = qnwnorm(5,1,0.05^2);

% PACK MODEL STRUCTURE
model.func   = 'mdsem01';                               % model functions
model.params = {alpha,k,delta,r};               % other parameters
model.dim    = {1,1,1};

% DEFINE APPROXIMATION SPACE
order         = 200;                                          % degree of approximation
smin          = min(model.e);
smax          = 2;
interp.fspace = fundefn('spli',order,smin,smax);                 % function space
snodes        = funnode(interp.fspace);                             % state collocaton nodes
interp.s      = gridmake(snodes);
interp.Phi    = funbas(interp.fspace);
n             = order;
s             = interp.s;

xinit         = zeros(n,1);
interp.cz     = ones(n,1);
interp.cx     = zeros(n,1);
optset('ncpsolve','type','smooth'); % 'minmax' / 'smooth' / 'spminmax'

m = 1;
p = 1;

global par
par = {'recsCompleteModel',interp,model};
X = [reshape(xinit',[n*m 1])
     reshape(interp.cz',[n*p 1])
     reshape(interp.cz',[n*p 1])];
[LB,UB] = feval(model.func,'b',s,[],[],[],[],[],model.params{:});
LB = [reshape(LB',[n*m 1])
      -inf(2*n*p,1)];
UB = [reshape(UB',[n*m 1])
      +inf(2*n*p,1)];

% $$$ tic
% $$$ [x,f] = pathmcp(X,...
% $$$                 LB, ...
% $$$                 UB,...
% $$$                 'pathtransform');
% $$$ toc
% $$$ return
tic
[X,f] = ncpsolve('ncpsolvetransform',...
                 LB, ...
                 UB,...
                 X,...
                 'recsCompleteModel',...
                 interp,...
                 model);
toc
% $$$ clear global par

x = reshape(X(1:n*m),m,n)';
z = reshape(X(n*m+1:n*(m+p)),p,n)';
interp.ch = reshape(X(n*(m+p)+1:n*(m+2*p)),p,n)';

interp.cx = interp.Phi\x;

return
tic
[interp.cz,x] = recsSolveREE(interp,model,[],s,xinit);
toc
interp.cx = funfitxy(interp.fspace,interp.Phi,x);

options.method  = 'expfunapprox';
interp.ch       = s.^(1/alpha);
interp.cz       = ones(n,1);
interp.cx       = zeros(n,1);
tic
[interp.ch,x,z] = recsSolveREE(interp,model,options,s,xinit);
toc
interp.cx = funfitxy(interp.fspace,interp.Phi,x);
interp.cz = funfitxy(interp.fspace,interp.Phi,z);

options.method    = 'resapprox-simple';
interp.cz     = ones(n,1);
interp.cx     = zeros(n,1);
tic
[interp.cx,x,z] = recsSolveREE(interp,model,options,s,xinit);
toc
interp.cz = funfitxy(interp.fspace,interp.Phi,z);

options.method    = 'resapprox-complete';
interp.cz     = ones(n,1);
interp.cx     = zeros(n,1);
tic
[interp.cx,x,z] = recsSolveREE(interp,model,options,s,xinit);
toc
interp.cz = funfitxy(interp.fspace,interp.Phi,z);

