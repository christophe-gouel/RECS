function [c,x,z,f,exitflag] = recsSolveREEFull(interp,model,s,x,options)

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

if nargin <=4, options = struct([]); end

defaultopt = struct(                  ...
    'eqsolver'          , 'lmmcp'    ,...
    'eqsolveroptions'   , struct([]) ,...
    'functional'        , 0          ,...
    'method'            , 'resapprox-complete');
warning('off','catstruct:DuplicatesFound')

options = catstruct(defaultopt,options);
eqsolver         = lower(options.eqsolver);
eqsolveroptions  = options.eqsolveroptions;
functional         = options.functional;
method           = lower(options.method);

e      = model.e;
params = model.params;
w      = model.w;
if isa(model.func,'char')
  func = str2func(model.func);
elseif isa(model.func,'function_handle')
  func = model.func;
else
  error('model.func must be either a string or a function handle')
end

switch method
 case 'expfunapprox'
  c      = interp.ch;
 case 'resapprox-complete'
  c      = interp.cx;
 otherwise
end
fspace = interp.fspace;
Phi    = interp.Phi;
if functional, params = [params fspace c]; end

[n,m]  = size(x);
output = struct('F',1,'Js',0,'Jx',0,'Jsn',0,'Jxn',0,'hmult',0);
p      = size(func('h',s(1,:),x(1,:),[],e(1,:),s(1,:),x(1,:),params,output),2);
K      = length(w);               % number of shock values
z      = zeros(n,0);

[~,grid] = spblkdiag(zeros(m,m,n),[],0);
X        = [reshape(x',[n*m 1]); reshape(c',[],1)];
[LB,UB]  = func('b',s,[],[],[],[],[],params);
LB       = [reshape(LB',[n*m 1]); -inf(n*size(c,2),1)];
UB       = [reshape(UB',[n*m 1]); +inf(n*size(c,2),1)];

% $$$ [f,J] = recsFullPb(X,s,func,params,grid,e,w,fspace,method,Phi,m,functional);
% $$$ Jnum = numjac(@(VAR) recsFullPb(VAR,s,func,params,grid,e,w,fspace,method,Phi,m,functional),X);
% $$$ spy(J)
% $$$ figure
% $$$ spy(Jnum)
% $$$ % $$$
% $$$ norm(full(J)-Jnum)
% $$$ norm(full(J(1:n*m,n*m+1:n*(m+size(c,2))))-Jnum(1:n*m,n*m+1:n*(m+size(c,2))))
% $$$ norm(full(J(n*m+1:n*(m+size(c,2)),n*m+1:n*(m+size(c,2))))-Jnum(n*m+1:n*(m+size(c,2)),n*m+1:n*(m+size(c,2))))
% $$$ z = [];
% $$$ return
exitflag = 1;

switch eqsolver
 case 'fsolve'
  options = optimset('Display','off',...
                     'Jacobian','on');
  options = optimset(options,eqsolveroptions);
  [X,f,exitflag] = fsolve(@(VAR) recsFullPb(VAR,s,func,params,grid,e,w,...
                                            fspace,method,Phi,m,functional),...
                          X,options);
  if exitflag~=1, disp('No convergence'); end
 case 'lmmcp'
  [X,f,exitflag] = lmmcp(@(VAR) recsFullPb(VAR,s,func,params,grid,e,w,...
                                           fspace,method,Phi,m,functional),...
                         X,LB,UB,eqsolveroptions);
  if exitflag~=1, disp('No convergence'); end
 case 'ncpsolve'
  [X,f] = ncpsolve(@(VAR) ncpsolvetransform(VAR,@recsFullPb,s,func,params,grid,e,w,...
                                            fspace,method,Phi,m,functional),...
                   LB,UB,X);
  f     = -f;
 case 'path'
  global par
  par   = {@recsFullPb,s,func,params,grid,e,w,fspace,method,Phi,m,functional};
  [X,f] = pathmcp(X,LB,UB,'pathtransform');
  clear global par
end

if exitflag~=1
  warning('recs:FailureREE','Failure to find a rational expectations equilibrium');
end

x     = reshape(X(1:n*m),m,n)';
c     = reshape(X(n*m+1:end),[],n)';

% Calculation of z
ind    = (1:n);
ind    = ind(ones(1,K),:);
ss     = s(ind,:);
xx     = x(ind,:);
output = struct('F',1,'Js',0,'Jx',0);
snext  = func('g',ss,xx,[],e(repmat(1:K,1,n),:),[],[],params,output);
switch method
 case 'expfunapprox'
  h   = funeval(c,fspace,snext);
  if nargout(func)==6
    output            = struct('F',0,'Js',0,'Jx',0,'Jsn',0,'Jxn',0,'hmult',1);
    [~,~,~,~,~,hmult] = func('h',[],[],[],e(repmat(1:K,1,n),:),snext,zeros(size(snext,1),m),params,output);
    h                 = h.*hmult;
  end

   case 'resapprox-complete'
    [LB,UB] = func('b',snext,[],[],[],[],[],params);
    xnext   = min(max(funeval(c,fspace,snext),LB),UB);
    output  = struct('F',1,'Js',0,'Jx',0,'Jsn',0,'Jxn',0,'hmult',1);
    if nargout(func)<6
      h                = func('h',ss,xx,[],e(repmat(1:K,1,n),:),snext,xnext,params,output);
    else
      [h,~,~,~,~,hmult] = func('h',ss,xx,[],e(repmat(1:K,1,n),:),snext,xnext,params,output);
      h               = h.*hmult;
    end
end
z     = reshape(w'*reshape(h,K,n*p),n,p);

