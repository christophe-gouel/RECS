function [x,s,z,F] = recsSolveDeterministicPb(model,s0,T,xss,zss,sss,options)
% RECSSOLVEDETERMINISTICPB 
  
% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

if nargin <=6, options = struct([]); end

defaultopt = struct(              ...
    'eqsolver'         , 'lmmcp' ,...
    'eqsolveroptions'  , struct([]));
warning('off','catstruct:DuplicatesFound')
options = catstruct(defaultopt,options);
eqsolveroptions  = options.eqsolveroptions;

m = size(xss,2);
p = size(zss,2);
d = size(sss,2);
n = T*(m+p+d);

e = model.w'*model.e;
if isa(model.func,'char')
  func = str2func(model.func);
elseif isa(model.func,'function_handle')
  func = model.func;
else
  error('model.func must be either a string or a function handle')
end
params = model.params;

[LB,UB] = func('b',s0,[],[],[],[],[],params);
LB = [LB -inf(1,p+d)];
UB = [UB +inf(1,p+d)];
LB = reshape(LB(ones(T,1),:)',n,1);
UB = reshape(UB(ones(T,1),:)',n,1);

X = [xss zss sss];
X = X(ones(T,1),:)';
X = reshape(X,n,1);

switch lower(options.eqsolver)
 case 'fsolve'
  options = optimset('Display','off',...
                     'Jacobian','on');
  options = optimset(options,eqsolveroptions);
  [X,F,exitflag] = fsolve(@(x) recsDeterministicPb(x,func,s0,xss,p,e,params),...
                          X,...
                          options);
  if exitflag~=1, disp('No convergence'); end
 case 'lmmcp'
  [X,F,exitflag] = lmmcp(@(x) recsDeterministicPb(x,func,s0,xss,p,e,params),...
                         X,...
                         LB,...
                         UB,...
                         eqsolveroptions);
  if exitflag~=1, disp('No convergence'); end
 case 'ncpsolve'
  [X,F] = ncpsolve(@(x) ncpsolvetransform(x,@recsDeterministicPb,func,s0,xss,p,e,params),...
                   LB, ...
                   UB,...
                   X);
  F     = -F;
 case 'path'
  global par
  par   = {@recsDeterministicPb,func,s0,xss,p,e,params};
  [X,F] = pathmcp(X,...
                  LB, ...
                  UB,...
                  'pathtransform');
  clear global par
end
X = reshape(X,m+p+d,T)';
x = X(:,1:m);
z = X(:,(m+1):(m+p));
s = [s0; X(1:end-1,(m+p+1):(m+p+d))];
F = reshape(F,m+p+d,T)';
