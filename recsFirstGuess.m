function [interp,x,z] = recsFirstGuess(interp,model,s,sss,xss,T,options)
% RECSFIRSTGUESS 
  
% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

if nargin <=5 || isempty(T), T = 10; end
if nargin <=6, options = struct([]); end
  
[sss,xss,zss] = recsSS(model,sss,xss,options);

n = size(s,1);
x = zeros(n,size(xss,2));
z = zeros(n,size(zss,2));

for i=1:n
  [X,~,Z] = recsSolveDeterministicPb(model,s(i,:),T,xss,zss,sss,options);
  x(i,:) = X(1,:);
  z(i,:) = Z(1,:);
end

interp.cx = funfitxy(interp.fspace,interp.Phi,x);
interp.cz = funfitxy(interp.fspace,interp.Phi,z);

e = model.e;
params = model.params;
func = model.func;
K     = size(e,1);
ind   = (1:n);
ind   = ind(ones(1,K),:);
ss    = s(ind,:);
xx    = x(ind,:);
ee    = e(repmat(1:K,1,n),:);

output = struct('F',1,'Js',0,'Jx',0);
snext = func('g',ss,xx,[],ee,[],[],params,output);
xnext = funeval(interp.cx,interp.fspace,snext);

output    = struct('F',0,'Js',1,'Jx',1,'Jsn',0,'Jxn',0,'hmult',0);
[~,hs,hx] = func('h',ss,xx,[],ee,snext,xnext,params,output);

output = struct('F',1,'Js',0,'Jx',0,'Jsn',0,'Jxn',0,'hmult',0);
he     = numjac(@(E) reshape(func('h',ss,xx,[],E,snext,xnext,params,output),[],1),ee);

if abs(norm(hs(:),Inf)+norm(hx(:),Inf)+norm(he(:),Inf))<eps
  output    = struct('F',1,'Js',0,'Jx',0,'Jsn',0,'Jxn',0,'hmult',0);
  h         = func('h',[],[],[],[],s,x,params,output);
  interp.ch = funfitxy(interp.fspace,interp.Phi,h);
end

