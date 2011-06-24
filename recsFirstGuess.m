function [interp,x,z] = recsFirstGuess(interp,model,s,sss,xss,T,options)
% RECSFIRSTGUESS

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

if nargin <=5 || isempty(T), T = 10; end
if nargin <=6, options = struct([]); end

if isa(model.func,'char')
  model.func = str2func(model.func);
elseif isa(model.func,'function_handle')
  model.func = model.func;
else
  error('model.func must be either a string or a function handle')
end

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

params = model.params;
func   = model.func;
if ~isfield(interp,'ch')
  e      = model.e;
  K      = size(e,1);
  ind    = (1:n);
  ind    = ind(ones(1,K),:);
  ss     = s(ind,:);
  xx     = x(ind,:);
  ee     = e(repmat(1:K,1,n),:);

  output = struct('F',1,'Js',0,'Jx',0);
  snext = func('g',ss,xx,[],ee,[],[],params,output);
  xnext = funeval(interp.cx,interp.fspace,snext);

  output    = struct('F',0,'Js',1,'Jx',1,'Jsn',0,'Jxn',0,'hmult',0);
  [~,hs,hx] = func('h',ss,xx,[],ee,snext,xnext,params,output);

  output = struct('F',1,'Js',0,'Jx',0,'Jsn',0,'Jxn',0,'hmult',0);
  if size(ss,1)>100;
    i = unique(randi(size(ss,1),100,1));
  else
    i = 1:size(ss,1);
  end
  he     = numjac(@(E) reshape(func('h',ss(i,:),xx(i,:),[],E,snext(i,:),xnext(i,:),...
                                    params,output),[],1),ee(i,:));
end

if isfield(interp,'ch') || abs(norm(hs(:),Inf)+norm(hx(:),Inf)+norm(he(:),Inf))<eps
  output    = struct('F',1,'Js',0,'Jx',0,'Jsn',0,'Jxn',0,'hmult',0);
  h         = func('h',[],[],[],[],s,x,params,output);
  interp.ch = funfitxy(interp.fspace,interp.Phi,h);
end

