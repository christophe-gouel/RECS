function [F,J] = DeterministicProblem(X,func,s0,xss,p,e,params)

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

d = size(s0,2);
m = size(xss,2);

n = size(X,1);
T = n/(m+p+d);

e = e(ones(T,1),:);

X = reshape(X,m+p+d,T)';

x     = X(:,1:m);
xnext = [x(2:end,:); xss];
z     = X(:,(m+1):(m+p));
snext = X(:,(m+p+1):(m+p+d));
s     = [s0; snext(1:end-1,:)];

if nargout==2
  output = struct('F',1,'Js',1,'Jx',1,'Jz',1,'Jsn',1,'Jxn',1);
  [f,fs,fx,fz] = feval(func,'f',s,x,z,[],[],[],params,output);
  [h,hs,hx,hsnext,hxnext] = feval(func,'h',s,x,[],e,snext,xnext,params,output);
  [g,gs,gx] = feval(func,'g',s,x,[],e,snext,xnext,params,output);

  identityz = eye(p);
  identityz = permute(identityz(:,:,ones(1,T)),[3 1 2]);
  identitys = eye(d);
  identitys = permute(identitys(:,:,ones(1,T)),[3 1 2]);
  Jmd = cat(3,...
            cat(2,fx,-hx,-gx),...
            cat(2,fz,identityz,zeros(T,d,p)),...
            cat(2,zeros(T,m,d),-hsnext,identitys));
  Jsub = cat(3,...
             zeros(T-1,m+p+d,m+p),...
             cat(2,fs(2:end,:,:),-hs(2:end,:,:),gs(2:end,:,:)));
  Jsup = cat(3,...
             cat(2,zeros(T-1,m,m),-hxnext(1:end-1,:,:),zeros(T-1,d,m)),...
             zeros(T-1,m+p+d,p+d));

  J = blktridiag(permute(Jmd,[2 3 1]),permute(Jsub,[2 3 1]),permute(Jsup,[2 3 1]));

  else
  output = struct('F',1,'Js',0,'Jx',0,'Jz',0,'Jsn',0,'Jxn',0);
  f = feval(func,'f',s,x,z,[],[],[],params,output);
  h = feval(func,'h',s,x,[],e,snext,xnext,params,output);
  g = feval(func,'g',s,x,[],e,snext,xnext,params,output);
end
F = [f z-h snext-g]';
F = reshape(F,n,1);