function [F,J] = recsDeterministicPb(X,func,s0,xss,p,e,params,ix,nx)
% RECSDETERMINISTICPB Evaluates the equations and Jacobian of the deterministic problem.
%
% RECSDETERMINISTICPB is called by recsSolveDeterministicPb. It is not meant to be
% called directly by the user.
%
% See also RECSSOLVEDETERMINISTICPB.

% Copyright (C) 2011-2012 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
d = size(s0,2);
m = size(xss,2);

n = size(X,1);
T = n/(m+p+d+sum(nx));

e = e(ones(T,1),:);

X = reshape(X,m+p+d+sum(nx),T)';

x     = X(:,1:m);
xnext = [x(2:end,:); xss];
w     = X(:,(m+1):(m+nx(1)));
v     = X(:,(m+nx(1)+1):(m+nx(1)+nx(2)))
z     = X(:,(m+nx(1)+nx(2)+1):(m+nx(1)+nx(2)+p));
snext = X(:,(m+nx(1)+nx(2)+p+1):(m+nx(1)+nx(2)+p+d));
s     = [s0; snext(1:end-1,:)];

%% Computation of equations and Jacobian
if nargout==2  
  %% With Jacobian
  output = struct('F',1,'Js',1,'Jx',1,'Jz',1,'Jsn',1,'Jxn',1,'hmult',1);
  [f,fs,fx,fz]                    = func('f',s,x,z ,[],[]   ,[]   ,params,output);
  if nargout(func)<6
    [h,hs,hx,hsnext,hxnext]       = func('h',s,x,[],e ,snext,xnext,params,output);
  else
    [h,hs,hx,hsnext,hxnext,hmult] = func('h',s,x,[],e ,snext,xnext,params,output);
    h      = h.*hmult;
    hs     = hs.*hmult(:,:,ones(d,1));
    hx     = hx.*hmult(:,:,ones(m,1));
    hsnext = hsnext.*hmult(:,:,ones(d,1));
    hxnext = hxnext.*hmult(:,:,ones(m,1));
  end
  [g,gs,gx]               = func('g',s,x,[],e ,snext,xnext,params,output);

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
             cat(2,fs(2:end,:,:),-hs(2:end,:,:),-gs(2:end,:,:)));
  Jsup = cat(3,...
             cat(2,zeros(T-1,m,m),-hxnext(1:end-1,:,:),zeros(T-1,d,m)),...
             zeros(T-1,m+p+d,p+d));

  J = blktridiag(permute(Jmd,[2 3 1]),permute(Jsub,[2 3 1]),permute(Jsup,[2 3 1]));

else
  %% Without Jacobian
  output = struct('F',1,'Js',0,'Jx',0,'Jz',0,'Jsn',0,'Jxn',0,'hmult',1);
  f                   = func('f',s,x,z ,[],[]   ,[]   ,params,output);
  if nargout(func)<6
    h                 = func('h',s,x,[],e ,snext,xnext,params,output);
  else
    [h,~,~,~,~,hmult] = func('h',s,x,[],e ,snext,xnext,params,output);
    h                 = h.*hmult;
  end
  g                   = func('g',s,x,[],e ,snext,xnext,params,output);
end
F = [f z-h snext-g]';
F = reshape(F,n,1);