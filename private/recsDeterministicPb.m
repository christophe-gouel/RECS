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
M = m+sum(nx);
T = n/(M+p+d);

e = e(ones(T,1),:);

X = reshape(X,M+p+d,T)';

x     = X(:,1:m);
xnext = [x(2:end,:); xss];
w     = X(:,(m+1):(m+nx(1)));
v     = X(:,(m+nx(1)+1):M);
z     = X(:,(M+1):(M+p));
snext = X(:,(M+p+1):(M+p+d));
s     = [s0; snext(1:end-1,:)];

%% Computation of equations and Jacobian
if nargout==2  
  %% With Jacobian
  output                  = struct('F',1,'Js',1,'Jx',1,'Jz',1,'Jsn',1,'Jxn',1);
  [f,fs,fx,fz]            = mcptransform(func,'f',s,x,z,[],[],[],params,output,ix,nx,w,v);
  [h,hs,hx,hsnext,hxnext] = mcptransform(func,'h',s,x,[],e ,snext,xnext,params,output,ix,nx,[],[]);
  [g,gs,gx]               = mcptransform(func,'g',s,x,[],e ,[],[],params,output,ix,nx,[],[]);

  identityz = eye(p);
  identityz = permute(identityz(:,:,ones(1,T)),[3 1 2]);
  identitys = eye(d);
  identitys = permute(identitys(:,:,ones(1,T)),[3 1 2]);
  Jmd = cat(3,...
            cat(2,fx,-hx,-gx),...
            cat(2,fz,identityz,zeros(T,d,p)),...
            cat(2,zeros(T,M,d),-hsnext,identitys));
  Jsub = cat(3,...
             zeros(T-1,M+p+d,M+p),...
             cat(2,fs(2:end,:,:),-hs(2:end,:,:),-gs(2:end,:,:)));
  Jsup = cat(3,...
             cat(2,zeros(T-1,M,M),-hxnext(1:end-1,:,:),zeros(T-1,d,M)),...
             zeros(T-1,M+p+d,p+d));

  J = blktridiag(permute(Jmd,[2 3 1]),permute(Jsub,[2 3 1]),permute(Jsup,[2 3 1]));

else
  %% Without Jacobian
  output = struct('F',1,'Js',0,'Jx',0,'Jz',0,'Jsn',0,'Jxn',0);
  h = mcptransform(func,'h',s,x,[],e ,snext,xnext,params,output,ix,nx,[],[]);
  f  = mcptransform(func,'f',s,x,z,[],[],[],params,output,ix,nx,w,v);
  g  = mcptransform(func,'g',s,x,[],e ,[],[],params,output,ix,nx,[],[]);

end
F = [f z-h snext-g]';
F = reshape(F,n,1);