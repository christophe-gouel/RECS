function [G,J] = recsFullPb(X,s,func,params,grid,e,w,fspace,method,Phi,m,functional)
% RECSFULLPB

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

n     = size(s,1);
x     = X(1:m*n);
x     = reshape(x,m,n)';
c     = reshape(X(m*n+1:end),[],n)';
if functional, params{end} = c; end

switch method
 case 'expfunapprox'
  p  = size(c,2);
  if nargout==1
    output = struct('F',1,'Js',0,'Jx',0,'Jsn',0,'Jxn',0,'hmult',0);
    R = funeval(c,fspace,Phi)-func('h',[],[],[],[],s,x,params,output);
    R = reshape(R',n*p,1);
  end
 case 'resapprox-complete'
  R = funeval(c,fspace,Phi)-x;
  R = reshape(R',n*m,1);
end

if nargout==2 % With Jacobian
  [F,Fx,Fc] = recsEquilibrium(reshape(x',[n*m 1]),...
                              s,zeros(n,0),func,params,grid,c,e,w,fspace,method);

  B = funbas(fspace,s);

  switch method
   case 'expfunapprox'
    output = struct('F',1,'Js',0,'Jx',0,'Jsn',0,'Jxn',1,'hmult',0);
    [h,~,~,~,hxnext] = func('h',[],[],[],[],s,x,params,output);
    R  = funeval(c,fspace,Phi)-h;
    R  = reshape(R',n*p,1);
    Rx = -spblkdiag(permute(hxnext,[2 3 1]));

    Rc = kron(B,speye(p));
   
   case 'resapprox-complete'
    Rx = -speye(n*m);

    Rc = kron(B,speye(m));

  end

  J = [Fx Fc;
       Rx Rc];
else % Without Jacobian
  F = recsEquilibrium(reshape(x',[n*m 1]),...
                       s,zeros(n,0),func,params,grid,c,e,w,fspace,method);
end

G = [F; R];
