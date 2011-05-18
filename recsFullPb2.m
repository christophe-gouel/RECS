function [F,J] = recsFullPb2(x,s,func,params,grid,c,e,w,fspace,Phi)
% RECSFULLPB 
%
% This code is not optimized. Some formulas are calculated both in this file and
% in recsEquilibrium. It could be improved.
  
% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt
  
[n,d] = size(s);
x     = reshape(x,[],n)';
m     = size(x,2);

F2 = funeval(c,fspace,Phi)-x;
F2 = reshape(F2',n*m,1);

if nargout==2 % With Jacobian
  [F1,J11] = recsEquilibrium(x,s,zeros(n,0),func,params,grid,c,e,w,fspace,'resapprox-complete');  

  B = funbas(fspace,s);

  J21 = -speye(n*m);
  J22 = kron(speye(m),B);
  
  % Calculation of J12
  J12   = sparse(n*m,n*m);
  k     = size(e,1);
  ind   = (1:n);
  ind   = ind(ones(1,k),:);
  ss    = s(ind,:);
  xx    = x(ind,:);
  ee    = e(repmat(1:k,1,n),:);

  output   = struct('F',1,'Js',0,'Jx',0);
  snext    = func('g',ss,xx,[],ee,[],[],params,output);
  Bsnext   = funbas(fspace,snext);

  [LB,UB]              = func('b',snext,[],[],[],[],[],params);
%  xnext                = min(max(funeval(c,fspace,Bsnext),LB),UB);
  xnext                = min(max(funeval(c,fspace,snext),LB),UB);

  output = struct('F',1,'Js',0,'Jx',0,'Jsn',0,'Jxn',1,'hmult',1);
  if nargout(func)<6
    [h,~,~,~,hxnext]       = func('h',ss,xx,[],ee,snext,xnext,params,output);
  else
    [h,~,~,~,hxnext,hmult] = func('h',ss,xx,[],ee,snext,xnext,params,output);
    h      = h.*hmult;
    hxnext = hxnext.*hmult(:,:,ones(m,1));
  end
  size(hxnext)
  size(h)
  p           = size(h,2);
  z           = reshape(w'*reshape(h,k,n*p),n,p);
  output      = struct('F',1,'Js',0,'Jx',1,'Jz',1);
  [~,~,~,fz] = func('f',s,x,z,[],[],[],params,output);

  Jtmp = reshape(w'*reshape(hxnext,k,n*p*m),n,p,m);
  Jtmp = arraymult(fz,repmat(Jtmp,[1 1 k]),n,m,p,m*k);
  Jtmp = spblkdiag(permute(Jtmp,[2 3 1]));
  J12  = Jtmp*kron(speye(m),Bsnext);
  
  J = [J11 J12;
       J21 J22];
else % Without Jacobian
  F1 = recsEquilibrium(x,s,zeros(n,0),func,params,grid,c,e,w,fspace,'resapprox-complete');  
end

F = [F1; F2];
