function [F,J] = recsFullPb(X,s,func,params,grid,e,w,fspace,method,Phi,m)
% RECSFULLPB
%
% This code is not optimized. Some formulas are calculated both in this file and
% in recsEquilibrium. It could be improved.

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

n     = size(s,1);
x     = X(1:m*n);
x     = reshape(x,m,n)';
c     = reshape(X(m*n+1:end),[],n)';

switch method
 case 'expfunapprox'
  p  = size(c,2);
  if nargout==1
    output = struct('F',1,'Js',0,'Jx',0,'Jsn',0,'Jxn',0,'hmult',0);
    F2 = funeval(c,fspace,Phi)-func('h',[],[],[],[],s,x,params,output);
    F2 = reshape(F2',n*p,1);
  end
 case 'resapprox-complete'
  F2 = funeval(c,fspace,Phi)-x;
  F2 = reshape(F2',n*m,1);
end

if nargout==2 % With Jacobian
  [F1,J11] = recsEquilibrium(reshape(x',[n*m 1]),...
                             s,zeros(n,0),func,params,grid,c,e,w,fspace,method);

  B = funbas(fspace,s);


  switch method
   case 'expfunapprox'
    output = struct('F',1,'Js',0,'Jx',0,'Jsn',0,'Jxn',1,'hmult',0);
    [h,~,~,~,hxnext] = func('h',[],[],[],[],s,x,params,output);
    F2  = funeval(c,fspace,Phi)-h;
    F2  = reshape(F2',n*p,1);
    J21 = -spblkdiag(permute(hxnext,[2 3 1]));

    J22 = kron(B,speye(p));

    sizeJ12 = p;
   case 'resapprox-complete'
    J21 = -speye(n*m);
    J22 = kron(B,speye(m));

    sizeJ12 = m;
  end

  % Calculation of J12
  J12   = sparse([],[],[],n*m,n*sizeJ12,0.1*(n*m*n*sizeJ12)); % TODO: improve the space allocation
  k     = size(e,1);
  ind   = (1:n);
  ind   = ind(ones(1,k),:);
  ss    = s(ind,:);
  xx    = x(ind,:);
  ee    = e(repmat(1:k,1,n),:);
  % ss, xx and ee are organised as
  %   s_i   e_1
  %   s_i   e_k
  %   s_i   e_K
  %   s_i+1 e_1

  output   = struct('F',1,'Js',0,'Jx',0);
  snext    = func('g',ss,xx,[],ee,[],[],params,output);
  Bsnext   = funbas(fspace,snext);

  switch method
   case 'expfunapprox'
    h = Bsnext*c;
    if nargout(func)==6
      output = struct('F',0,'Js',0,'Jx',0,'Jsn',0,'Jxn',0,'hmult',1);
      [~,~,~,~,~,hmult] = func('h',[],[],[],ee,snext,zeros(size(snext,1),m),params,output);
      h               = h.*hmult;
    end
   case 'resapprox-complete'
    [LB,UB]              = func('b',snext,[],[],[],[],[],params);
    xnext                = min(max(Bsnext*c,LB),UB);

    output = struct('F',1,'Js',0,'Jx',0,'Jsn',0,'Jxn',1,'hmult',1);
    if nargout(func)<6
      [h,~,~,~,hxnext]       = func('h',ss,xx,[],ee,snext,xnext,params,output);
    else
      [h,~,~,~,hxnext,hmult] = func('h',ss,xx,[],ee,snext,xnext,params,output);
      h      = h.*hmult;
      hxnext = hxnext.*hmult(:,:,ones(m,1));
    end
  end
  p           = size(h,2);
  z           = reshape(w'*reshape(h,k,n*p),n,p);
  output      = struct('F',1,'Js',0,'Jx',1,'Jz',1);
  [~,~,~,fz] = func('f',s,x,z,[],[],[],params,output);

% $$$   iter = 0;
  for i=1:n % The kronecker products with identity matrix could be remplaced by
            % matrix repetition and block diagonal matrix construction

   ira = (i-1)*k+1:i*k;
  switch method
   case 'expfunapprox'
    tmp = kron(w',eye(p))*kron(Bsnext(ira,:),speye(p));
   case 'resapprox-complete'
    tmp = kron(w',eye(p))*spblkdiag(permute(hxnext(ira,:,:),[2 3 1]))*kron(Bsnext(ira,:),speye(m));
% $$$    tmp  = zeros(p,n*m);
% $$$    for j=1:k
% $$$      iter = iter+1;
% $$$      tmp = tmp+w(j)*permute(hxnext(iter,:,:),[2 3 1])*kron(Bsnext(iter,:),speye(m));
% $$$    end
  end
   J12((i-1)*m+1:i*m,:) = permute(fz(i,:,:),[2 3 1])*tmp;
  end

  J = [J11 J12;
       J21 J22];
else % Without Jacobian
  F1 = recsEquilibrium(reshape(x',[n*m 1]),...
                       s,zeros(n,0),func,params,grid,c,e,w,fspace,method);
end

F = [F1; F2];
