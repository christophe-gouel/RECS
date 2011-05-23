function [F,J] = recsFullPb2(X,s,func,params,grid,e,w,fspace,Phi)
% RECSFULLPB 
%
% This code is not optimized. Some formulas are calculated both in this file and
% in recsEquilibrium. It could be improved.
  
% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt
  
[n,d] = size(s);
x     = X(1:length(X)/2);
x     = reshape(x,[],n)';
m     = size(x,2);
c     = reshape(X(length(X)/2+1:end),m,n)';

F2 = funeval(c,fspace,Phi)-x;
F2 = reshape(F2',n*m,1);

if nargout==2 % With Jacobian
  [F1,J11] = recsEquilibrium(reshape(x',[n*m 1]),...
                             s,zeros(n,0),func,params,grid,c,e,w,fspace,'resapprox-complete');  

  B = funbas(fspace,s);

  J21 = -speye(n*m);

  Tmn = TvecMat(m,n);
  J22 = kron(B,speye(m));

  % Calculation of J12
  J12   = sparse(n*m,n*m);
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
  p           = size(h,2);
  z           = reshape(w'*reshape(h,k,n*p),n,p);
  output      = struct('F',1,'Js',0,'Jx',1,'Jz',1);
  [~,~,~,fz] = func('f',s,x,z,[],[],[],params,output);

% $$$   for i=1:n % The kronecker products with identity matrix could be remplaced by
% $$$             % matrix repetition and block diagonal matrix construction
% $$$     hxnextk = spblkdiag(permute(hxnext((i-1)*k+1:i*k,:,:),[2 3 1]));
% $$$ %    J12((i-1)*m+1:i*m,:) = kron(permute(fz(i,:,:),[2 3 1]),w')*hxnextk*kron(speye(m),Bsnext((i-1)*k+1:i*k,:))*Tmn';
% $$$    Btmp = Bsnext((i-1)*k+1:i*k,:);
% $$$    tmp  = zeros(0,0);
% $$$    for j=1:k
% $$$ %     size(kron(speye(m),Btmp(j,:)))
% $$$      tmp = [tmp; kron(speye(m),Btmp(j,:))*Tmn];
% $$$    end
% $$$    J12((i-1)*m+1:i*m,:) = kron(permute(fz(i,:,:),[2 3 1]),w')*hxnextk*tmp;
% $$$   end

  iter = 0;
  for i=1:n % The kronecker products with identity matrix could be remplaced by
            % matrix repetition and block diagonal matrix construction
   tmp  = zeros(p,n*m);
   for j=1:k
     iter = iter+1;
     tmp = tmp+w(j)*permute(hxnext(iter,:,:),[2 3 1])*kron(speye(m),Bsnext(iter,:))*Tmn;
   end
   J12((i-1)*m+1:i*m,:) = permute(fz(i,:,:),[2 3 1])*tmp;
  end

  disp('it') 
% $$$   Jtmp = reshape(w'*reshape(hxnext,k,n*p*m),n,p,m);
% $$$   Jtmp = arraymult(fz,repmat(Jtmp,[1 1 k]),n,m,p,m*k);
% $$$   Jtmp = spblkdiag(permute(Jtmp,[2 3 1]));
% $$$   J12  = Jtmp*kron(speye(m),Bsnext);

% $$$   
% $$$   norm(full(J12-J12test))
% $$$   spy(J12)
% $$$   figure()
% $$$   spy(J12test)
  
  J = [J11 J12;
       J21 J22];
% $$$   nnz(J)
% $$$   nzmax(J)
else % Without Jacobian
  F1 = recsEquilibrium(reshape(x',[n*m 1]),...
                       s,zeros(n,0),func,params,grid,c,e,w,fspace,'resapprox-complete');  
end

F = [F1; F2];
