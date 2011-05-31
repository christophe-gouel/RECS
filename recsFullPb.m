function [F,J] = recsFullPb(X,s,func,params,grid,e,w,fspace,method,Phi,m,functional)
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
   
   case 'resapprox-complete'
    J21 = -speye(n*m);

    J22 = kron(B,speye(m));

  end

  % Calculation of J12
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
      h                 = h.*hmult;
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

    p           = size(h,2);

  end
  z          = reshape(w'*reshape(h,k,n*p),n,p);
  [~,~,~,fz] = func('f',s,x,z,[],[],[],params,struct('F',0,'Js',0,'Jx',0,'Jz',1));

  Bsnext = mat2cell(Bsnext',n,k*ones(n,1))';
  switch method
   case 'expfunapprox'     
    [~,gridJ12] = spblkdiag(zeros(1,n,p),[],0);
    J12    = cellfun(@(X) full(spblkdiag((X*w)',gridJ12,1,p)),...
                     Bsnext,'UniformOutput',false);
   case 'resapprox-complete'
    [~,gridJ12] = spblkdiag(zeros(p,m,k),[],0);
    kw     = kron(w',eye(p));
    hxnext = num2cell(reshape(hxnext,[n k p m]),[2 3 4]);
    J12    = cellfun(@(X,Y) kw*spblkdiag(permute(X,[3 4 2 1]),gridJ12)*kron(Y',speye(m)),...
                     hxnext,Bsnext,'UniformOutput',false);
  end
  J12    = reshape(cat(1,J12{:}),[n p numel(c)]);
    
  J12 = arraymult(fz,J12,n,m,p,numel(c));
  J12 = reshape(permute(J12,[2 1 3]),[n*m numel(c)]);

  J = [J11 J12;
       J21 J22];
else % Without Jacobian
  F1 = recsEquilibrium(reshape(x',[n*m 1]),...
                       s,zeros(n,0),func,params,grid,c,e,w,fspace,method);
end

F = [F1; F2];
