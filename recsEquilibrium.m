function [F,J] = recsEquilibrium(x,s,z,func,params,grid,c,e,w,fspace,method)
% RECSEQUILIBRIUM evaluate the equilibrium equations and Jacobian

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

[n,d] = size(s);
x     = reshape(x,[],n)';
m     = size(x,2);
checkjacobian = 0; % Check analytical derivatives against numerical ones

switch method
 case {'expapprox','resapprox-simple'}
  if nargout==2
    [F,J] = feval(func,'f',s,x,z,[],[],[],params);
    J     = permute(J,[2 3 1]);
    if isempty(J) || (checkjacobian) % Numerical derivatives
      Jnum   = zeros(m,m,n);
      for i=1:n
        Jnum(:,:,i) = fdjac(@equilibriumfunction,x(i,:),s(i,:),z(i,:),func, ...
                            params,[],[],[],[],method);
      end
      if (checkjacobian)
        norm(J(:)-Jnum(:))
      else
        J = Jnum;
      end
    end
  else
    F     = feval(func,'f',s,x,z,[],[],[],params);
  end
  F       = reshape(F',n*m,1);
 case {'expfunapprox','resapprox-complete'}
  K     = size(e,1);
  ind   = (1:n);
  ind   = ind(ones(1,K),:);
  ss    = s(ind,:);
  xx    = x(ind,:);
  ee    = e(repmat(1:K,1,n),:);

  if (nargout==2) && (((nargout(func)>=4) && isequal(method,'resapprox-complete')) ...
                      || ((nargout(func)>=3) && isequal(method,'expfunapprox'))) % Analytic derivatives
    [snext,gx]           = feval(func,'g',ss,xx,[],ee,[],[],params);

    switch method
     case 'expfunapprox'
      H                 = funeval(c,fspace,snext,[zeros(1,d); eye(d)]);
      if nargout(func)==5
        [~,~,~,~,hmult] = feval(func,'h',[],[],[],ee,snext,zeros(size(snext,1),m),params);
        H               = H.*hmult(:,:,ones(1+d,1));
      end
      h   = H(:,:,1);
      hds = H(:,:,2:end);
     case 'resapprox-complete'
      [LB,UB]              = feval(func,'b',snext,[],[],[],[],[],params);
      Xnext                = funeval(c,fspace,snext,[zeros(1,d); eye(d)]);
      xnext                = min(max(Xnext(:,:,1),LB),UB);
      xnextds              = Xnext(:,:,2:end);

      if nargout(func)<5
        [h,hx,hsnext,hxnext] = feval(func,'h',ss,xx,[],ee,snext,xnext,params);
      else
        [h,hx,hsnext,hxnext,hmult] = feval(func,'h',ss,xx,[],ee,snext,xnext,params);
        h      = h.*hmult;
        hx     = hx.*hmult(:,:,ones(m,1));
        hsnext = hsnext.*hmult(:,:,ones(d,1));
        hxnext = hxnext.*hmult(:,:,ones(m,1));
      end
    end
    p         = size(h,2);
    z         = reshape(w'*reshape(h,K,n*p),n,p);
    [F,fx,fz] = feval(func,'f',s,x,z,[],[],[],params);
    F         = reshape(F',n*m,1);

    switch method
     case 'expfunapprox'
      Jtmp = arraymult(hds,gx,K*n,p,d,m);
% $$$       Jtmp = multiprod(hds,gx,[2 3],[2 3]);
% $$$       Jtmp = mtimesx(permute(hds,[2 3 1]),permute(gx,[2 3 1]));
% $$$       Jtmp = permute(Jtmp,[3 1 2]);
     case 'resapprox-complete'
      Jtmp = hx+arraymult(hsnext+arraymult(hxnext,xnextds,K*n,p,m,d),gx,K*n,p,d,m);
% $$$       Jtmp = hx+multiprod(hsnext+multiprod(hxnext,xnextds,[2 3]),gx,[2 3],[2 3]);
% $$$       Jtmp = hx+permute(mtimesx(permute(hsnext,[2 3 1])+mtimesx(permute(hxnext,[2 3 1]),permute(xnextds,[2 ...
% $$$                           3 1])),permute(gx,[2 3 1])),[3 1 2]);
    end
    Jtmp = reshape(w'*reshape(Jtmp,K,n*p*m),n,p,m);
    J    = fx+arraymult(fz,Jtmp,n,m,p,m);
% $$$     J    = fx+multiprod(fz,Jtmp,[2 3],[2 3]);
% $$$     J    = fx+permute(mtimesx(permute(fz,[2 3 1]),permute(Jtmp,[2 3 1])),[3 1 2]);
    J    = permute(J,[2 3 1]);
    if (checkjacobian)
      Jnum                   = zeros(m,m,n);
      for i=1:n
        Jnum(:,:,i) = fdjac(@equilibriumfunction,x(i,:),s(i,:),[],func, ...
                            params,c,e,w,fspace,method);
      end
      norm(J(:)-Jnum(:))
    end
  elseif (nargout==2) && (((nargout(func)<4) && isequal(method,'resapprox-complete')) ...
                          || ((nargout(func)<3) && isequal(method,'expfunapprox'))) % Numerical derivatives
    F = equilibriumfunction(reshape(x',[n*m 1]),s,[],func,params,c,e,w,fspace,method);
    J                   = zeros(m,m,n);
    for i=1:n
      J(:,:,i) = fdjac(@equilibriumfunction,x(i,:),s(i,:),[],func,params,c,e,w, ...
                       fspace,method);
    end
  else
    F = equilibriumfunction(reshape(x',[n*m 1]),s,[],func,params,c,e,w,fspace,method);
  end
end

if nargout==2, J = spblkdiag(J,grid); end


function F = equilibriumfunction(x,s,z,func,params,c,e,w,fspace,method)

n = size(s,1);
x = reshape(x,[],n)';
m = size(x,2);

if ~isempty(z) % Cases expapprox and resapprox-simple
  F = feval(func,'f',s,x,z,[],[],[],params);
else
  K       = size(e,1);
  ind     = (1:n);
  ind     = ind(ones(1,K),:);
  ss      = s(ind,:);
  xx      = x(ind,:);
  ee      = e(repmat(1:K,1,n),:);

  snext   = feval(func,'g',ss,xx,[],ee,[],[],params);

  switch method
   case 'expfunapprox'
    h                 = funeval(c,fspace,snext);
    if nargout(func)==5
      [~,~,~,~,hmult] = feval(func,'h',[],[],[],ee,snext,zeros(size(snext,1),m),params);
      h               = h.*hmult;
    end

   case 'resapprox-complete'
    [LB,UB] = feval(func,'b',snext,[],[],[],[],[],params);
    xnext   = min(max(funeval(c,fspace,snext),LB),UB);
    if nargout(func)<5
      h = feval(func,'h',ss,xx,[],ee,snext,xnext,params);
    else
      [h,~,~,~,hmult] = feval(func,'h',ss,xx,[],ee,snext,xnext,params);
      h               = h.*hmult;
    end

  end
  p       = size(h,2);
  z       = reshape(w'*reshape(h,K,n*p),n,p);
  F       = feval(func,'f',s,x,z,[],[],[],params);
end
F = reshape(F',n*m,1);

return