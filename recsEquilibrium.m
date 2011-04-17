function [F,J] = recsEquilibrium(x,s,z,func,params,grid,c,e,w,fspace,method)
% RECSEQUILIBRIUM evaluate the equilibrium equations and Jacobian
%
% RECSEQUILIBRIUM is called by RECSSOLVEEQUILIBRIUM. It is not meant to be called
% directly by the user.
%
% See also RECSSOLVEEQUILIBRIUM.

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

[n,d] = size(s);
x     = reshape(x,[],n)';
m     = size(x,2);
checkjacobian = 0; % Check analytical derivatives against numerical ones

switch method
 case {'expapprox','resapprox-simple'}
  if nargout==2
    output  = struct('F',1,'Js',0,'Jx',1,'Jz',0);
    [F,~,J] = func('f',s,x,z,[],[],[],params,output);
    J       = permute(J,[2 3 1]);
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
    output = struct('F',1,'Js',0,'Jx',0,'Jz',0);
    F      = func('f',s,x,z,[],[],[],params,output);
  end
  F        = reshape(F',n*m,1);
 case {'expfunapprox','resapprox-complete'}
  K     = size(e,1);
  ind   = (1:n);
  ind   = ind(ones(1,K),:);
  ss    = s(ind,:);
  xx    = x(ind,:);
  ee    = e(repmat(1:K,1,n),:);

  if (nargout==2) && (((nargout(func)>=4) && isequal(method,'resapprox-complete')) ...
                      || ((nargout(func)>=3) && isequal(method,'expfunapprox'))) % Analytic derivatives
    output              = struct('F',1,'Js',0,'Jx',1);
    [snext,~,gx]        = func('g',ss,xx,[],ee,[],[],params,output);

    switch method
     case 'expfunapprox'
      H                 = funeval(c,fspace,snext,[zeros(1,d); eye(d)]);
      if nargout(func)==6
        output = struct('F',0,'Js',0,'Jx',0,'Jsn',0,'Jxn',0,'hmult',1);
        [~,~,~,~,~,hmult] = func('h',[],[],[],ee,snext,zeros(size(snext,1),m),params,output);
        H               = H.*hmult(:,:,ones(1+d,1));
      end
      h   = H(:,:,1);
      hds = H(:,:,2:end);
     case 'resapprox-complete'
      [LB,UB]              = func('b',snext,[],[],[],[],[],params);
      Xnext                = funeval(c,fspace,snext,[zeros(1,d); eye(d)]);
      xnext                = min(max(Xnext(:,:,1),LB),UB);
      xnextds              = Xnext(:,:,2:end);

      output = struct('F',1,'Js',0,'Jx',1,'Jsn',1,'Jxn',1,'hmult',1);
      if nargout(func)<6
        [h,~,hx,hsnext,hxnext]       = func('h',ss,xx,[],ee,snext,xnext,params,output);
      else
        [h,~,hx,hsnext,hxnext,hmult] = func('h',ss,xx,[],ee,snext,xnext,params,output);
        h      = h.*hmult;
        hx     = hx.*hmult(:,:,ones(m,1));
        hsnext = hsnext.*hmult(:,:,ones(d,1));
        hxnext = hxnext.*hmult(:,:,ones(m,1));
      end
    end
    p           = size(h,2);
    z           = reshape(w'*reshape(h,K,n*p),n,p);
    output      = struct('F',1,'Js',0,'Jx',1,'Jz',1);
    [F,~,fx,fz] = func('f',s,x,z,[],[],[],params,output);
    F           = reshape(F',n*m,1);

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
  output = struct('F',1,'Js',0,'Jx',0,'Jz',0);
  F      = func('f',s,x,z,[],[],[],params,output);
else
  K       = size(e,1);
  ind     = (1:n);
  ind     = ind(ones(1,K),:);
  ss      = s(ind,:);
  xx      = x(ind,:);
  ee      = e(repmat(1:K,1,n),:);

  output  = struct('F',1,'Js',0,'Jx',0);
  snext   = func('g',ss,xx,[],ee,[],[],params,output);

  switch method
   case 'expfunapprox'
    h                 = funeval(c,fspace,snext);
    if nargout(func)==6
      output  = struct('F',0,'Js',0,'Jx',0,'Jsn',0,'Jxn',0,'hmult',1);
      [~,~,~,~,~,hmult] = func('h',[],[],[],ee,snext,zeros(size(snext,1),m),params,output);
      h                 = h.*hmult;
    end

   case 'resapprox-complete'
    [LB,UB] = func('b',snext,[],[],[],[],[],params);
    xnext   = min(max(funeval(c,fspace,snext),LB),UB);
    output  = struct('F',1,'Js',0,'Jx',0,'Jsn',0,'Jxn',0,'hmult',1);
    if nargout(func)<6
      h                 = func('h',ss,xx,[],ee,snext,xnext,params,output);
    else
      [h,~,~,~,~,hmult] = func('h',ss,xx,[],ee,snext,xnext,params,output);
      h                 = h.*hmult;
    end

  end
  p       = size(h,2);
  z       = reshape(w'*reshape(h,K,n*p),n,p);
  output  = struct('F',1,'Js',0,'Jx',0,'Jz',0);
  F       = func('f',s,x,z,[],[],[],params,output);
end
F = reshape(F',n*m,1);

return