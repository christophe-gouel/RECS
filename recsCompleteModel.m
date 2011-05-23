function [F,J] = recsCompleteModel(X,interp,model)

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

e       = model.e;
[d,m,p] = model.dim{:};
func    = model.func;
params  = model.params;
w       = model.w;

fspace = interp.fspace;
Phi    = interp.Phi;
s      = interp.s;

n = size(s,1);

if (nargout==2) && (nargout(func)==4) % Analytic derivatives

  x = reshape(X(1:n*m),m,n)';
  z = reshape(X(n*m+1:n*(m+p)),p,n)';
  c = reshape(X(n*(m+p)+1:n*(m+2*p)),p,n)';

  K   = size(e,1);
  ind = (1:n);
  ind = ind(ones(1,K),:);
  ss  = s(ind,:);
  xx  = x(ind,:);
  ee  = e(repmat(1:K,1,n),:);

  [f.eq,f.eqx,f.eqz] = feval(func,'f',s,x,z,[],[],[],params{:}); % (n,m) (n,m,m) (n,m,p)
  [snext,gx]         = feval(func,'g',ss,xx,[],ee,[],[],params{:}); % (n*K,p) (n*K,p,m)
  Phisnext           = funbas(fspace,snext);
  f.exp              = z-reshape(w'*reshape(Phisnext*c,K,n*p),n,p); % (n,p)
  zds                = funeval(c,fspace,snext,eye(d));
  [h,~,~,hxnext]     = feval(func,'h',[],[],[],[],s,x,params{:}); % (n,p) (n,p,m)
  f.ratexp           = h-Phi*c; % (n,p)

  F = [reshape(f.eq',n*m,1)
       reshape(f.exp',n*p,1)
       reshape(f.ratexp',n*p,1) ]; % (n*(m+2p),1)

  J21          = arraymult(zds,gx,K*n,p,d,m);
  J21          = reshape(w'*reshape(J21,K,n*p*m),n,p,m);
  [J21,gridpm] = spblkdiag(permute(-J21,[2 3 1]));

  J23 = kron(speye(p),Phisnext);
  J23 = -reshape(w'*reshape(J23,K,n^2*p^2),n*p,n*p);

  J = [spblkdiag(permute(f.eqx,[2 3 1]))         spblkdiag(permute(f.eqz,[2 3 1])) sparse(n*m,n*p)
       J21                                       speye(n*p)                        J23
       spblkdiag(permute(hxnext,[2 3 1]),gridpm) sparse(n*p,n*p)                   -kron(speye(p),Phi)];
  if 0 % Check analytical derivatives against numerical ones
    Jnum = fdjac(@modelfunction,X,s,n,m,p,func,e,w,params,fspace,Phi);
    normest(J-Jnum)
% $$$     normest(J(1+n*m:n*(m+p),1:m*n)-Jnum(1+n*m:n*(m+p),1:m*n))
% $$$     normest(J(1+n*m:n*(m+p),1+m*n:(m+p)*n)-Jnum(1+n*m:n*(m+p),1+m*n:(m+p)*n))
% $$$     normest(J(1+n*m:n*(m+p),1+(m+p)*n:(m+2*p)*n)-Jnum(1+n*m:n*(m+p),1+(m+p)*n:(m+2*p)*n))
% $$$     normest(J(1+n*(m+p):n*(m+2*p),1:m*n)-Jnum(1+n*(m+p):n*(m+2*p),1:m*n))
% $$$     normest(J(1+n*(m+p):n*(m+2*p),1+n*m:(m+p)*n)-Jnum(1+n*(m+p):n*(m+2*p),1+n*m:(m+p)*n))
% $$$     normest(J(1+n*(m+p):n*(m+2*p),1+n*(m+p):(m+2*p)*n)-Jnum(1+n*(m+p):n*(m+2*p),1+n*(m+p):(m+2*p)*n))
  end
elseif (nargout==2) && (nargout(func)<4) % Numerical derivatives
  F = modelfunction(X,s,n,m,p,func,e,w,params,fspace,Phi);
  J = sparse(fdjac(@modelfunction,X,s,n,m,p,func,e,w,params,fspace,Phi));
else
  F = modelfunction(X,s,n,m,p,func,e,w,params,fspace,Phi);
end

function F = modelfunction(X,s,n,m,p,func,e,w,params,fspace,Phi)

x = reshape(X(1:n*m),m,n)';
z = reshape(X(n*m+1:n*(m+p)),p,n)';
c = reshape(X(n*(m+p)+1:n*(m+2*p)),p,n)';

K   = size(e,1);
ind = (1:n);
ind = ind(ones(1,K),:);
ss  = s(ind,:);
xx  = x(ind,:);
ee  = e(repmat(1:K,1,n),:);

f.eq     = feval(func,'f',s,x,z,[],[],[],params{:}); % (n,m)
snext    = feval(func,'g',ss,xx,[],ee,[],[],params{:}); % (n*K,p)
f.exp    = z-reshape(w'*reshape(funeval(c,fspace,snext),K,n*p),n,p); % (n,p)
f.ratexp = feval(func,'h',[],[],[],[],s,x,params{:})-Phi*c; % (n,p)

F = [reshape(f.eq',n*m,1)
     reshape(f.exp',n*p,1)
     reshape(f.ratexp',n*p,1) ]; % (n*(m+2p),1)
return