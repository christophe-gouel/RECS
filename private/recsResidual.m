function [R,Rx,Rc] = recsResidual(s,x,h,params,c,fspace,funapprox,Phi,ixforward,NewtonMethod)
% RECSRESIDUAL evaluates the residual of rational expectations and its Jacobians
%
% RECSRESIDUAL is called by recsSolveREEFull and recsSolveREEIterNewton. It is not
% meant to be called directly by the user.
%
% If called by a Newton solver, the residual equation is defined by
%  R = Phi*c-x,
% because the Jacobian of this equation is easy to calculate. If called by a
% quasi-Newton or a successive approximation solver the residual equation is
%  R = Phi\x-c.
%
% See also RECSSOLVEREEFULL, RECSSOLVEREEITERNEWTON.

% Copyright (C) 2011-2013 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%%
[n,m] = size(x);
mf    = sum(ixforward); % Number of forward response variables

%% Evaluate residual
switch funapprox
  case 'resapprox'
    if NewtonMethod, R = funeval(c,fspace,Phi)-x(:,ixforward);
    else             R = funfitxy(fspace,Phi,x(:,ixforward))-c;
    end
    R = reshape(R',n*mf,1);

  case 'expfunapprox'
    p  = size(c,2);
    if nargout==1
      if NewtonMethod, R = funeval(c,fspace,Phi)-h(zeros(n,0),[],[],s,x,params);
      else             R = funfitxy(fspace,Phi,h(zeros(n,0),[],[],s,x,params))-c;
      end
      R = reshape(R',n*p,1);
    end
    
end

%% Evaluate Jacobian if required
if nargout>=2
  B = funbas(fspace,s);
  
  switch funapprox
    case 'resapprox'
      if m==mf
        Rx = -speye(n*mf);
      else
        indx = 1:n*m;
        Rx = sparse(1:n*mf,indx(repmat(ixforward,1,n)),-ones(n*mf,1),n*mf,n*m,n*mf);
      end
      
      Rc = kron(B,speye(mf));
      
    case 'expfunapprox'
      [hv,~,~,~,~,hxnext] = h(zeros(n,0),[],[],s,x,params,[1 0 0 0 0 1]);
      R  = funeval(c,fspace,Phi)-hv;
      R  = reshape(R',n*p,1);
      Rx = -spblkdiag(permute(hxnext,[2 3 1]));
      
      Rc = kron(B,speye(p));
      
  end
end