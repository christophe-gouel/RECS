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
  case 'expfunapprox'
    p  = size(c,2);
    if nargout==1
      output = struct('F',1,'Js',0,'Jx',0,'Jsn',0,'Jxn',0,'hmult',0);
    if NewtonMethod, R = funeval(c,fspace,Phi)-h([],[],[],s,x,params,output);
    else             R = funfitxy(fspace,Phi,h([],[],[],s,x,params,output))-c;
    end
      R = reshape(R',n*p,1);
    end
    
  case 'resapprox'
    if NewtonMethod, R = funeval(c,fspace,Phi)-x(:,ixforward);
    else             R = funfitxy(fspace,Phi,x(:,ixforward))-c;
    end
    R = reshape(R',n*mf,1);

end

%% Evaluate Jacobian if required
if nargout>=2
  B = funbas(fspace,s);
  
  switch funapprox
    case 'expfunapprox'
      output = struct('F',1,'Js',0,'Jx',0,'Jsn',0,'Jxn',1,'hmult',0);
      [hv,~,~,~,hxnext] = h([],[],[],s,x,params,output);
      R  = funeval(c,fspace,Phi)-hv;
      R  = reshape(R',n*p,1);
      Rx = -spblkdiag(permute(hxnext,[2 3 1]));
      
      Rc = kron(B,speye(p));
      
    case 'resapprox'
      if m==mf
        Rx = -speye(n*mf);
      else
        indx = 1:n*m;
        Rx = sparse(1:n*mf,indx(repmat(ixforward,1,n)),-ones(n*mf,1),n*mf,n*m,n*mf);
      end
      
      Rc = kron(B,speye(mf));
      
  end
end