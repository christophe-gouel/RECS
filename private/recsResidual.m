function [R,Rx,Rc] = recsResidual(s,x,func,params,c,fspace,funapprox,Phi)
% RECSRESIDUAL evaluates the residual of rational expectations and its Jacobians
%
% RECSRESIDUAL is called by recsSolveREEFull and recsSolveREEIterNewton. It is not
% meant to be called directly by the user.
%
% See also RECSSOLVEREEFULL, RECSSOLVEREEITERNEWTON.

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%%
[n,m] = size(x);

switch funapprox
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

if nargout>=2
  %% With Jacobian
  B = funbas(fspace,s);

  switch funapprox
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
end