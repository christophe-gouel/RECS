function [F,J] = SSResidual(X,func,params,e,d,m)

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

if nargout==2
  J = sparse(fdjac(@residual_function,X,func,params,e,d,m));
end
F = residual_function(X,func,params,e,d,m);

function F = residual_function(X,func,params,e,d,m)
ss = X(1:d)';
xx = X(d+1:d+m)';
output = struct('F',1,'Js',0,'Jx',0,'Jz',0,'Jsn',0,'Jxn',0);
zz = feval(func,'h',ss,xx,[],e,ss,xx,params,output);
g  = feval(func,'g',ss,xx,[],e,[],[],params,output);
f  = feval(func,'f',ss,xx,zz,[],[],[],params,output);
F  = [ss-g f]';
return
