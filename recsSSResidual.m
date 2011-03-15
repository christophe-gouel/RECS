function [F,J] = recsSSResidual(X,func,params,e,d,m)
% RECSSSRESIDUAL

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

if nargout==2
  J = sparse(fdjac(@residual_function,X,func,params,e,d,m));
end
F = residual_function(X,func,params,e,d,m);

function F = residual_function(X,func,params,e,d,m)
ss = X(1:d)';
xx = X(d+1:d+m)';
zz = func('h',ss,xx,[],e,ss,xx,params);
g  = func('g',ss,xx,[],e,[],[],params);
f  = func('f',ss,xx,zz,[],[],[],params);
F  = [ss-g f]';
return
