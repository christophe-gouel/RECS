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

switch method
 case {'expapprox','resapprox-simple'}
  if nargout==2 % With Jacobian
    output  = struct('F',1,'Js',0,'Jx',1,'Jz',0);
    [F,~,J] = func('f',s,x,z,[],[],[],params,output);
  else % Without Jacobian
    output = struct('F',1,'Js',0,'Jx',0,'Jz',0);
    F      = func('f',s,x,z,[],[],[],params,output);
  end
 case {'expfunapprox','resapprox-complete'}
  k     = size(e,1);
  ind   = (1:n);
  ind   = ind(ones(1,k),:);
  ss    = s(ind,:);
  xx    = x(ind,:);
  ee    = e(repmat(1:k,1,n),:);

  if nargout==2 % With Jacobian
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
      h  = H(:,:,1);
      hs = H(:,:,2:end);
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
    z           = reshape(w'*reshape(h,k,n*p),n,p);
    output      = struct('F',1,'Js',0,'Jx',1,'Jz',1);
    [F,~,fx,fz] = func('f',s,x,z,[],[],[],params,output);

    switch method
     case 'expfunapprox'
      Jtmp = arraymult(hs,gx,k*n,p,d,m);
     case 'resapprox-complete'
      Jtmp = hx+arraymult(hsnext+arraymult(hxnext,xnextds,k*n,p,m,d),gx,k*n,p,d,m);
    end
    Jtmp = reshape(w'*reshape(Jtmp,k,n*p*m),n,p,m);
    J    = fx+arraymult(fz,Jtmp,n,m,p,m);
  else % Without Jacobian
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
    z       = reshape(w'*reshape(h,k,n*p),n,p);
    output  = struct('F',1,'Js',0,'Jx',0,'Jz',0);
    F       = func('f',s,x,z,[],[],[],params,output);
  end
end

F = reshape(F',n*m,1);
if nargout==2
  J = permute(J,[2 3 1]);
  J = spblkdiag(J,grid); 
end
