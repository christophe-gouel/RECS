function [F,Jx,Jc] = recsEquilibrium(x,s,z,func,params,gridJx,c,e,w,fspace,method,extrapolate)
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
    [F,~,Jx] = func('f',s,x,z,[],[],[],params,output);
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

  if nargout>=2 % With Jacobian
    output              = struct('F',1,'Js',0,'Jx',1);
    [snext,~,gx]        = func('g',ss,xx,[],ee,[],[],params,output);
    if extrapolate, snextinterp = snext;
    else            snextinterp = max(min(snext,fspace.b),fspace.a); end
    Bsnext = funbasx(fspace,snextinterp,[zeros(1,d); eye(d)]);

    switch method
     case 'expfunapprox'
      H                 = funeval(c,fspace,Bsnext,[zeros(1,d); eye(d)]);
      if nargout(func)==6
        output = struct('F',0,'Js',0,'Jx',0,'Jsn',0,'Jxn',0,'hmult',1);
        [~,~,~,~,~,hmult] = func('h',[],[],[],ee,snext,zeros(size(snext,1),m),params,output);
        H               = H.*hmult(:,:,ones(1+d,1));
      end
      h  = H(:,:,1);
      hs = H(:,:,2:end);
     case 'resapprox-complete'
      [LB,UB]              = func('b',snext,[],[],[],[],[],params);
      Xnext                = funeval(c,fspace,Bsnext,[zeros(1,d); eye(d)]);
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
      Jxtmp = arraymult(hs,gx,k*n,p,d,m);
     case 'resapprox-complete'
      Jxtmp = hx+arraymult(hsnext+arraymult(hxnext,xnextds,k*n,p,m,d),gx,k*n,p,d,m);
    end
    Jxtmp = reshape(w'*reshape(Jxtmp,k,n*p*m),n,p,m);
    Jx    = fx+arraymult(fz,Jxtmp,n,m,p,m);
    
    if nargout==3
      if ~strcmp(Bsnext.format,'expanded'), Bsnext = funbconv(Bsnext,zeros(1,d)); end
      Bsnext = mat2cell(Bsnext.vals{1}',n,k*ones(n,1))';
      switch method % The product with fz is vectorized for 'expfunapprox' and
                    % executed in a loop for resapprox-complete', because it
                    % seems to be the fastest ways to do it.
       case 'expfunapprox'     
        [~,gridJc] = spblkdiag(zeros(1,n,p),[],0);
        Jc    = cellfun(@(X) full(spblkdiag((X*w)',gridJc,1,p)),...
                        Bsnext,'UniformOutput',false);
        Jc    = permute(reshape(cat(1,Jc{:}),[p n numel(c)]),[2 1 3]);     
        Jc    = arraymult(fz,Jc,n,m,p,numel(c));
        Jc    = reshape(permute(Jc,[2 1 3]),[n*m numel(c)]);
       case 'resapprox-complete'
        [~,gridJc] = spblkdiag(zeros(p,m,k),[],0);
        kw     = kron(w',eye(p));
        hxnext = num2cell(reshape(hxnext,[k n p m]),[1 3 4])';
        Jctmp  = cellfun(@(X,Y) kw*spblkdiag(permute(X,[3 4 2 1]),gridJc)*kron(Y',speye(m)),...
                            hxnext,Bsnext,'UniformOutput',false);
        Jc   = zeros(n*m,numel(c));
        for i=1:n
          Jc((i-1)*m+1:i*m,:) = permute(fz(i,:,:),[2 3 1])*Jctmp{i};
        end
      end
    end
  else % Without Jacobian
    output  = struct('F',1,'Js',0,'Jx',0);
    snext   = func('g',ss,xx,[],ee,[],[],params,output);

    switch method
     case 'expfunapprox'
      if extrapolate, snextinterp = snext;
      else,           snextinterp = max(min(snext,fspace.b),fspace.a); end
      h                   = funeval(c,fspace,snextinterp);
      if nargout(func)==6
        output  = struct('F',0,'Js',0,'Jx',0,'Jsn',0,'Jxn',0,'hmult',1);
        [~,~,~,~,~,hmult] = func('h',[],[],[],ee,snext,zeros(size(snext,1),m),params,output);
        h                 = h.*hmult;
      end

     case 'resapprox-complete'
      [LB,UB] = func('b',snext,[],[],[],[],[],params);
      if extrapolate, snextinterp = snext;
      else,           snextinterp = max(min(snext,fspace.b),fspace.a); end
      xnext   = min(max(funeval(c,fspace,snextinterp),LB),UB);
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
if nargout>=2
  Jx = permute(Jx,[2 3 1]);
  Jx = spblkdiag(Jx,gridJx);
end
