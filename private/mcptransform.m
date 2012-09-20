function [out1,out2,out3,out4,out5] = mcptransform(func,flag,s,x,z,e,snext,xnext,params,output,ix,nx,w,v)
% MCPTRANSFORM Transform a MCP problem in which bounds are endogenous into a problem with exogenous bounds

% Copyright (C) 2011-2012 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

voidcell                   = cell(1,5);
[out1,out2,out3,out4,out5] = voidcell{:};

switch flag

  case 'b'
    [LBx,UBx] = func('b',s,[],[],[],[],[],params);
    n                    = size(s,1);
    LBxmod               = -inf(size(LBx));
    LBxmod(:,~(ix(:,1))) = LBx(:,~(ix(:,1)));
    UBxmod               = +inf(size(UBx));
    UBxmod(:,~(ix(:,2))) = UBx(:,~(ix(:,2)));
    out1                 = [LBxmod zeros(n,nx(1)+nx(2))];
    out2                 = [UBxmod  +inf(n,nx(1)+nx(2))];
    
  case 'f'
    [n,d]     = size(s);
    [LBx,UBx] = func('b',s,[],[],[],[],[],params);

    % F
    [f,fs,fx,fz] = func('f',s,x,z,[],[],[],params,output);

    f(:,ix(:,1)) = f(:,ix(:,1))-w;
    f(:,ix(:,2)) = f(:,ix(:,2))+v;
    out1 = [f x(:,ix(:,1))-LBx(:,ix(:,1)) UBx(:,ix(:,2))-x(:,ix(:,2))];
    
    % dF/ds
    if output.Js
      dLBxds = zeros(n,nx(1),d);
      dUBxds = zeros(n,nx(2),d);
      if sum(nx)
        for i=1:n
          if nx(1), dLBxds(i,:,:) = numjac(@(S) Bounds(func,S,params,1,ix),s(i,:)); end
          if nx(2), dUBxds(i,:,:) = numjac(@(S) Bounds(func,S,params,2,ix),s(i,:)); end
        end
      end
      out2 = cat(2,fs,-dLBxds,dUBxds);
    end
    
    % dF/dX
    if output.Jx
      m    = size(x,2);
      dfdw = zeros(m,nx(1));
      if nx(1), dfdw(sub2ind(size(dfdw),find(ix(:,1)),1:nx(1))) = 1; end
      dfdw = permute(dfdw(:,:,ones(n,1)),[3 1 2]);
      dfdv = zeros(m,nx(2));
      if nx(2), dfdv(sub2ind(size(dfdv),find(ix(:,2)),1:nx(2))) = 1; end
      dfdv = permute(dfdv(:,:,ones(n,1)),[3 1 2]);
      out3 = cat(2,...
                 cat(3,fx                    ,-dfdw               ,dfdv                ),...
                 cat(3, permute(dfdw,[1 3 2]),zeros(n,nx(1),nx(1)+nx(2))),...
                 cat(3,-permute(dfdv,[1 3 2]),zeros(n,nx(2),nx(1)+nx(2))));
    end
    
    % dF/dz
    if output.Jz
      p    = size(fz,3);
      out4 = cat(2,fz,zeros(n,nx(1)+nx(2),p));
    end

  case 'g'
    [n,d]             = size(s);
    [out1,out2,gx]    = func('g',s,x,[],e ,[],[],params,output);
    
    if output.Jx
      out3 = cat(3,gx,zeros(n,d,nx(1)+nx(2)));
    end
    
  case 'h'
    [n,d]        = size(snext);
    m            = size(xnext,2);
    output.hmult = 1;
    if nargout(func)<6
      [h,hs,hx,hsnext,hxnext]       = func('h',s,x,[],e,snext,xnext,params,output);
    else
      [h,hs,hx,hsnext,hxnext,hmult] = func('h',s,x,[],e,snext,xnext,params,output);
      h      = h.*hmult;
      if output.Js, hs     = hs.*hmult(:,:,ones(d,1)); end
      if output.Jx, hx     = hx.*hmult(:,:,ones(m,1)); end
      if output.Jsn, hsnext = hsnext.*hmult(:,:,ones(d,1)); end
      if output.Jxn, hxnext = hxnext.*hmult(:,:,ones(m,1)); end
    end
    out1 = h;
    out2 = hs;
    out4 = hsnext;
    p    = size(h,2);
    if output.Jx,  out3 = cat(3,hx,zeros(n,p,nx(1)+nx(2))); end
    if output.Jxn, out5 = cat(3,hxnext,zeros(n,p,nx(1)+nx(2))); end

end
