function model = mcptransform(model,ix,nx)
% MCPTRANSFORM Transform a MCP problem in which bounds are endogenous into a problem with exogenous bounds

% Copyright (C) 2011-2013 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
b = model.b;
f = model.f;
g = model.g;
h = model.h;

model.bp = @(s,params) bp(s,params,b,ix,nx);
model.fp = @(s,x,w,v,z,params,output) fp(s,x,w,v,z,params,output,b,f,ix,nx);
model.gp = @(s,x,e,params,output) gp(s,x,e,params,output,g,nx);
model.hp = @(s,x,e,snext,xnext,params,output) hp(s,x,e,snext,xnext,params,output,h,nx);

function [LB,UB] = bp(s,params,b,ix,nx)
%% BP

[LBx,UBx]            = b(s,params);
n                    = size(s,1);
LBxmod               = -inf(size(LBx));
LBxmod(:,~(ix(:,1))) = LBx(:,~(ix(:,1)));
UBxmod               = +inf(size(UBx));
UBxmod(:,~(ix(:,2))) = UBx(:,~(ix(:,2)));
LB                   = [LBxmod zeros(n,nx(1)+nx(2))];
UB                   = [UBxmod  +inf(n,nx(1)+nx(2))];

return

function [fpval,fps,fpx,fpz] = fp(s,x,w,v,z,params,output,b,f,ix,nx)
%% FP
[n,d]     = size(s);

% F
[fval,fs,fx,fz] = f(s,x,z,params,output);

% dF/ds 
if output.Js
  if sum(nx)
    [LBx,UBx,dLBxds,dUBxds] = b(s,params);
    dLBxds = dLBxds(:,ix(:,1),:);
    dUBxds = dUBxds(:,ix(:,2),:); 
  else
    [LBx,UBx] = b(s,params);
    dLBxds    = zeros(n,0,d);
    dUBxds    = zeros(n,0,d);
  end
  fps = cat(2,fs,-dLBxds,dUBxds);
else
  [LBx,UBx] = b(s,params);
end % if output.Js

fval(:,ix(:,1)) = fval(:,ix(:,1))-w;
fval(:,ix(:,2)) = fval(:,ix(:,2))+v;
fpval = [fval x(:,ix(:,1))-LBx(:,ix(:,1)) UBx(:,ix(:,2))-x(:,ix(:,2))];

% dF/dX
if output.Jx
  m    = size(x,2);
  dfdw = zeros(m,nx(1));
  if nx(1), dfdw(sub2ind(size(dfdw),int16(find(ix(:,1)))',1:nx(1))) = 1; end
  dfdw = permute(dfdw(:,:,ones(n,1)),[3 1 2]);
  dfdv = zeros(m,nx(2));
  if nx(2), dfdv(sub2ind(size(dfdv),int16(find(ix(:,2)))',1:nx(2))) = 1; end
  dfdv = permute(dfdv(:,:,ones(n,1)),[3 1 2]);
  fpx = cat(2,...
            cat(3,fx                    ,-dfdw               ,dfdv                ),...
            cat(3, permute(dfdw,[1 3 2]),zeros(n,nx(1),nx(1)+nx(2))),...
            cat(3,-permute(dfdv,[1 3 2]),zeros(n,nx(2),nx(1)+nx(2))));
end % if output.Jx
    
% dF/dz
if output.Jz
  p   = size(fz,3);
  fpz = cat(2,fz,zeros(n,nx(1)+nx(2),p));
end % if output.Jz

return

function [gpval,gps,gpx] = gp(s,x,e,params,output,g,nx)
%% GP

[n,d]             = size(s);
[gpval,gps,gx]    = g(s,x,e,params,output);

if output.Jx
  gpx = cat(3,gx,zeros(n,d,nx(1)+nx(2)));
end

return

function [hpval,hps,hpx,hpsn,hpxn] = hp(s,x,e,snext,xnext,params,output,h,nx)
%% HP

[n,d]        = size(snext);
m            = size(xnext,2);
output.hmult = 1;
if nargout(h)<6
  [hv,hs,hx,hsnext,hxnext]       = h(s,x,e,snext,xnext,params,output);
else
  [hv,hs,hx,hsnext,hxnext,hmult] = h(s,x,e,snext,xnext,params,output);
  hv      = hv.*hmult;
  if output.Js, hs     = hs.*hmult(:,:,ones(d,1)); end
  if output.Jx, hx     = hx.*hmult(:,:,ones(m,1)); end
  if output.Jsn, hsnext = hsnext.*hmult(:,:,ones(d,1)); end
  if output.Jxn, hxnext = hxnext.*hmult(:,:,ones(m,1)); end
end
hpval = hv;
hps   = hs;
hpsn  = hsnext;
p     = size(hv,2);
if output.Jx,  hpx = cat(3,hx,zeros(n,p,nx(1)+nx(2))); end
if output.Jxn, hpxn = cat(3,hxnext,zeros(n,p,nx(1)+nx(2))); end

return