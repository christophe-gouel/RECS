function [out1,out2,out3,out4,out5,out6] = mtest(flag,s,x,z,e,snext,xnext,params,output)

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

voidcell                       = cell(1,6);
[ou1,out2,out3,out4,out5,out6] = voidcell{:};

[alpha,k,delta,r,mu] = params{:};

n = size(s,1);
d = 1;
m = 3;
p = 2;

A = s;

% x = [S H P]
iS = 1;
iH = 2;
iP = 3;
if ~isempty(x)
  S  = x(:,iS);
  H  = x(:,iH);
  P  = x(:,iP);
end

switch flag
 case 'f' % EQUILIBRIUM FUNCTION
  % f
  if output.F
    out1(:,iS) = P+k-((1-delta)/(1+r))*z(:,1);
    out1(:,iH) = z(:,2)-H.^mu;
    out1(:,iP) = A-P.^alpha-S;
  end

  % df/ds
  if output.Js
    out2 = zeros(n,m,d);
    out2(:,iP,1) = ones(n,1,1);
  end

  % df/dx
  if output.Jx
    out3          = zeros(n,m,m);
    out3(:,iS,iP) = ones(n,1,1);
    out3(:,iH,iH) = -mu*H.^(mu-1);
    out3(:,iP,iS) = -ones(n,1,1);
    out3(:,iP,iP) = -alpha*P.^(alpha-1);
  end

  % df/dz
  if output.Jz
    out4         = zeros(n,m,p);
    out4(:,iS,1) = -((1-delta)/(1+r))*ones(n,1,1);
    out4(:,iH,2) = ones(n,1,1);
  end
 case 'g' % STATE TRANSITION FUNCTION
  % g
  if output.F
    out1 = (1-delta)*S + H.*e;
  end

  % dg/ds
  if output.Js, out2 = zeros(n,d,d); end

  % dg/dx
  if output.Jx
    out3         = zeros(n,d,m);
    out3(:,1,iS) = (1-delta)*ones(n,1,1);
    out3(:,1,iH) = e;
  end
 case 'h' % EXPECTATION FUNCTION
  % h
  if output.F
    out1      = zeros(size(snext,1),p);
    out1(:,1) =    xnext(:,iP);
    out1(:,2) = e.*xnext(:,iP);
  end

  % dh/ds
  if output.Js, out2 = zeros(size(snext,1),p,d); end

  % dh/dx
  if output.Jx, out3 = zeros(size(snext,1),p,m); end

  % dh/dsnext
  if output.Jsn, out4 = zeros(size(snext,1),p,d); end

  % dh/dxnext
  if output.Jxn
    out5         = zeros(size(snext,1),p,m);
    out5(:,1,iP) = ones(size(snext,1),1,1);
    out5(:,2,iP) = e;
  end

 case 'b' % BOUND FUNCTION
  out1 = [zeros(n,1) -inf(n,2)];
  out2 = inf(n,3);
end