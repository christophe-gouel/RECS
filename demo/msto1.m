function [out1,out2,out3,out4,out5] = mdsem03(flag,s,x,z,e,snext,xnext,alpha,k,delta,r,mu)
% Equations of a competitive storage model with supply reaction with an explicit
% market equilibrium equation

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

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
  out1(:,iS) = P+k-((1-delta)/(1+r))*z(:,1);
  out1(:,iH) = z(:,2)-H.^mu;
  out1(:,iP) = A-P.^alpha-S;

  % df/dx
  if nargout>=2
    out2          = zeros(n,m,m);
    out2(:,iS,iP) = ones(n,1,1);
    out2(:,iH,iH) = -mu*H.^(mu-1);
    out2(:,iP,iS) = -ones(n,1,1);
    out2(:,iP,iP) = -alpha*P.^(alpha-1);
  end

  % df/dz
  if nargout>=3
    out3         = zeros(n,m,p);
    out3(:,iS,1) = -((1-delta)/(1+r))*ones(n,1,1);
    out3(:,iH,2) = ones(n,1,1);
  end
 case 'g' % STATE TRANSITION FUNCTION
  % g
  out1 = (1-delta)*S + H.*e;

  % dg/dx
  if nargout>=2
    out2         = zeros(n,d,m);
    out2(:,1,iS) = (1-delta)*ones(n,1,1);
    out2(:,1,iH) = e;
  end
 case 'h' % EXPECTATION FUNCTION
  % h
  out1      = zeros(size(snext,1),p);
  out1(:,1) = xnext(:,iP);
  out1(:,2) = xnext(:,iP);

  % dh/dx
  if nargout>=2, out2 = zeros(size(snext,1),p,m); end

  % dh/dsnext
  if nargout>=3, out3 = zeros(size(snext,1),p,d); end

  % dh/dxnext
  if nargout>=4
    out4         = zeros(size(snext,1),p,m);
    out4(:,1,iP) = ones(size(snext,1),1,1);
    out4(:,2,iP) = ones(size(snext,1),1,1);
  end

  % hmult
  if nargout==5
    out5      = ones(size(snext,1),p);
    out5(:,2) = e;
  end
 case 'b' % BOUND FUNCTION
  out1 = [zeros(n,1) -inf(n,2)];
  out2 = inf(n,3);
 case 'e' % EULER EQUATION ERROR
  out1 = zeros(n,2);
  out1(:,iS) = ones(n,1)-(max(A.^(1/alpha),((1-delta)/(1+r))*z(:,1)-k).^alpha)./(A-S);
  out1(:,iH) = 1-(z(:,2).^(1/mu))./H;
end