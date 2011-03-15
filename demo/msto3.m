function [out1,out2,out3,out4] = mdsem05(flag,s,x,z,e,snext,xnext,params)

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

[delta,r,alpha,k,mu,PF,PC,Sgbar] = params{:};

n = size(s,1);
d = 2;
m = 5;
p = 2;

% s = [A Sg]
iA  = 1;
iSg = 2;
if ~isempty(s)
  A  = s(:,iA);
  Sg = s(:,iSg);
end

% x = [S H P dSgp dSgm]
iS    = 1;
iH    = 2;
iP    = 3;
idSgp = 4;
idSgm = 5;
if ~isempty(x)
  S    = x(:,iS);
  H    = x(:,iH);
  P    = x(:,iP);
  dSgp = x(:,idSgp);
  dSgm = x(:,idSgm);
end

switch flag
 case 'f' % EQUILIBRIUM FUNCTION
  % f
  out1          = zeros(n,m);
  out1(:,iS)    = P+k-((1-delta)/(1+r))*z(:,1);
  out1(:,iH)    = z(:,2)-H.^mu;
  out1(:,iP)    = A-P.^alpha-S-dSgp+dSgm;
  out1(:,idSgp) = P-PF;
  out1(:,idSgm) = PC-P;

  % df/dx
  if nargout>=2
    out2             = zeros(n,m,m);
    out2(:,iS,iP)    = ones(n,1,1);
    out2(:,iH,iH)    = -mu*H.^(mu-1);
    out2(:,iP,iS)    = -ones(n,1,1);
    out2(:,iP,iP)    = -alpha*P.^(alpha-1);
    out2(:,iP,idSgp) = -ones(n,1,1);
    out2(:,iP,idSgm) = ones(n,1,1);
    out2(:,idSgp,iP) = ones(n,1,1);
    out2(:,idSgm,iP) = -ones(n,1,1);
  end

  % df/dz
  if nargout>=3
    out3         = zeros(n,m,p);
    out3(:,iS,1) = -((1-delta)/(1+r))*ones(n,1,1);
    out3(:,iH,2) = ones(n,1,1);
  end
 case 'g' % STATE TRANSITION FUNCTION
  % g
  out1        = zeros(n,d);
  out1(:,iA)  = (1-delta)*S + H.*e;
  out1(:,iSg) = (1-delta)*Sg+dSgp-dSgm;

  % dg/dx
  if nargout>=2
    out2              = zeros(n,d,m);
    out2(:,iA,iS)     = (1-delta)*ones(n,1,1);
    out2(:,iA,iH)     = e;
    out2(:,iSg,idSgp) = ones(n,1,1);
    out2(:,iSg,idSgm) = -ones(n,1,1);
  end
 case 'h' % EXPECTATION FUNCTION
  % h
  out1      = zeros(size(snext,1),p);
  out1(:,1) = xnext(:,iP);
  out1(:,2) = e.*xnext(:,iP);

  % dh/dx
  if nargout>=2, out2 = zeros(size(snext,1),p,m); end

  % dh/dsnext
  if nargout>=3, out3 = zeros(size(snext,1),p,d); end

  % dh/dxnext
  if nargout==4
    out4         = zeros(size(snext,1),p,m);
    out4(:,1,iP) = ones(size(snext,1),1,1);
    out4(:,2,iP) = e;
  end
 case 'b' % BOUND FUNCTION
  out1 = [zeros(n,1) -inf(n,2) zeros(n,2)];
  out2 = [inf(n,3) (Sgbar-(1-delta)*Sg) (1-delta)*Sg];
end
