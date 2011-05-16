function [out1,out2,out3,out4,out5] = sto3model(flag,s,x,z,e,snext,xnext,params,output)

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

voidcell                        = cell(1,6);
[out1,out2,out3,out4,out5,out6] = voidcell{:};

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
  if output.F
    out1          = zeros(n,m);
    out1(:,iS)    = P+k-((1-delta)/(1+r))*z(:,1);
    out1(:,iH)    = z(:,2)-H.^mu;
    out1(:,iP)    = A-P.^alpha-S-dSgp+dSgm;
    out1(:,idSgp) = P-PF;
    out1(:,idSgm) = PC-P;
  end
  
  % df/ds
  if output.Js
    out2         = zeros(n,m,d);
    out2(:,iP,1) = ones(n,1,1);
  end

  % df/dx
  if output.Jx
    out3             = zeros(n,m,m);
    out3(:,iS,iP)    = ones(n,1,1);
    out3(:,iH,iH)    = -mu*H.^(mu-1);
    out3(:,iP,iS)    = -ones(n,1,1);
    out3(:,iP,iP)    = -alpha*P.^(alpha-1);
    out3(:,iP,idSgp) = -ones(n,1,1);
    out3(:,iP,idSgm) = ones(n,1,1);
    out3(:,idSgp,iP) = ones(n,1,1);
    out3(:,idSgm,iP) = -ones(n,1,1);
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
    out1        = zeros(n,d);
    out1(:,iA)  = (1-delta)*S + H.*e;
    out1(:,iSg) = (1-delta)*Sg+dSgp-dSgm;
  end
  
  % dg/ds
  if output.Js
    out2 = zeros(n,d,d);
    out2(:,iSg,iSg) = (1-delta)*ones(n,1,1);
  end

  % dg/dx
  if output.Jx
    out3              = zeros(n,d,m);
    out3(:,iA,iS)     = (1-delta)*ones(n,1,1);
    out3(:,iA,iH)     = e;
    out3(:,iSg,idSgp) = ones(n,1,1);
    out3(:,iSg,idSgm) = -ones(n,1,1);
  end
  
 case 'h' % EXPECTATION FUNCTION
  n = size(snext,1);
  % h
  if output.F
    out1      = zeros(n,p);
    out1(:,1) =    xnext(:,iP);
    out1(:,2) = e.*xnext(:,iP);
  end

  % dh/ds
  if output.Js,  out2 = zeros(n,p,d); end

  % dh/dx
  if output.Jx,  out3 = zeros(n,p,m); end

  % dh/dsnext
  if output.Jsn, out4 = zeros(n,p,d); end

  % dh/dxnext
  if output.Jxn
    out5         = zeros(n,p,m);
    out5(:,1,iP) =  ones(n,1,1);
    out5(:,2,iP) = e;
  end
  
 case 'b' % BOUND FUNCTION
  out1 = [zeros(n,1) -inf(n,2) zeros(n,2)];
  out2 = [inf(n,3) (Sgbar-(1-delta)*Sg) (1-delta)*Sg];
end
