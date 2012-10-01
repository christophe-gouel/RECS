function [out1,out2,out3,out4,out5,out6] = sto2model(flag,s,x,z,e,snext,xnext,params,output)
% STO2MODEL Equations of a competitive storage model with floor-price backed by public storage

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

voidcell                        = cell(1,6);
[out1,out2,out3,out4,out5,out6] = voidcell{:};

[k,delta,r,mu,alpha,PF,Sgbar] = params{:};

n = size(s,1);
d = 1;
m = 4;
p = 2;

A = s;

% x = [S H P Sg]
iS  = 1; % Speculative storage
iH  = 2;
iP  = 3;
iSg = 4; % Public storage
if ~isempty(x)
  S  = x(:,iS);
  H  = x(:,iH);
  P  = x(:,iP);
  Sg = x(:,iSg);
end

switch flag
 case 'f' % EQUILIBRIUM FUNCTION
  % f
  if output.F
    out1        = zeros(n,m);
    out1(:,iS)  = P+k-((1-delta)/(1+r))*z(:,1);
    out1(:,iH)  = z(:,2)-H.^mu;
    out1(:,iP)  = A-P.^alpha-S-Sg;
    out1(:,iSg) = P-PF;
  end
  
  % df/ds
  if output.Js
    out2         = zeros(n,m,d);
    out2(:,iP,1) = ones(n,1,1);
  end

  % df/dx
  if output.Jx
    out3           = zeros(n,m,m);
    out3(:,iS,iP)  = ones(n,1,1);
    out3(:,iH,iH)  = -mu*H.^(mu-1);
    out3(:,iP,iS)  = -ones(n,1,1);
    out3(:,iP,iP)  = -alpha*P.^(alpha-1);
    out3(:,iP,iSg) = -ones(n,1,1);
    out3(:,iSg,iP) = ones(n,1,1);
  end

  % df/dz
  if output.Jz
    out4         = zeros(n,m,p);
    out4(:,iS,1) = -((1-delta)/(1+r))*ones(n,1,1);
    out4(:,iH,2) = ones(n,1,1);
  end
  
 case 'g' % STATE TRANSITION FUNCTION
  % g
  if output.F, out1 = (1-delta)*(S+Sg) + H.*e; end
  
  % dg/ds
  if output.Js, out2 = zeros(n,d,d); end

  % dg/dx
  if output.Jx
    out3          = zeros(n,d,m);
    out3(:,1,iS)  = (1-delta)*ones(n,1,1);
    out3(:,1,iH)  = e;
    out3(:,1,iSg) = (1-delta)*ones(n,1,1);
  end
  
 case 'h' % EXPECTATION FUNCTION
  n = size(snext,1);
  % h
  if output.F
    out1      = zeros(n,p);
    out1(:,1) = xnext(:,iP);
    out1(:,2) = xnext(:,iP);
  end

  % dh/ds
  if output.Js, out2 = zeros(n,p,d); end

  % dh/dx
  if output.Jx, out3 = zeros(n,p,m); end

  % dh/dsnext
  if output.Jsn, out4 = zeros(n,p,d); end

  % dh/dxnext
  if output.Jxn
    out5         = zeros(n,p,m);
    out5(:,1,iP) = ones(n,1,1);
    out5(:,2,iP) = ones(n,1,1);
  end

  % hmult
  if output.hmult
    out6      = ones(n,p);
    out6(:,2) = e;
  end

 case 'b' % BOUND FUNCTION
  out1 = [zeros(n,1) -inf(n,2) zeros(n,1)];
  out2 = [inf(n,3) Sgbar*ones(n,1)];
end
