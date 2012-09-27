function [out1,out2,out3,out4,out5] = sto5model(flag,s,x,z,e,snext,xnext,params,output)
% STO5MODEL Equations of a storage-trade model of one small-country

% Copyright (C) 2011-2012 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

voidcell                        = cell(1,6);
[out1,out2,out3,out4,out5,out6] = voidcell{:};

[delta,r,k,alpha,tau,rho,sigma] = params{:};

n = size(s,1);
d = 2;
m = 4;
p = 1;

% s = [A Pw]
iA   = 1;
iPw  = 2;
if ~isempty(s)
  A  = s(:,iA);
  Pw = s(:,iPw);
end

% x = [S P M X]
iS = 1;
iP = 2;
iM = 3;
iX = 4;
if ~isempty(x)
  S = x(:,iS);
  P = x(:,iP);
  M = x(:,iM);
  X = x(:,iX);
end

switch flag
 case 'f'
  if output.F
    out1       = zeros(n,m);
    out1(:,iS) = P+k-((1-delta)/(1+r))*z;
    out1(:,iP) = A+M-S-X-P.^alpha;
    out1(:,iM) = Pw+tau-P;
    out1(:,iX) = tau+P-Pw;
  end
  
  % df/ds
  if output.Js
    out2           = zeros(n,m,d);
    out2(:,iP,iA)  =  ones(n,1,1);
    out2(:,iM,iPw) =  ones(n,1,1);
    out2(:,iX,iPw) = -ones(n,1,1);
  end

  % fx
  if output.Jx
    out3          = zeros(n,m,m);
    out3(:,iS,iP) = ones(n,1,1);
    out3(:,iP,iS) = -ones(n,1,1);
    out3(:,iP,iP) = -alpha*P.^(alpha-1);
    out3(:,iP,iM) = ones(n,1,1);
    out3(:,iP,iX) = -ones(n,1,1);
    out3(:,iM,iP) = -ones(n,1,1);
    out3(:,iX,iP) = ones(n,1,1);
  end

  % fz
  if output.Jz
    out4          = zeros(n,m,p);
    out4(:,iS,1)  = -((1-delta)/(1+r))*ones(n,1,1);
  end

 case 'g'
  % g
  if output.F
    out1        = zeros(n,d);
    out1(:,iA)  = (1-delta)*S+e(:,1);
    out1(:,iPw) = Pw.^(rho).*exp(e(:,2))/exp(sigma^2/2);
  end

  % dg/ds
  if output.Js
    out2            = zeros(n,d,d); 
    out2(:,iPw,iPw) = rho*Pw.^(rho-1).*exp(e(:,2))/exp(sigma^2/2);
  end

  % gx
  if output.Jx
    out3          = zeros(n,d,m);
    out3(:,iA,iS) = (1-delta)*ones(n,1,1);
  end
 case 'h'
  n = size(snext,1);
  % h
  if output.F, out1 = xnext(:,iP); end

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
  end

 case 'b'
  out1 = [zeros(n,1) -inf(n,1) zeros(n,2)];
  out2 = inf(n,4);
end
