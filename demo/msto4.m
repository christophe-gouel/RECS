function [out1,out2,out3,out4] = msto4(flag,s,x,z,e,snext,xnext,delta,r,k,alpha,tau,rho,sigma)
% Equations of a storage-trade model of one small-country

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

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
  out1(:,iS) = P+k-((1-delta)/(1+r))*z;
  out1(:,iP) = A+M-S-X-P.^alpha;
  out1(:,iM) = Pw+tau-P;
  out1(:,iX) = tau+P-Pw;

  % fx
  if nargout>=2,
    out2          = zeros(n,m,m);
    out2(:,iS,iP) = ones(n,1,1);
    out2(:,iP,iS) = -ones(n,1,1);
    out2(:,iP,iP) = -alpha*P.^(alpha-1);
    out2(:,iP,iM) = ones(n,1,1);
    out2(:,iP,iX) = -ones(n,1,1);
    out2(:,iM,iP) = -ones(n,1,1);
    out2(:,iX,iP) = ones(n,1,1);
  end

  % fz
  if nargout>=3,
    out3          = zeros(n,m,p);
    out3(:,iS,1)  = -((1-delta)/(1+r))*ones(n,1,1);
  end

 case 'g'
  % g
  out1        = zeros(n,d);
  out1(:,iA)  = (1-delta)*S+e(:,1);
  out1(:,iPw) = Pw.^(rho).*exp(e(:,2))/exp(sigma^2/2);

  % gx
  if nargout>=2
    out2          = zeros(n,d,m);
    out2(:,iA,iS) = (1-delta)*ones(n,1,1);
  end
 case 'h'
  % h
  out1 = xnext(:,iP);

  % dh/dx
  if nargout>=2, out2 = zeros(size(snext,1),p,m); end

  % dh/dsnext
  if nargout>=3, out3 = zeros(size(snext,1),p,d); end

  % dh/dxnext
  if nargout==4
    out4         = zeros(size(snext,1),p,m);
    out4(:,1,iP) = ones(size(snext,1),1,1);
  end

 case 'b'
  out1 = [zeros(n,1) -inf(n,1) zeros(n,2)];
  out2 = inf(n,4);
end
