function [out1,out2,out3,out4] = mdsem13(flag,s,x,z,e,snext,xnext,params)
% Equations of a two-country storage-trade model with supply reaction

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

delta = params(1);
r = params(2);
k = params(3);
alpha = params(4);
theta = params(5);
lambda = params(6);
eta = params(7);
mu = params(8);
taua = params(9);
taub = params(10);

n = size(s,1);
d = 2;
m = 8;
p = 4;

% s = [Aa Ab]
iAa   = 1;
iAb   = 2;
if ~isempty(s)
  Aa  = s(:,iAa);
  Ab  = s(:,iAb);
end

% x = [Sa Sb Pa Pb Xa Xb Ha Hb]
iSa = 1;
iSb = 2;
iPa = 3;
iPb = 4;
iXa = 5;
iXb = 6;
iHa = 7;
iHb = 8;
if ~isempty(x)
  Sa = x(:,iSa);
  Sb = x(:,iSb);
  Pa = x(:,iPa);
  Pb = x(:,iPb);
  Xa = x(:,iXa);
  Xb = x(:,iXb);
  Ha = x(:,iHa);
  Hb = x(:,iHb);
end

switch flag
 case 'f'     % Equilibrium
  out1(:,iSa) = Pa+k-((1-delta)/(1+r))*z(:,1);
  out1(:,iSb) = Pb+k-((1-delta)/(1+r))*z(:,2);
  out1(:,iPa) = Aa-Sa-Xa+Xb-Pa.^alpha;
  out1(:,iPb) = Ab-Sb-Xb+Xa-lambda*Pb.^alpha;
  out1(:,iXa) = Pa-Pb+theta+taub;
  out1(:,iXb) = Pb-Pa+theta+taua;
  out1(:,iHa) = z(:,3)-Ha.^mu;
  out1(:,iHb) = z(:,4)-(eta^(-mu))*Hb.^mu;

  % fx
  if nargout>=2,
    out2          = zeros(n,m,m);
    out2(:,iSa,iPa) = ones(n,1,1);
    out2(:,iSb,iPb) = ones(n,1,1);
    out2(:,iPa,iSa) = -ones(n,1,1);
    out2(:,iPa,iXa) = -ones(n,1,1);
    out2(:,iPa,iXb) = ones(n,1,1);
    out2(:,iPa,iPa) = -alpha*Pa.^(alpha-1);
    out2(:,iPb,iSb) = -ones(n,1,1);
    out2(:,iPb,iXb) = -ones(n,1,1);
    out2(:,iPb,iXa) = ones(n,1,1);
    out2(:,iPb,iPb) = -lambda*alpha*Pb.^(alpha-1);
    out2(:,iXa,iPa) = ones(n,1,1);
    out2(:,iXa,iPb) = -ones(n,1,1);
    out2(:,iXb,iPb) = ones(n,1,1);
    out2(:,iXb,iPa) = -ones(n,1,1);
    out2(:,iHa,iHa) = -mu*Ha.^(mu-1);
    out2(:,iHb,iHb) = -mu*(eta^(-mu))*Hb.^(mu-1);
  end

  % fz
  if nargout>=3,
    out3          = zeros(n,m,p);
    out3(:,iSa,1)  = -((1-delta)/(1+r))*ones(n,1,1);
    out3(:,iSb,2)  = -((1-delta)/(1+r))*ones(n,1,1);
    out3(:,iHa,3) = ones(n,1,1);
    out3(:,iHb,4) = ones(n,1,1);
  end


 case 'g'         % State transition
  % g
  out1        = zeros(n,d);
  out1(:,iAa)  = (1-delta)*Sa+Ha.*e(:,1);
  out1(:,iAb)  = (1-delta)*Sb+Hb.*e(:,2);

  % gx
  if nargout>=2
    out2          = zeros(n,d,m);
    out2(:,iAa,iSa) = (1-delta)*ones(n,1,1);
    out2(:,iAb,iSb) = (1-delta)*ones(n,1,1);
    out2(:,iAa,iHa) = e(:,1);
    out2(:,iAb,iHb) = e(:,2);
  end


 case 'h'     % Expectations
  % h
  out1(:,1) = xnext(:,iPa);
  out1(:,2) = xnext(:,iPb);
  out1(:,3) = e(:,1).*xnext(:,iPa);
  out1(:,4) = e(:,2).*xnext(:,iPb);

  % dh/dx
  if nargout>=2, out2 = zeros(size(snext,1),p,m); end

  % dh/dsnext
  if nargout>=3, out3 = zeros(size(snext,1),p,d); end

  % dh/dxnext
  if nargout==4
    out4         = zeros(size(snext,1),p,m);
    out4(:,1,iPa) = ones(size(snext,1),1,1);
    out4(:,2,iPb) = ones(size(snext,1),1,1);
    out4(:,3,iPa) = e(:,1);
    out4(:,4,iPb) = e(:,2);
  end


 case 'b'   % Bounds
  out1 = [zeros(n,2) -inf(n,2) zeros(n,2) -inf(n,2)];
  out2 = inf(n,8);

end
