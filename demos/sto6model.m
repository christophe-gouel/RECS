function [out1,out2,out3,out4,out5,out6] = sto6model(flag,s,x,z,e,snext,xnext,params,output)
% Equations of a two-country storage-trade model with supply reaction

% Copyright (C) 2011-2012 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

voidcell                        = cell(1,6);
[out1,out2,out3,out4,out5,out6] = voidcell{:};

[delta,r,k,alpha,theta,lambda,eta,mu,taua,taub] = params{:};

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
    if output.F
      out1        = zeros(n,m);
      out1(:,iSa) = Pa+k-((1-delta)/(1+r))*z(:,1);
      out1(:,iSb) = Pb+k-((1-delta)/(1+r))*z(:,2);
      out1(:,iPa) = Aa-Sa-Xa+Xb-Pa.^alpha;
      out1(:,iPb) = Ab-Sb-Xb+Xa-lambda*Pb.^alpha;
      out1(:,iXa) = Pa-Pb+theta+taub;
      out1(:,iXb) = Pb-Pa+theta+taua;
      out1(:,iHa) = z(:,3)-Ha.^mu;
      out1(:,iHb) = z(:,4)-(eta^(-mu))*Hb.^mu;
    end
    
    % df/ds
    if output.Js
      out2 = zeros(n,m,d);
      out2(:,iPa,iAa) = ones(n,1,1);
      out2(:,iPb,iAb) = ones(n,1,1);
    end
    
    % df/dx
    if output.Jx
      out3            = zeros(n,m,m);
      out3(:,iSa,iPa) = ones(n,1,1);
      out3(:,iSb,iPb) = ones(n,1,1);
      out3(:,iPa,iSa) = -ones(n,1,1);
      out3(:,iPa,iXa) = -ones(n,1,1);
      out3(:,iPa,iXb) = ones(n,1,1);
      out3(:,iPa,iPa) = -alpha*Pa.^(alpha-1);
      out3(:,iPb,iSb) = -ones(n,1,1);
      out3(:,iPb,iXb) = -ones(n,1,1);
      out3(:,iPb,iXa) = ones(n,1,1);
      out3(:,iPb,iPb) = -lambda*alpha*Pb.^(alpha-1);
      out3(:,iXa,iPa) = ones(n,1,1);
      out3(:,iXa,iPb) = -ones(n,1,1);
      out3(:,iXb,iPb) = ones(n,1,1);
      out3(:,iXb,iPa) = -ones(n,1,1);
      out3(:,iHa,iHa) = -mu*Ha.^(mu-1);
      out3(:,iHb,iHb) = -mu*(eta^(-mu))*Hb.^(mu-1);
    end
    
    % df/dz
    if output.Jz
      out4          = zeros(n,m,p);
      out4(:,iSa,1) = -((1-delta)/(1+r))*ones(n,1,1);
      out4(:,iSb,2) = -((1-delta)/(1+r))*ones(n,1,1);
      out4(:,iHa,3) = ones(n,1,1);
      out4(:,iHb,4) = ones(n,1,1);
    end
    
  case 'g'         % State transition
    if output.F
      out1         = zeros(n,d);
      out1(:,iAa)  = (1-delta)*Sa+Ha.*e(:,1);
      out1(:,iAb)  = (1-delta)*Sb+Hb.*e(:,2);
    end
    
    % dg/ds
    if output.Js, out2 = zeros(n,d,d); end
    
    % dg/dx
    if output.Jx
      out3            = zeros(n,d,m);
      out3(:,iAa,iSa) = (1-delta)*ones(n,1,1);
      out3(:,iAb,iSb) = (1-delta)*ones(n,1,1);
      out3(:,iAa,iHa) = e(:,1);
      out3(:,iAb,iHb) = e(:,2);
    end
    
  case 'h'     % Expectations
    n = size(snext,1);
    
    if output.F
      out1      = zeros(n,p);
      out1(:,1) = xnext(:,iPa);
      out1(:,2) = xnext(:,iPb);
      out1(:,3) = xnext(:,iPa);
      out1(:,4) = xnext(:,iPb);
    end
    
    % dh/ds
    if output.Js,  out2 = zeros(n,p,d); end
    
    % dh/dx
    if output.Jx,  out3 = zeros(n,p,m); end
    
    % dh/dsnext
    if output.Jsn, out4 = zeros(n,p,d); end
    
    % dh/dxnext
    if output.Jxn
      out5          = zeros(n,p,m);
      out5(:,1,iPa) = ones(n,1,1);
      out5(:,2,iPb) = ones(n,1,1);
      out5(:,3,iPa) = ones(n,1,1);
      out5(:,4,iPb) = ones(n,1,1);
    end
    
    % hmult
    if output.hmult
      out6      = ones(n,p);
      out6(:,3) = e(:,1);
      out6(:,4) = e(:,2);
    end
    
  case 'b'   % Bounds
    out1 = [zeros(n,2) -inf(n,2) zeros(n,2) -inf(n,2)];
    out2 = inf(n,8);
    
end
