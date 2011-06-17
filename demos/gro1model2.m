function [out1,out2,out3,out4,out5] = gro1model2(flag,s,x,z,e,snext,xnext,params,output)
% GRO1MODEL Equations of the stochastic growth model

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

voidcell                   = cell(1,5);
[out1,out2,out3,out4,out5] = voidcell{:};

[tau,beta,alpha,delta,rho,a] = params{:};

n = size(s,1);
d = 2;
m = 1;
p = 1;

iK = 1;
iZ = 2;
if ~isempty(s)
  K = s(:,iK);
  Z = s(:,iZ);
end

C = x;

switch flag
 case 'f' % EQUILIBRIUM FUNCTION
  % f
  if output.F,  out1 = C.^(-1)-beta*z; end

  % df/ds
  if output.Js, out2 = zeros(n,m,d); end

  % df/dx
  if output.Jx, out3 = -C.^(-2); end

  % df/dz
  if output.Jz, out4 = -beta*ones(n,1,1); end
  
 case 'g' % STATE TRANSITION FUNCTION
  % g
  if output.F
    out1       = zeros(n,d);
    out1(:,iK) = a*exp(Z).*K.^alpha+(1-delta)*K-C;
    out1(:,iZ) = rho*Z+e;
  end
  
  % dg/ds
  if output.Js
    out2          = zeros(n,d,d); 
    out2(:,iK,iK) = a*alpha*exp(Z).*K.^(alpha-1)+(1-delta)*ones(n,1,1);
    out2(:,iK,iZ) = a*exp(Z).*K.^alpha;
    out2(:,iZ,iZ) = rho*ones(n,1,1);
  end

  % dg/dx
  if output.Jx
    out3         = zeros(n,d,m);
    out3(:,iK,1) = -ones(n,1,1);
  end
  
 case 'h' % EXPECTATION FUNCTION
  n = size(snext,1);
  % h
  if output.F
    out1 = xnext.^(-1).*(1-delta+a*alpha*exp(snext(:,iZ)).*snext(:,iK).^(alpha-1)); 
  end

  % dh/ds
  if output.Js,  out2 = zeros(n,p,d); end

  % dh/dx
  if output.Jx,  out3 = zeros(n,p,m); end

  % dh/dsnext
  if output.Jsn
    out4         = zeros(n,p,d);
    out4(:,1,iK) = a*alpha*(alpha-1)*xnext.^(-1).*exp(snext(:,iZ)).*...
        snext(:,iK).^(alpha-2); 
    out4(:,1,iZ) = a*alpha*xnext.^(-1).*exp(snext(:,iZ)).*snext(:,iK).^(alpha-1); 
  end

  % dh/dxnext
  if output.Jxn
    out5 = -xnext.^(-2).*(1-delta+a*alpha*exp(snext(:,iZ)).*snext(:,iK).^(alpha-1)); ; 
  end
  
 case 'b' % BOUND FUNCTION
  out1 = -inf(n,1);
  out2 =  inf(n,1);

 case 'ss' % Deterministic steady-state
  Kstar = ((1/beta-1+delta)/(a*alpha))^(1/(alpha-1));
  %%  'State variables'
  out1 = [Kstar 0];
  %%  'Response variables'
  out2 = Kstar^alpha-delta*Kstar;
  %%  'Expectations'
  out3 = out2^(-1)*(1-delta+a*alpha*Kstar^(alpha-1));
end
