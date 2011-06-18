function [out1,out2,out3,out4,out5] = gro1model(flag,s,x,z,e,snext,xnext,params,output)
% GRO1MODEL Equations of the stochastic growth model

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

voidcell                   = cell(1,5);
[out1,out2,out3,out4,out5] = voidcell{:};

[tau,beta,alpha,delta,rho,a,theta] = params{:};

n = size(s,1);
d = 2;
m = 2;
p = 1;

iK = 1;
iZ = 2;
if ~isempty(s)
  K = s(:,iK);
  Z = s(:,iZ);
end

iC = 1;
iL = 2;
if ~isempty(x)
  C = x(:,iC);
  L = x(:,iL);
end

switch flag
 case 'f' % EQUILIBRIUM FUNCTION
  % f
  if output.F
    out1       = zeros(n,m);
    out1(:,iC) = ((C.^theta.*(1-L).^(1-theta)).^(1-tau))./C-beta*z;
    out1(:,iL) = (1-theta)/(1-L)-theta*(1-alpha)*exp(Z).*K.^alpha.*L.^(1-alpha)./C;
  end

  % df/ds
  if output.Js
    out2          = zeros(n,m,d); 
    out2(:,iL,iK) = -theta*(1-alpha)*alpha*exp(Z).*K.^(alpha-1).*L.^(1-alpha)./C;
    out2(:,iL,iZ) = -theta*(1-alpha)*exp(Z).*K.^alpha.*L.^(1-alpha)./C;
  end

  % df/dx
  if output.Jx
    out3          = zeros(n,m,m);
    out3(:,iC,iC) = -tau*theta*((C.^theta.*(1-L).^(1-theta)).^(1-tau))./(C.^2); 
    out2(:,iC,iL) = -(1-tau)*(1-theta)*((C.^theta.*(1-L).^(1-theta)).^(1-tau))./(C.*(1-L));
  end

  % df/dz
  if output.Jz
    out4         = zeros(n,m,p);
    out4(:,iC,1) = -beta*ones(n,1,1);
  end
  
 case 'g' % STATE TRANSITION FUNCTION
  % g
  if output.F
    out1       = zeros(n,d);
% $$$     disp('Investment:')
% $$$     a*exp(Z).*K.^alpha-C
% $$$     disp('Residual capital:')
% $$$     (1-delta)*K
    out1(:,iK) = a*exp(Z).*K.^alpha.*L.^(1-alpha)+(1-delta)*K-C;
    out1(:,iZ) = rho*Z+e;
  end
  
  % dg/ds
  if output.Js
    out2          = zeros(n,d,d); 
    out2(:,iK,iK) = a*alpha*exp(Z).*K.^(alpha-1).*L.^(1-alpha)+(1-delta)*ones(n,1,1);
    out2(:,iK,iZ) = a*exp(Z).*K.^alpha;
    out2(:,iZ,iZ) = rho*ones(n,1,1);
  end

  % dg/dx
  if output.Jx
    out3          = zeros(n,d,m);
    out3(:,iK,iC) = -ones(n,1,1);
    out3(:,iK,iL) = a*(1-alpha)*exp(Z).*K.^alpha.*L.^(-alpha)
  end
  
 case 'h' % EXPECTATION FUNCTION
  n = size(snext,1);
  % h
  if output.F
    out1 = ((xnext(:,iC).^theta.*(1-xnext(:,iL)).^(1-theta)).^(1-tau))./xnext(:,iC).*...
           (1-delta+a*alpha*exp(snext(:,iZ)).*snext(:,iK).^(alpha-1).*xnext(:,iL).^(1-alpha)); 
  end

  % dh/ds
  if output.Js,  out2 = zeros(n,p,d); end

  % dh/dx
  if output.Jx,  out3 = zeros(n,p,m); end

  % dh/dsnext
  if output.Jsn
    out4         = zeros(n,p,d);
    out4(:,1,iK) = ((xnext(:,iC).^theta.*(1-xnext(:,iL)).^(1-theta)).^(1-tau))./xnext(:,iC).*...
        a*alpha*(alpha-1).*exp(snext(:,iZ)).*...
        snext(:,iK).^(alpha-2).*xnext(:,iL).^(1-alpha); 
    out4(:,1,iZ) = ((xnext(:,iC).^theta.*(1-xnext(:,iL)).^(1-theta)).^(1-tau))./xnext(:,iC).*...
        a*alpha.*exp(snext(:,iZ)).*...
        snext(:,iK).^(alpha-1).*xnext(:,iL).^(1-alpha); 
  end

  % dh/dxnext
  if output.Jxn
    out5         = zeros(n,p,m);
    out5(:,1,iC) = -tau*theta*((xnext(:,iC).^theta.*(1-xnext(:,iL)).^(1-theta)).^(1-tau))./(xnext(:,iC).^2).*...
        (1-delta+a*alpha*exp(snext(:,iZ)).*snext(:,iK).^(alpha-1).*xnext(:,iL).^(1-alpha));
    out5(:,1,iL) = ((xnext(:,iC).^theta.*(1-xnext(:,iL)).^(1-theta)).^(1-tau))./xnext(:,iC).*...
        a*(1-alpha)*alpha*exp(snext(:,iZ)).*snext(:,iK).^(alpha-1).*xnext(:,iL).^(-alpha)+...
        -(1-tau)*(1-theta)*((xnext(:,iC).^theta.*(1-xnext(:,iL)).^(1-theta)).^(1-tau))./(xnext(:,iC).*(1-xnext(:,iL))).*...
        (1-delta+a*alpha*exp(snext(:,iZ)).*snext(:,iK).^(alpha-1).*xnext(:,iL).^(1-alpha)); 
  end
  
 case 'b' % BOUND FUNCTION
  out1 = -inf(n,2) ;
  out2 =  inf(n,2);

 case 'ss' % Deterministic steady-state
  phi   = ((1/beta-1+delta)/alpha)^(1/(1-alpha));
  Omega = phi^(1-alpha)-delta;
  Psi   = (1-alpha)*phi^(-alpha)*theta/(1-theta);
  
% $$$   Kstar = ((1/beta-1+delta)/(a*alpha))^(1/(alpha-1));
  Kstar = Psi/(Omega+phi*Psi);
  %%  'State variables'
  out1 = [Kstar 0];
  %%  'Response variables'
  out2 = [Omega*Kstar phi*Kstar];
  %%  'Expectations'
  out3 = out2(:,iC).^(-tau)*(1-delta+a*alpha*Kstar^(alpha-1).*out2(:,iL).^(1-alpha));
end
