function [out1,out2,out3,out4,out5] = gro2model(flag,s,x,z,e,snext,xnext,p,output)

voidcell                   = cell(1,5);
[out1,out2,out3,out4,out5] = voidcell{:};

switch flag

  case 'b'
    n = size(s,1);
    out1 = -inf(n,3);
    out1(:,2) = 0;
    out2 = inf(n,3);

  case 'f'
    n = size(s,1);

    % f
    if output.F
      out1 = zeros(n,3);
      out1(:,1) = -p(1).*s(:,1).^p(6).*exp(s(:,2)) + x(:,1) + x(:,2);
      out1(:,2) = x(:,3) + x(:,1).^(-p(2));
      out1(:,3) = p(4).*z(:,1) + x(:,3);
    end

    % df/ds
    if output.Js
      out2 = zeros(n,3,2);
      out2(:,1,1) = -p(1).*p(6).*s(:,1).^p(6).*exp(s(:,2))./s(:,1); % d eq_1 w.r.t. K
      out2(:,1,2) = -p(1).*s(:,1).^p(6).*exp(s(:,2)); % d eq_1 w.r.t. Z
    end

    % df/dx
    if output.Jx
      out3 = zeros(n,3,3);
      out3(:,1,1) = 1; % d eq_1 w.r.t. C
      out3(:,1,2) = 1; % d eq_1 w.r.t. I
      out3(:,2,1) = -p(2).*x(:,1).^(-p(2))./x(:,1); % d eq_2 w.r.t. C
      out3(:,2,3) = 1; % d eq_2 w.r.t. Mu
      out3(:,3,3) = 1; % d eq_3 w.r.t. Mu
    end

    % df/dz
    if output.Jz
      out4 = zeros(n,3,1);
      out4(:,3,1) = p(4); % d eq_3 w.r.t. E
    end
  
  case 'g'
    n = size(s,1);

    % g
    if output.F
      out1 = zeros(n,2);
      out1(:,1) = x(:,2) + s(:,1).*(-p(3) + 1);
      out1(:,2) = p(5).*s(:,2) + e(:,1);      
    end

    if output.Js
      out2 = zeros(n,2,2);
      out2(:,1,1) = -p(3) + 1; % d eq_1 w.r.t. K
      out2(:,2,2) = p(5); % d eq_2 w.r.t. Z
    end

    if output.Jx
      out3 = zeros(n,2,3);
      out3(:,1,2) = 1; % d eq_1 w.r.t. I
    end
  
  case 'h'
    n = size(snext,1);

    %h
    if output.F
      out1 = zeros(n,1);
      out1(:,1) = p(1).*p(6).*xnext(:,1).^(-p(2)).*snext(:,1).^(p(6) - 1).*exp(snext(:,2)) - xnext(:,3).*(-p(3) + 1);
    end

    if output.Js
      out2 = zeros(n,1,2);
    end

    if output.Jx
      out3 = zeros(n,1,3);
    end

    if output.Jsn
      out4 = zeros(n,1,2);
      out4(:,1,1) = p(1).*p(6).*xnext(:,1).^(-p(2)).*snext(:,1).^(p(6) - 1).*(p(6) - 1).*exp(snext(:,2))./snext(:,1); % d eq_1 w.r.t. K(1)
      out4(:,1,2) = p(1).*p(6).*xnext(:,1).^(-p(2)).*snext(:,1).^(p(6) - 1).*exp(snext(:,2)); % d eq_1 w.r.t. Z(1)
    end

    if output.Jxn
      out5 = zeros(n,1,3);
      out5(:,1,1) = -p(1).*p(6).*p(2).*xnext(:,1).^(-p(2)).*snext(:,1).^(p(6) - 1).*exp(snext(:,2))./xnext(:,1); % d eq_1 w.r.t. C(1)
      out5(:,1,3) = p(3) - 1; % d eq_1 w.r.t. Mu(1)
    end
  
  case 'e'
    out1 = [];

  case 'params'
    out1 = [ 0.24077193  1.          0.0196      0.95        0.9         0.3       ];

  case 'ss'
    out1 = [1.0  0.0];
    out2 = [0.221171929825  0.0  0.0];

  case 'J'
    out1.fs = [1 1; 0 0; 0 0];
    out1.fx = [1 1 0; 1 0 1; 0 0 1];
    out1.fz = [0; 0; 1];
    out1.gs = [1 0; 0 1];
    out1.gx = [0 1 0; 0 0 0];
    out1.hs = [0 0];
    out1.hx = [0 0 0];
    out1.hsnext = [1 1];
    out1.hxnext = [1 0 1];

end