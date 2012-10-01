function [out1,out2,out3,out4,out5] = gro1model(flag,s,x,z,e,snext,xnext,p,output)

voidcell                   = cell(1,5);
[out1,out2,out3,out4,out5] = voidcell{:};

switch flag

  case 'b'
    n = size(s,1);
    out1 = -inf(n,1);
    out2 = inf(n,1);

  case 'f'
    n = size(s,1);

    % f
    if output.F
      out1 = zeros(n,1);
      out1(:,1) = -p(4).*z(:,1) + x(:,1).^(-p(2));
    end

    % df/ds
    if output.Js
      out2 = zeros(n,1,2);
    end

    % df/dx
    if output.Jx
      out3 = zeros(n,1,1);
      out3(:,1,1) = -p(2).*x(:,1).^(-p(2))./x(:,1); % d eq_1 w.r.t. C
    end

    % df/dz
    if output.Jz
      out4 = zeros(n,1,1);
      out4(:,1,1) = -p(4); % d eq_1 w.r.t. E
    end
  
  case 'g'
    n = size(s,1);

    % g
    if output.F
      out1 = zeros(n,2);
      out1(:,1) = p(1).*s(:,1).^p(6).*exp(s(:,2)) - x(:,1) + s(:,1).*(-p(3) + 1);
      out1(:,2) = p(5).*s(:,2) + e(:,1);      
    end

    if output.Js
      out2 = zeros(n,2,2);
      out2(:,1,1) = p(1).*p(6).*s(:,1).^p(6).*exp(s(:,2))./s(:,1) - p(3) + 1; % d eq_1 w.r.t. K
      out2(:,1,2) = p(1).*s(:,1).^p(6).*exp(s(:,2)); % d eq_1 w.r.t. Z
      out2(:,2,2) = p(5); % d eq_2 w.r.t. Z
    end

    if output.Jx
      out3 = zeros(n,2,1);
      out3(:,1,1) = -1; % d eq_1 w.r.t. C
    end
  
  case 'h'
    n = size(snext,1);

    %h
    if output.F
      out1 = zeros(n,1);
      out1(:,1) = xnext(:,1).^(-p(2)).*(p(1).*p(6).*snext(:,1).^(p(6) - 1).*exp(snext(:,2)) - p(3) + 1);
    end

    if output.Js
      out2 = zeros(n,1,2);
    end

    if output.Jx
      out3 = zeros(n,1,1);
    end

    if output.Jsn
      out4 = zeros(n,1,2);
      out4(:,1,1) = p(1).*p(6).*xnext(:,1).^(-p(2)).*snext(:,1).^(p(6) - 1).*(p(6) - 1).*exp(snext(:,2))./snext(:,1); % d eq_1 w.r.t. K(1)
      out4(:,1,2) = p(1).*p(6).*xnext(:,1).^(-p(2)).*snext(:,1).^(p(6) - 1).*exp(snext(:,2)); % d eq_1 w.r.t. Z(1)
    end

    if output.Jxn
      out5 = zeros(n,1,1);
      out5(:,1,1) = -p(2).*xnext(:,1).^(-p(2)).*(p(1).*p(6).*snext(:,1).^(p(6) - 1).*exp(snext(:,2)) - p(3) + 1)./xnext(:,1); % d eq_1 w.r.t. C(1)
    end
  
  case 'e'
    out1 = [];

  case 'params'
    out1 = [ 0.21888357  2.          0.0196      0.95        0.9         0.33      ];

  case 'ss'
    out1 = [1.0  0.0];
    out2 = [0.199283572568];

  case 'J'
    out1.fs = [0 0];
    out1.fx = [1];
    out1.fz = [1];
    out1.gs = [1 1; 0 1];
    out1.gx = [1; 0];
    out1.hs = [0 0];
    out1.hx = [0];
    out1.hsnext = [1 1];
    out1.hxnext = [1];

end