function [out1,out2,out3,out4,out5] = mf_gro2(flag,s,x,z,e,snext,xnext,p,out);

output = struct('F',1,'Js',0,'Jx',0,'Jsn',0,'Jxn',0,'Jz',0,'hmult',0);

if nargin == 9
  output                     = catstruct(output,out);
  voidcell                   = cell(1,5);
  [out1,out2,out3,out4,out5] = voidcell{:};
else
  if nargout >= 2, output.Js = 1; end
  if nargout >= 3, output.Jx = 1; end
  if nargout >= 4
    if strcmpi(flag, 'f')
      output.Jz = 1;
    else
      output.Jsn = 1;
    end
  end
  if nargout >= 5, output.Jxn = 1; end
end


switch flag

  case 'b';
    n = size(s,1);
    out1 = zeros(n,3);
    out1(:,1) = -inf;
    out1(:,2) = 0;
    out1(:,3) = -inf;
    out2 = zeros(n,3);
    out2(:,1) = inf;
    out2(:,2) = inf;
    out2(:,3) = inf;

  case 'f';
    n = size(s,1);

    % f
    if output.F
      out1 = zeros(n,3);
      out1(:,1) = x(:,1) + x(:,2) - s(:,1).^p(6).*p(1).*exp(s(:,2));
      out1(:,2) = x(:,3) + x(:,1).^(-p(2));
      out1(:,3) = p(4).*z(:,1) + x(:,3);
    end

    % df/ds
    if output.Js
      out2 = zeros(n,3,2);
      out2(:,1,1) = -s(:,1).^p(6).*p(1).*p(6).*exp(s(:,2))./s(:,1); % d eq_1 w.r.t. K
      out2(:,1,2) = -s(:,1).^p(6).*p(1).*exp(s(:,2)); % d eq_1 w.r.t. Z
      out2(:,2,1) = 0; % d eq_2 w.r.t. K
      out2(:,2,2) = 0; % d eq_2 w.r.t. Z
      out2(:,3,1) = 0; % d eq_3 w.r.t. K
      out2(:,3,2) = 0; % d eq_3 w.r.t. Z
    end

    % df/dx
    if output.Jx
      out3 = zeros(n,3,3);
      out3(:,1,1) = 1; % d eq_1 w.r.t. C
      out3(:,1,2) = 1; % d eq_1 w.r.t. I
      out3(:,1,3) = 0; % d eq_1 w.r.t. Mu
      out3(:,2,1) = -x(:,1).^(-p(2)).*p(2)./x(:,1); % d eq_2 w.r.t. C
      out3(:,2,2) = 0; % d eq_2 w.r.t. I
      out3(:,2,3) = 1; % d eq_2 w.r.t. Mu
      out3(:,3,1) = 0; % d eq_3 w.r.t. C
      out3(:,3,2) = 0; % d eq_3 w.r.t. I
      out3(:,3,3) = 1; % d eq_3 w.r.t. Mu
    end

    % df/dz
    if output.Jz
      out4 = zeros(n,3,1);
      out4(:,1,1) = 0; % d eq_1 w.r.t. E
      out4(:,2,1) = 0; % d eq_2 w.r.t. E
      out4(:,3,1) = p(4); % d eq_3 w.r.t. E
    end
        
  case 'g';
    n = size(s,1);

    % g
    if output.F
      out1 = zeros(n,2);
      out1(:,1) = (1 - p(3)).*s(:,1) + x(:,2);
      out1(:,2) = p(5).*s(:,2) + e(:,1);      
    end

    if output.Js
      out2 = zeros(n,2,2);
      out2(:,1,1) = 1 - p(3); % d eq_1 w.r.t. K
      out2(:,1,2) = 0; % d eq_1 w.r.t. Z
      out2(:,2,1) = 0; % d eq_2 w.r.t. K
      out2(:,2,2) = p(5); % d eq_2 w.r.t. Z
    end

    if output.Jx
      out3 = zeros(n,2,3);
      out3(:,1,1) = 0; % d eq_1 w.r.t. C
      out3(:,1,2) = 1; % d eq_1 w.r.t. I
      out3(:,1,3) = 0; % d eq_1 w.r.t. Mu
      out3(:,2,1) = 0; % d eq_2 w.r.t. C
      out3(:,2,2) = 0; % d eq_2 w.r.t. I
      out3(:,2,3) = 0; % d eq_2 w.r.t. Mu
    end
        
  case 'h';
    n = size(snext,1);

    %h
    if output.F
      out1 = zeros(n,1);
      out1(:,1) = -(1 - p(3)).*xnext(:,3) + xnext(:,1).^(-p(2)).*snext(:,1).^(-1 + p(6)).*p(1).*p(6).*exp(snext(:,2));
    end

    if output.Js
      out2 = zeros(n,1,2);
      out2(:,1,1) = 0; % d eq_1 w.r.t. K
      out2(:,1,2) = 0; % d eq_1 w.r.t. Z
    end

    if output.Jx
      out3 = zeros(n,1,3);
      out3(:,1,1) = 0; % d eq_1 w.r.t. C
      out3(:,1,2) = 0; % d eq_1 w.r.t. I
      out3(:,1,3) = 0; % d eq_1 w.r.t. Mu
    end

    if output.Jsn
      out4 = zeros(n,1,2);
      out4(:,1,1) = -xnext(:,1).^(-p(2)).*snext(:,1).^(-1 + p(6)).*(1 - p(6)).*p(1).*p(6).*exp(snext(:,2))./snext(:,1); % d eq_1 w.r.t. K(1)
      out4(:,1,2) = xnext(:,1).^(-p(2)).*snext(:,1).^(-1 + p(6)).*p(1).*p(6).*exp(snext(:,2)); % d eq_1 w.r.t. Z(1)
    end

    if output.Jxn
      out5 = zeros(n,1,3);
      out5(:,1,1) = -xnext(:,1).^(-p(2)).*snext(:,1).^(-1 + p(6)).*p(1).*p(6).*p(2).*exp(snext(:,2))./xnext(:,1); % d eq_1 w.r.t. C(1)
      out5(:,1,2) = 0; % d eq_1 w.r.t. I(1)
      out5(:,1,3) = -1 + p(3); % d eq_1 w.r.t. Mu(1)
    end
        
  case 'e';
    warning('Euler equation errors not implemented in Dolo');

  case 'params';
    out1 = [0.240771929825,1,0.0196,0.95,0.9,0.3];

  case 'solution';


end
