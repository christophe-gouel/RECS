function [out1,out2,out3,out4,out5] = mf_cs1(flag,s,x,z,e,snext,xnext,p,out);

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
    out1 = zeros(n,1);
    out1(:,1) = -inf;
    out2 = zeros(n,1);
    out2(:,1) = s(:,1);

  case 'f';
    n = size(s,1);

    % f
    if output.F
      out1 = zeros(n,1);
      out1(:,1) = z(:,1).*(p(1) + 1)./(p(2) + 1) - x(:,1).^(-p(3));
    end

    % df/ds
    if output.Js
      out2 = zeros(n,1,1);
      out2(:,1,1) = 0; % d eq_1 w.r.t. X
    end

    % df/dx
    if output.Jx
      out3 = zeros(n,1,1);
      out3(:,1,1) = p(3).*x(:,1).^(-p(3))./x(:,1); % d eq_1 w.r.t. C
    end

    % df/dz
    if output.Jz
      out4 = zeros(n,1,1);
      out4(:,1,1) = (p(1) + 1)./(p(2) + 1); % d eq_1 w.r.t. E
    end
        
  case 'g';
    n = size(s,1);

    % g
    if output.F
      out1 = zeros(n,1);
      out1(:,1) = e(:,1) + (p(1) + 1).*(-x(:,1) + s(:,1));      
    end

    if output.Js
      out2 = zeros(n,1,1);
      out2(:,1,1) = p(1) + 1; % d eq_1 w.r.t. X
    end

    if output.Jx
      out3 = zeros(n,1,1);
      out3(:,1,1) = -p(1) - 1; % d eq_1 w.r.t. C
    end
        
  case 'h';
    n = size(snext,1);

    %h
    if output.F
      out1 = zeros(n,1);
      out1(:,1) = xnext(:,1).^(-p(3));
    end

    if output.Js
      out2 = zeros(n,1,1);
      out2(:,1,1) = 0; % d eq_1 w.r.t. X
    end

    if output.Jx
      out3 = zeros(n,1,1);
      out3(:,1,1) = 0; % d eq_1 w.r.t. C
    end

    if output.Jsn
      out4 = zeros(n,1,1);
      out4(:,1,1) = 0; % d eq_1 w.r.t. X(1)
    end

    if output.Jxn
      out5 = zeros(n,1,1);
      out5(:,1,1) = -p(3).*xnext(:,1).^(-p(3))./xnext(:,1); % d eq_1 w.r.t. C(1)
    end
        
  case 'e';
    warning('Euler equation errors not implemented in Dolo');

  case 'params';
    out1 = [0.05,0.1,2];

  case 'solution';


end
