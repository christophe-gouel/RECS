
function [out1,out2,out3,out4,out5,out6] = mf_optimal_growth(flag,s,x,z,e,snext,xnext,p,out);

% p is the vector of parameters

output = struct('F',1,'Js',0,'Jx',0,'Jsn',0,'Jxn',0,'Jz',0, 'hmult',0);

if nargin == 9
    output = catstruct(output,out);
    voidcell                        = cell(1,6);
    [out1,out2,out3,out4,out5,out6] = voidcell{:};
else
    if nargout >= 2
        output.Js = 1;
    end
    if nargout >= 3
        outpout.Jx = 1;
    end
    if nargout >= 4
        if strcompi(flag, 'f')
            output.Jz = 1;
        else
            output.Jsn = 1;
        end
    end
    if nargout >= 5
        output.Jxn = 1;
    end
    if nargout >= 6
        output.hmult = 1;
    end
end
    

switch flag

case 'b';
	n = size(s,1);
	out1 = zeros(n,2);
	out1(:,1) = 0;
	out1(:,2) = -inf;
	out2 = zeros(n,2);
	out2(:,1) = s(:,2);
	out2(:,2) = inf;


case 'f';
	n = size(s,1);
    ep = z;

        % f
        if output.F
		out1 = zeros(n,2);
		out1(:,1) = -1 + (-x(:,1) + s(:,2)).^(-1 + p(4)).*x(:,1).^p(2).*p(1).*p(4).*ep(:,1);
		out1(:,2) = -x(:,2) + p(1).*ep(:,2) + x(:,1).^(1 - p(2))./(1 - p(2));

        end
        % df/ds
        if output.Js
		out2 = zeros(n,2,2);
		out2(:,1,1) = (-x(:,1) + s(:,2)).^(-1 + p(4)).*x(:,1).^p(2).*(1 - p(4)).*p(1).*p(4).*ep(:,1)./(-x(:,1) + s(:,2)) + (-x(:,1) + s(:,2)).^(-1 + p(4)).*x(:,1).^p(2).*p(1).*p(2).*p(4).*ep(:,1)./x(:,1); % d eq_1 w.r.t. c
		out2(:,1,2) = 0; % d eq_1 w.r.t. Welf
		out2(:,2,1) = x(:,1).^(1 - p(2))./x(:,1); % d eq_2 w.r.t. c
		out2(:,2,2) = -1; % d eq_2 w.r.t. Welf

        end
        % df/dx
        if output.Jx
		out3 = zeros(n,2,2);
		out3(:,1,1) = (-x(:,1) + s(:,2)).^(-1 + p(4)).*x(:,1).^p(2).*p(1).*p(4); % d eq_1 w.r.t. eh
		out3(:,1,2) = 0; % d eq_1 w.r.t. eW
		out3(:,2,1) = 0; % d eq_2 w.r.t. eh
		out3(:,2,2) = p(1); % d eq_2 w.r.t. eW

        end
        % df/dz
        if output.Jz
		out4 = zeros(n,2,2);
		out4(:,1,1) = 0; % d eq_1 w.r.t. A
		out4(:,1,2) = -(-x(:,1) + s(:,2)).^(-1 + p(4)).*x(:,1).^p(2).*(1 - p(4)).*p(1).*p(4).*ep(:,1)./(-x(:,1) + s(:,2)); % d eq_1 w.r.t. w
		out4(:,2,1) = 0; % d eq_2 w.r.t. A
		out4(:,2,2) = 0; % d eq_2 w.r.t. w

        end
        

case 'g';
	n = size(s,1);

        % g
        if output.F
		out1 = zeros(n,2);
		out1(:,1) = 1 - (1 - s(:,1)).*p(3) + e(:,1);
		out1(:,2) = 1 + (-x(:,1) + s(:,2)).^p(4).*(1 - (1 - s(:,1)).*p(3) + e(:,1));

        end
        if output.Js
		out2 = zeros(n,2,2);
		out2(:,1,1) = 0; % d eq_1 w.r.t. c
		out2(:,1,2) = 0; % d eq_1 w.r.t. Welf
		out2(:,2,1) = -(-x(:,1) + s(:,2)).^p(4).*(1 - (1 - s(:,1)).*p(3) + e(:,1)).*p(4)./(-x(:,1) + s(:,2)); % d eq_2 w.r.t. c
		out2(:,2,2) = 0; % d eq_2 w.r.t. Welf

        end
        if output.Jx
		out3 = zeros(n,2,2);
		out3(:,1,1) = p(3); % d eq_1 w.r.t. A
		out3(:,1,2) = 0; % d eq_1 w.r.t. w
		out3(:,2,1) = (-x(:,1) + s(:,2)).^p(4).*p(3); % d eq_2 w.r.t. A
		out3(:,2,2) = (-x(:,1) + s(:,2)).^p(4).*(1 - (1 - s(:,1)).*p(3) + e(:,1)).*p(4)./(-x(:,1) + s(:,2)); % d eq_2 w.r.t. w

        end
        

case 'h';
	n = size(snext,1);

        %h
        if output.F
		out1 = zeros(n,2);
		out1(:,1) = xnext(:,1).^(-p(2)).*snext(:,1);
		out1(:,2) = xnext(:,2);

        end
        if output.Js
		out2 = zeros(n,2,2);
		out2(:,1,1) = 0; % d eq_1 w.r.t. A
		out2(:,1,2) = 0; % d eq_1 w.r.t. w
		out2(:,2,1) = 0; % d eq_2 w.r.t. A
		out2(:,2,2) = 0; % d eq_2 w.r.t. w

        end
        if output.Jx
		out3 = zeros(n,2,2);
		out3(:,1,1) = 0; % d eq_1 w.r.t. c
		out3(:,1,2) = 0; % d eq_1 w.r.t. Welf
		out3(:,2,1) = 0; % d eq_2 w.r.t. c
		out3(:,2,2) = 0; % d eq_2 w.r.t. Welf

        end
        if output.Jsn
		out4 = zeros(n,2,2);
		out4(:,1,1) = xnext(:,1).^(-p(2)); % d eq_1 w.r.t. A(1)
		out4(:,1,2) = 0; % d eq_1 w.r.t. w(1)
		out4(:,2,1) = 0; % d eq_2 w.r.t. A(1)
		out4(:,2,2) = 0; % d eq_2 w.r.t. w(1)

        end
        if output.Jxn
		out5 = zeros(n,2,2);
		out5(:,1,1) = -xnext(:,1).^(-p(2)).*p(2).*snext(:,1)./xnext(:,1); % d eq_1 w.r.t. c(1)
		out5(:,1,2) = 0; % d eq_1 w.r.t. Welf(1)
		out5(:,2,1) = 0; % d eq_2 w.r.t. c(1)
		out5(:,2,2) = 1; % d eq_2 w.r.t. Welf(1)

        end
        if output.hmult
            out6 = ones(size(out1));
            %warning('hmult not implemented')
        end
        

case 'e';
	warning('Euler equation errors not implemented in Dolo');

case 'params';
	out1 = [0.96,4.0,0.9,0.3];

case 'solution';
	sol = struct;
        sol.s_ss = [ 1.          1.58655814];
        sol.x_ss = [ 1.4176294  -2.92503201];
        sol.X = [[  9.00000000e-01  -1.23916284e-16 ;   5.08244059e-01   3.27261015e-01]];
        sol.Y = [[-0.12452152 ;  1.00375066]];
        sol.Z = [[ 0.01887194  0.68582943 ;  0.92264591  0.24759916]];
        out1 = sol;



end;
