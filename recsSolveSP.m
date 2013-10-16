function [x,exitflag] = recsSolveSP(s,x,b,bnext,f,g,h,params,dim,cxnext,e,w, ...
                                    fspacenext,options)

eqsolver         = lower(options.eqsolver);
eqsolveroptions  = options.eqsolveroptions;
extrapolate      = options.extrapolate;
funapprox        = lower(options.funapprox);
loop_over_s      = options.loop_over_s;

[~,m,p] = dim{:};

n                = size(x,1);
[LB,UB]          = b(s,params);

x              = reshape(x',[n*m 1]);
LB             = reshape(LB',[n*m 1]);
UB             = reshape(UB',[n*m 1]);
 
[x,~,exitflag] = runeqsolver(@equilibriumsp,x,LB,UB,eqsolver,eqsolveroptions,...
                             s,bnext,f,g,h,e,w,m,p,params,cxnext,fspacenext);

if exitflag~=1, disp('No convergence'); end

x              = reshape(x,m,n)';


function F = equilibriumsp(x,s,bnext,f,g,h,e,w,m,p,params,cxnext,fspacenext)

n     = size(s,1);
x     = reshape(x,m,n)';

k     = size(e,1);
ind   = (1:n);
ind   = ind(ones(1,k),:);
ss    = s(ind,:);
xx    = x(ind,:);
ee    = e(repmat(1:k,1,n),:);

snext   = g(ss,xx,ee,params);

[LBnext,UBnext] = bnext(snext,params);
xnext = min(max(funeval(cxnext,fspacenext,snext),LBnext),UBnext);
hv                 = h(ss,xx,ee,snext,xnext,params);

z       = reshape(w'*reshape(hv,k,n*p),n,p);
F       = f(s,x,z,params);
F = reshape(F',n*m,1);

% disp('Min')
% 
% [min(snext);
%  fspacenext.a]
% 
% disp('Max')
% [max(snext);
%  fspacenext.b]

return