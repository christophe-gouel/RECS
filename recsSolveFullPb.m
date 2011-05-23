function [F,J] = recsSolveFullPb(s,x,func,params,c,e,w,fspace,Phi)

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

[n,m]   = size(x);

[~,grid] = spblkdiag(zeros(m,m,n),[],0);


X = [reshape(x',[n*m 1]); reshape(c',[n*m 1])];
J=[];
%[F,J] = recsFullPb2(X,s,func,params,grid,e,w,fspace,Phi);
[LB,UB] = feval(func,'b',s,[],[],[],[],[],params);
LB = [reshape(LB',[n*m 1]); -inf(n*m,1)];
UB = [reshape(UB',[n*m 1]); +inf(n*m,1)];

% $$$ Jnum = numjac(@(VAR) recsFullPb2(VAR,s,func,params,grid,e,w,fspace,Phi),X);
% $$$ spy(J)
% $$$ figure
% $$$ spy(Jnum)
% $$$ size(Jnum)
% $$$ size(J)
% $$$ 
% $$$ norm(full(J)-Jnum)
% $$$ norm(full(J(1:200,201:400))-Jnum(1:200,201:400))
% $$$ norm(full(J(201:400,201:400))-Jnum(201:400,201:400))
% $$$ return

[X,F,exitflag] = lmmcp(@(VAR) recsFullPb2(VAR,s,func,params,grid,e,w,fspace,Phi),...
                       X,LB,UB);
% $$$ [X,F,exitflag] = lmmcp(@equation,...
% $$$                        X,LB,UB);
if exitflag~=1, disp('No convergence'); end

% $$$ global par
% $$$ par   = {@recsFullPb2,s,func,params,grid,e,w,fspace,Phi};
% $$$ %par   = {@equation};
% $$$ [X,f] = pathmcp(X,LB,UB,'pathtransform');
% $$$ clear global par

x     = reshape(X(1:n*m),m,n)';
c     = reshape(X(n*m+1:end),m,n)';
return
function [F,J] = equation(Variable)
  if nargout==2
    J = numjac(@(VAR) recsFullPb2(VAR,s,func,params,grid,e,w,fspace,Phi), ...
               Variable);
    J = sparse(J);
  end      
  F = recsFullPb2(Variable,s,func,params,grid,e,w,fspace,Phi);
end
end
  
