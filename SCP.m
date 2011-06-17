function [x,F,exitflag] = SCP(x,z1,z0,problem,N)
% SCP Solves a problem through simple continuation method
%
% SCP does not increase the number of steps in case of failure when the problem
% is solved with ncpsolve, since ncpsolve does not output an exitflag.
%
% X = SCP(X,Z1,Z0,PROBLEM,N)  
%
% [X,F] = SCP(X,Z1,Z0,PROBLEM,N)  
%
% [X,F,EXITFLAG] = SCP(X,Z1,Z0,PROBLEM,N) returns EXITFLAG that describes the exit
% conditions. Possible values are
%      1         : SCP converged
%      0         : Failure to converge

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt
  
exitflag = 0;
it       = 0;
x0       = x;
maxit    = 10;
while(exitflag~=1 && it < maxit)
  it = it+1;
  x  = x0;
  for n=1:N
    z = ((N-n)*z0 + n*z1)/N;
    [x,F,exitflag] = problem(x,z);
    if exitflag~=1, break, end
  end
  N = N+1;
end

if exitflag~=1 && it==maxit
  exitflag = 0;
else
  exitflag = 1;
end
