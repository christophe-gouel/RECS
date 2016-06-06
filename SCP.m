function [x,F,exitflag,N] = SCP(x,z1,z0,problem,N)
% SCP Solves a problem through simple continuation method
%
% SCP is to be used over another solver. It makes successive call to the
% underlying solver while adjusting its parameters. The idea is to start from
% parameters, Z0, with which the problem is easy to solve in order to converge
% to the final parameters, Z1, with which the problem may be more difficult to
% solve.
%
% X = SCP(X,Z1,Z0,PROBLEM) attempts to solve for X the problem PROBLEM(X,Z1),
% where Z1 designates parameters of the problem and X a first guess of the
% solution. If it fails, it will try in N=2 steps by solving iteratively
% PROBLEM(X,Z) with Z=((N-n)*Z0+n*Z1)/N, with Z0 a set of parameters for which
% the problem is easy to solve or for which the first guess X is the
% solution. In case, it fails, N is increased until it converges or until 10
% attempts have been tried.
%
% X = SCP(X,Z1,Z0,PROBLEM,N) starts by slicing the problem in N slices.
%
% [X,F] = SCP(X,Z1,Z0,PROBLEM,...) returns F, the value returned by PROBLEM at
% the solution point.
%
% [X,F,EXITFLAG] = SCP(X,Z1,Z0,PROBLEM,...) returns EXITFLAG that describes the exit
% conditions. Possible values are
%      1         : SCP converged
%      0         : Failure to converge
%
% [X,F,EXITFLAG,N] = SCP(X,Z1,Z0,PROBLEM,...) returns N the number of slices in
% which the problem was sliced to achieve convergence.

% Copyright (C) 2011-2016 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt
  
if nargin<=4, N = 1; end
exitflag  = 0;
it        = 0;
maxit     = 10;
Nlist     = 1:maxit;
lastgoodx = x;
lastgoodz = z0;
while(exitflag~=1 && it < maxit)
  it = it + 1;
  x  = lastgoodx;
  z0 = lastgoodz;
  for n=1:N
    z = ((N - n)*z0 + n*z1)/N;
    [x,F,exitflag] = problem(x,z);
    if exitflag==1
      lastgoodx = x;
      lastgoodz = z;
    else
      break
    end
  end
  N = N + Nlist(it);
end

N = N - Nlist(it);
if exitflag~=1 && it==maxit
  exitflag = 0;
else
  exitflag = 1;
end
