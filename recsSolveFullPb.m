function [F,J] = recsSolveFullPb(s,x,func,params,c,e,w,fspace,Phi)

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

[n,m]   = size(x);

[~,grid] = spblkdiag(zeros(m,m,n),[],0);

[F,J] = recsFullPb2(x,s,func,params,grid,c,e,w,fspace,Phi);

