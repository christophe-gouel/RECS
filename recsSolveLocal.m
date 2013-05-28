function X = recsSolveLocal(model)
% RECSSOLVELOCAL Calculates the first-order perturbation solution of a rational expectations model

% Copyright (C) 2011-2013 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
e      = model.w'*model.e;
params = model.params;
sss    = model.sss;
xss    = model.xss;
zss    = model.zss;
[d,m]  = model.dim{1:2};

%% Linearization at steady state
output            = struct('F',0,'Js',1,'Jx',1,'Jz',1,'Jsn',1,'Jxn',1,'hmult',0);
[~,fs,fx,fz]      = model.f(sss,xss,zss,params,output);
[~,gs,gx]         = model.g(sss,xss,e,params,output);
[~,hs,hx,hsn,hxn] = model.h(sss,xss,e,sss,xss,params,output);

% Reshape derivatives
fs  = permute(fs,[2 3 1]);
fx  = permute(fx,[2 3 1]);
fz  = permute(fz,[2 3 1]);
gs  = permute(gs,[2 3 1]);
gx  = permute(gx,[2 3 1]);
hs  = permute(hs,[2 3 1]);
hx  = permute(hx,[2 3 1]);
hsn = permute(hsn,[2 3 1]);
hxn = permute(hxn,[2 3 1]);

%% Definition of matrices A and B, such that A*y_{t+1}=B*y_t with y = [s-sss; x-xss]
A = [-fz*[hsn hxn];
     eye(d) zeros(d,m)];
B = [fs+fz*hs fx+fz*hx;
     gs       gx       ];

%% Complex generalized Schur decomposition
[AA,BB,Q,Z] = qz(A,B);
[AA,BB,Q,Z] = ordqz(AA,BB,Q,Z,'udo');

Z11 = Z(1:d    ,1:d);
Z21 = Z(d+1:end,1:d);

X   = real(Z21/Z11);
