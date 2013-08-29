function model = recsSolveLocal(model)
% RECSSOLVELOCAL Calculates the first-order perturbation solution of a rational expectations model
%
% RECSSOLVELOCAL implementes a slightly modified version of Klein's (2000) method.
%
% MODEL = RECSSOLVELOCAL(MODEL) calculates the first-order perturbation solution
% of the rational expectations model defined in the object MODEL. RECSSOLVELOCAL
% add to the properties of the object MODEL the structure LinearSolution, which
% contains the linear solution. To be used, RECSSOLVELOCAL requires the
% deterministic steady state to have been found before.
%
% References
% Klein, P. (2000), Using the generalized Schur form to solve a multivariate
% linear rational expectations model, Journal of Economic Dynamics and Control,
% 24(10), 1405-1423.
%
% See also RECSFIRSTGUESS, RECSSS.

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
output            = struct('F',1,'Js',0,'Jx',0);
ge                = numjac(@(E) model.g(sss,xss,E,params,output),e);

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
[~ ,~ ,~,Z] = ordqz(AA,BB,Q,Z,'udo');

%% Calculations of first-order perturbation
% Decision rule
Xs = real(Z(d+1:end,1:d)/Z(1:d,1:d));

% Law of motion
%{
Klein's (2000) formula:
AA11 = AA(1:d,1:d);
BB11 = BB(1:d,1:d);
Ss   = real(Z11*(AA11\BB11)/Z11);
%}
% It is simpler using RECS specific transition rule:
Ss = gs+gx*Xs;
Se = ge;

% Expectations
Zs = hs+hx*Xs+(hsn+hxn*Xs)*Ss;

%% Export to model object
% Matrices
model.LinearSolution   = struct('Xs',Xs,'Ss',Ss,'Se',Se,'Zs',Zs);

% Functions
model.LinearSolution.X = @(s) xss(ones(size(s,1),1),:)+...
                              (Xs*(s-sss(ones(size(s,1),1),:))')';
model.LinearSolution.S = @(s,E) sss(ones(size(s,1),1),:)+...
                                (Ss*(s-sss(ones(size(s,1),1),:))'+...
                                 Se*(E-e(ones(size(s,1),1),:))')';
model.LinearSolution.Z = @(s) zss(ones(size(s,1),1),:)+...
                              (Zs*(s-sss(ones(size(s,1),1),:))')';

