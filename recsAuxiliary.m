function [interp,xa,za] = recsAuxiliary(model,interp,s,x,z,xa,options)
% RECSAUXILIARY calculates auxiliary variables not included in the core model
%
% INTERP = RECSAUXILIARY(MODEL,INTERP,S,X,Z) 
% For the definition of MODEL and INTERP structure, see recsSolveREE
% documentation. RECSAUXILIARY returns the interpolation structure INTERP, in
% which the coefficients matrices cxa (for auxiliary variables) and cza (for
% auxiliary expectations) have been added as new fields.
%
% INTERP = RECSAUXILIARY(MODEL,INTERP,S,X,Z,XA) uses XA as first guess for the
% value of the auxiliary variables on the grid. This is only useful when
% auxiliary variables are function of auxiliary expectations. In this case, it
% can reduce a lot the time to reach the solution to start from a good first
% guess.
%
% INTERP = RECSAUXILIARY(MODEL,INTERP,S,X,Z,XA,OPTIONS) solves for auxiliary
% variables with the parameters defined by the structure OPTIONS. The fields of
% the structure are
%    extrapolate      : 1 or 2 if extrapolation is allowed outside the
%                       interpolation space, 0 or -1 to forbid it (default: 1).
%    saoptions        : options structure to be passed to solve SA
%
% [INTERP,XA] = RECSAUXILIARY(MODEL,INTERP,S,X,Z,...) returns the value of the
% auxiliary variables on the grid.
%
% [INTERP,XA,ZA] = RECSAUXILIARY(MODEL,INTERP,S,X,Z,...) returns the value of
% the auxiliary expectations on the grid.
%
% See also RECSSIMUL, RECSSOLVEREE.

% Copyright (C) 2011-2012 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
if nargin <=6, options = struct([]); end

defaultopt = struct(         ...
    'extrapolate', 1        ,...
    'saoptions'  ,struct([]));
warning('off','catstruct:DuplicatesFound')

options = catstruct(defaultopt,options);

e      = model.e;
func   = model.func;
params = model.params;
w      = model.w;

cx     = interp.cx;
fspace = interp.fspace;
Phi    = interp.Phi;

n      = size(s,1);
k      = size(e,1);

%% Phinext
ind           = 1:n;
ind           = ind(ones(1,k),:);
ss            = s(ind,:);
xx            = x(ind,:);
ee            = e(repmat(1:k,1,n),:);
snext         = func('g',ss,xx,[],ee,[],[],params,struct('F',1,'Js',0,'Jx',0));
if options.extrapolate>=1
  snextinterp = snext;
else
  snextinterp = max(min(snext,fspace.b(ones(n*k,1),:)),fspace.a(ones(n*k,1),:));
end
Phinext       = funbasx(fspace,snextinterp);

%% xnext
xnext = funeval(cx,fspace,Phinext);

%% pa (dimension of auxiliary expectations)
za  = func('ha',s(1,:),x(1,:),[],e(1,:),snext(1,:),xnext(1,:),params);
pa  = size(za,2);

%%
if pa>0
  %% With auxiliary expectations function
  if nargin<=5 || isempty(xa)
    xa  = func('xa',s,x,[z zeros(n,pa)],[],[],[],params);
  end
  [xa,~,exitflag] = SA(@ResidualVFI, xa, options.saoptions);

  if exitflag~=1
    warning('RECS:FailureVFI','Failure to converge for auxiliary variables');
  end

  interp.cza = funfitxy(fspace,Phi,za);
else
  %% Without auxiliary expectations function
  xa  = func('xa',s,x,[z zeros(n,pa)],[],[],[],params);
  cxa = funfitxy(fspace,Phi,xa);
  za  = zeros(n,0);
end

%%
interp.cxa = cxa;

%% Nested function
function R = ResidualVFI(xa_old)
% RESIDUALVFI

  cxa    = funfitxy(fspace,Phi,xa_old);
  xanext = funeval(cxa,fspace,Phinext);
  ha     = func('ha',ss,[xx xa_old(ind,:)],[],ee,snext,[xnext xanext],params);
  za     = reshape(w'*reshape(ha,k,n*pa),n,pa);
  xa     = func('xa',s,x,[z za],[],[],[],params);
  R      = xa-xa_old;

end

end