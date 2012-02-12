function [interp,xa,za] = recsAuxilary(model,interp,s,x,z,xa,options)
% RECSAUXILARY calculates auxilary variables not included in the core model
%
% INTERP = RECSAUXILARY(MODEL,INTERP,S,X,Z)
%
% INTERP = RECSAUXILARY(MODEL,INTERP,S,X,Z,XA)
%
% INTERP = RECSAUXILARY(MODEL,INTERP,S,X,Z,XA,OPTIONS)
%
% [INTERP,XA] = RECSAUXILARY(MODEL,INTERP,S,X,Z,...)
%
% [INTERP,XA,ZA] = RECSAUXILARY(MODEL,INTERP,S,X,Z,...)

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

%% pa (dimension of auxilary expectations)
za  = func('ha',s(1,:),x(1,:),[],e(1,:),snext(1,:),xnext(1,:),params);
pa  = size(za,2);

%%
if pa>0
  %% With auxilary expectations function
  if nargin<=5 || isempty(xa)
    xa  = func('xa',s,x,[z zeros(n,pa)],[],[],[],params);
  end
  [xa,~,exitflag] = SA(@ResidualVFI, xa, options.saoptions);

  if exitflag~=1
    warning('RECS:FailureVFI','Failure to converge for auxilary variables');
  end

  interp.cza = funfitxy(fspace,Phi,za);
else
  %% Without auxilary expectations function
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