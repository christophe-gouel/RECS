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
%    eqsolver         : 'fsolve', 'krylov' (default), 'lmmcp', 'ncpsolve',
%                       'path' or 'sa'
%    extrapolate      : 1 or 2 if extrapolation is allowed outside the
%                       interpolation space, 0 or -1 to forbid it (default: 1).
%    eqsolveroptions  : options structure to be passed to eqsolver
%
% [INTERP,XA] = RECSAUXILIARY(MODEL,INTERP,S,X,Z,...) returns the value of the
% auxiliary variables on the grid.
%
% [INTERP,XA,ZA] = RECSAUXILIARY(MODEL,INTERP,S,X,Z,...) returns the value of
% the auxiliary expectations on the grid.
%
% See also RECSSIMUL, RECSSOLVEREE.

% Copyright (C) 2011-2013 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
defaultopt = struct(                                      ...
    'eqsolver'       , 'krylov'                          ,...
    'eqsolveroptions', struct('Diagnostics'    , 'off' ,...
                              'DerivativeCheck', 'off' ,...
                              'Jacobian'       , 'off')  ,...
    'extrapolate'    , 1);
if nargin <=6
  options = defaultopt;
else
  warning('off','catstruct:DuplicatesFound')
  if isfield(options,'eqsolveroptions')
    options.eqsolveroptions = catstruct(defaultopt.eqsolveroptions,...
                                        options.eqsolveroptions);
  end
  options = catstruct(defaultopt,options);
end
eqsolver        = lower(options.eqsolver);
eqsolveroptions = options.eqsolveroptions;

b      = model.functions.b;
e      = model.shocks.e;
g      = model.functions.g;
params = model.params;
w      = model.shocks.w;

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
snext         = g(ss,xx,ee,params,struct('F',1,'Js',0,'Jx',0));
if options.extrapolate>=1
  snextinterp = snext;
else
  snextinterp = max(min(snext,fspace.b(ones(n*k,1),:)),...
                    fspace.a(ones(n*k,1),:));
end
Phinext       = funbas(fspace,snextinterp);
% Phinext       = funbasx(fspace,snextinterp);

%% xnext
[LBnext,UBnext] = b(snext,params);
xnext           = min(max(Phinext*cx,LBnext),UBnext);
% xnext = min(max(funeval(cx,fspace,Phinext),LBnext),UBnext);

%%
[ma,pa] = model.dima{:};
fa      = model.functions.fa;
za      = zeros(n,pa);

%%
if pa>0
  %% With auxiliary expectations function
  if nargin<=5 || isempty(xa), xa  = fa(s,x,z,za,params); end
  [xa,~,exitflag] = runeqsolver(@(X) ResidualVFI(X,false),xa(:),...
                                -inf(numel(xa),1),inf(numel(xa),1),...
                                eqsolver,eqsolveroptions);
  xa              = reshape(xa,n,ma);

  if exitflag~=1
    warning('RECS:FailureVFI','Failure to converge for auxiliary variables');
  end

  interp.cza = funfitxy(fspace,Phi,za);
else
  %% Without auxiliary expectations function
  xa  = fa(s,x,z,za,params);
  cxa = funfitxy(fspace,Phi,xa);
end

%%
interp.cxa = cxa;

%% Nested function
function R = ResidualVFI(xa_old,serial)
% RESIDUALVFI

  if ~serial, xa_old = reshape(xa_old,n,ma); end

  cxa    = funfitxy(fspace,Phi,xa_old);
  xanext = Phinext*cxa;
%   xanext = funeval(cxa,fspace,Phinext);
  ha     = model.functions.ha(ss,xx,xa_old(ind,:),ee,snext,xnext,xanext,params);
  za     = reshape(w'*reshape(ha,k,n*pa),n,pa);
  xa     = fa(s,x,z,za,params);
  R      = xa-xa_old;
  if ~serial, R = R(:); end

end

end