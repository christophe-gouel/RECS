function [se,lEE,lEf] = recsAccuracy(model,interp,s,options)
% RECSACCURACY Evaluates approximation accuracy
%
% SE = RECSACCURACY(MODEL,INTERP,S) evaluates the accuracy of the approximation
% defined in the interpolition structure INTERP for the model defined in the
% object MODEL. The accuracy is assessed over the state variables contained
% in the n-by-d-by-y array S, as output by resSimul. RECSACCURACY returns a
% n*y-by-d matrix SE containing the state variables on which the accuracy was
% evaluated.
% MODEL is an object created by recsmodel.
% INTERP is a structure, which includes the following fields:
%    cx      : coefficient matrix of the interpolation of the response variables
%    fspace  : a definition structure for the interpolation family (created by
%              the function fundef)
%
% SE = RECSACCURACY(MODEL,INTERP,S,OPTIONS) evaluates the accuracy with the
% parameters defined by the structure OPTIONS. The fields of the
% structure are
%    extrapolate      : 1 if extrapolation is allowed outside the
%                       interpolation space or 0 to forbid it (default: 1)
%
% [SE,LEE] = RECSACCURACY(MODEL,INTERP,S,...) returns the matrix LEE containing
% the value of Euler equation error (in log10) evaluated on the grid points SE.
%
% [SE,LEE,lEF] = RECSACCURACY(MODEL,INTERP,S,...) returns the nxy-by-m matrix
% lEF containing the value of equilibrium equation error (in log10) evaluated on
% the grid points SE. For models with finite bounds, this error is assessed
% using a minmax formulation: Ef = abs(min(max(-f,LB-x),UB-x)).
%
% See also RECSSIMUL, RECSSOLVEREE.

% Copyright (C) 2011-2013 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
% Options
defaultopt  = struct('display'    , 1,...
                     'extrapolate',1);
if nargin <4
  options   = defaultopt;
else
  warning('off','catstruct:DuplicatesFound')
  options   = catstruct(defaultopt,options);
end
display     = options.display;
extrapolate = options.extrapolate;

b         = model.functions.b;
e         = model.e;
f         = model.functions.f;
g         = model.functions.g;
h         = model.functions.h;
ixforward = model.ixforward;
params    = model.params;
w         = model.w;

cx     = interp.cx;
fspace = interp.fspace;

[d,m,p] = model.dim{:};
validateattributes(s,{'numeric'},{'size',[NaN,d,NaN],'nonempty'},3)
[n,~,t] = size(s);
k       = size(e,1);

se      = permute(s,[1 3 2]);
se      = reshape(se,n*t,d);

[LB,UB]    = b(se,params);
if extrapolate, seinterp = se;
else
  seinterp = max(min(se,fspace.b(ones(n*t,1),:)),fspace.a(ones(n*t,1),:));
end
xe         = min(max(funeval(cx,fspace,seinterp),LB),UB);

%% Calculation of ze
ind       = (1:n*t);
ind       = ind(ones(1,k),:);
ss        = se(ind,:);
xx        = xe(ind,:);
ee        = e(repmat(1:k,1,n*t),:);
sen       = g(ss,xx,ee,params);
[LBn,UBn] = b(sen,params);
if extrapolate, seninterp = sen;
else
  seninterp = max(min(sen,fspace.b(ones(n*t*k,1),:)),fspace.a(ones(n*t*k,1),:));
end
xen              = zeros(n*t*k,m);
xen(:,ixforward) = min(max(funeval(cx(:,ixforward),fspace,seninterp),...
                           LBn(:,ixforward)),UBn(:,ixforward));
hv               = h(ss,xx,ee,sen,xen,params);
ze               = reshape(w'*reshape(hv,k,n*t*p),n*t,p);

if display==1, disp('Accuracy of the solution'); end

%% Euler equation error
EE      = model.functions.ee(se,xe,ze,params);
lEE     = log10(abs(EE));
lEE_res = [log10(max(abs(EE)));
           log10(sum(abs(EE))/size(EE,1))];

if display==1 && all(~isnan(lEE_res))
  disp(' Euler equation error (in log10)');
  disp('    Max       Mean');
  disp(lEE_res');
end

%% Equilibrium equation error
fe      = f(se,xe,ze,params);

Ef      = abs(min(max(-fe,LB-xe),UB-xe));
lEf     = log10(Ef);
Ef_res  = [log10(max(Ef));
           log10(sum(Ef)/size(Ef,1))];

if display==1
  disp(' Equilibrium equation error (in log10 units)');
  disp('    Max       Mean');
  disp(Ef_res');
end

