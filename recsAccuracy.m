function recsAccuracy(model,interp,s,options)
% RECSACCURACY Evaluates approximation accuracy
%
% RECSACCURACY(MODEL,INTERP,S) evaluates the accuracy of the approximation
% defined in the interpolition structure INTERP for the model defined in the
% structure MODEL. The accuracy is assessed over the state variables contained
% in the n-by-d-by-y array S.
% MODEL is a structure, which includes the following fields:
%    func    : function name or anonymous function that defines the model's equations
%    params  : model's parameters, it is preferable to pass them as a cell array
%              (compulsory with the functional option) but other formats are
%              acceptable
% INTERP is a structure, which includes the following fields:
%    cx      : coefficient matrix of the interpolation of the response variables
%    fspace  : a definition structure for the interpolation family (created by
%              the function fundef)
%
% RECSACCURACY(MODEL,INTERP,S,OPTIONS) evaluates the accuracy with the
% parameters defined by the structure OPTIONS. The fields of the
% structure are  
%    extrapolate      : 1 if extrapolation is allowed outside the
%                       interpolation space or 0 to forbid it (default: 1)
%
% See also RECSSIMUL, RECSSOLVEREE.

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
% Options
if nargin <=4, options = struct([]); end
defaultopt  = struct('extrapolate',1);
warning('off','catstruct:DuplicatesFound')
options     = catstruct(defaultopt,options);
extrapolate = options.extrapolate;

e      = model.e;
if isa(model.func,'char')
  func = str2func(model.func);
elseif isa(model.func,'function_handle')
  func = model.func;
else
  error('model.func must be either a string or a function handle')
end
params = model.params;
w      = model.w;

cx     = interp.cx;
fspace = interp.fspace;

[n,d,t] = size(s);
k       = size(e,1);

se      = permute(s,[1 3 2]);
se      = reshape(se,n*t,d);

[LB,UB]    = func('b',se,[],[],[],[],[],params);
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
output    = struct('F',1,'Js',0,'Jx',0,'Jz',0,'Jsn',0,'Jxn',0,'hmult',1);
sen       = func('g',ss,xx,[],e(repmat(1:k,1,n*t),:),[],[],params,output);
[LBn,UBn] = func('b',sen,[],[],[],[],[],params);
if extrapolate, seninterp = sen;
else      
  seninterp = max(min(sen,fspace.b(ones(n*t*k,1),:)),fspace.a(ones(n*t*k,1),:)); 
end
xen       = min(max(funeval(cx,fspace,seninterp),LBn),UBn);
if nargout(func)<6
  h                 = func('h',ss,xx,[],e(repmat(1:k,1,n*t),:),sen,xen,params,output);
else
  [h,~,~,~,~,hmult] = func('h',ss,xx,[],e(repmat(1:k,1,n*t),:),sen,xen,params,output);
  h                 = h.*hmult;
end
p         = size(h,2);
ze        = reshape(w'*reshape(h,k,n*t*p),n*t,p);

disp('Accuracy of the solution');

%% Euler equation error
EE      = func('e',se,xe,ze,[],[],[],params);
lEE     = [log10(max(abs(EE)));
           log10(sum(abs(EE))/size(EE,1))];

if ~isnan(lEE)
  disp(' Euler equation error (in log10)');
  disp('    Max       Mean');
  disp(lEE');
end

%% Equilibrium equation error
fe      = func('f',se,xe,ze,[],[],[],params,output);

Ef = min(max(-fe,LB-xe),UB-xe);
Ef = [max(abs(Ef));
      mean(abs(Ef))];

disp(' Equilibrium equation error (minmax formulation)');
disp('    Max       Mean');
disp(Ef');


