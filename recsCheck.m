function [err,discrepancy] = recsCheck(model,s,x,z,e,snext,xnext)
% RECSCHECK Checks analytical derivatives against numerical ones
% OUTPUT
%   err : a 1x6 vector containing the maximal differences between analytical and
%   numerical derivatives
%   discrepancy : a structure containing the detailed differences between
%   analytical and numerical derivatives

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt


if nargin < 6
  snext = s;
  xnext = x;
end

if nargin < 5 || isempty(e)
  e = model.w'*model.e;
end

if size(s,1)~=1 || size(x,1)~=1 || size(e,1)~=1 || size(snext,1)~=1 || size(xnext,1)~=1
  error('Derivatives can only be check at one location, not on a grid');
end

params = model.params;
if isa(model.func,'char')
  func = str2func(model.func);
elseif isa(model.func,'function_handle')
  func = model.func;
else
  error('model.func must be either a string or a function handle')
end

% Analytical derivatives
[~,fx,fz] = func('f',s,x,z,[],[],[],params);

[~,gx] = func('g',s,x,[],e,[],[],params);

[~,hx,hsnext,hxnext] = func('h',s,x,[],e,snext,xnext,params);

% Numerical derivatives
fxnum = numjac(@(X) func('f',s,X,z,[],[],[],params),x);
fznum = numjac(@(Z) func('f',s,x,Z,[],[],[],params),z);

gxnum = numjac(@(X) func('g',s,X,[],e,[],[],params),x);

hxnum = numjac(@(X) func('h',s,x,[],e,snext,xnext,params),x);
hsnextnum = numjac(@(SNEXT) func('h',s,x,[],e,SNEXT,xnext,params),snext);
hxnextnum = numjac(@(XNEXT) func('h',s,x,[],e,snext,XNEXT,params),xnext);

% Error
err = norm(fx(:)-fxnum(:),inf);
err = [err norm(fz(:)-fznum(:),inf)];

err = [err norm(gx(:)-gxnum(:),inf)];

err = [err norm(hx(:)-hxnum(:),inf)];
err = [err norm(hsnext(:)-hsnextnum(:),inf)];
err = [err norm(hxnext(:)-hxnextnum(:),inf)];

if max(err)>1.e-4
   disp('Possible Error in Derivatives')
   disp('Discrepancies in derivatives = ')
   disp('    fx        fz        gx             hx      hsnext hxnext');
   disp(err)
end

discrepancy    = struct(...
    'fx', reshape(fx,size(fxnum))-fxnum, ...
    'fz', reshape(fz,size(fznum))-fznum, ...
    'gx', reshape(gx,size(gxnum))-gxnum, ...
    'hx', reshape(hx,size(hxnum))-hxnum, ...
    'hsnext', reshape(hsnext,size(hsnextnum))-hsnextnum, ...
    'hxnext', reshape(hxnext,size(hxnextnum))-hxnextnum);
