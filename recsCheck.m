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

func   = model.func;
params = model.params;

% Analytical derivatives
[~,fx,fz] = feval(func,'f',s,x,z,[],[],[],params{:});

[~,gx] = feval(func,'g',s,x,[],e,[],[],params{:});

[~,hx,hsnext,hxnext] = feval(func,'h',s,x,[],e,snext,xnext,params{:});

% Numerical derivatives
fxnum = fjac(func,3,'f',s,x,z,[],[],[],params{:});
fznum = fjac(func,4,'f',s,x,z,[],[],[],params{:});

gxnum = fjac(func,3,'g',s,x,[],e,[],[],params{:});

hxnum     = fjac(func,3,'h',s,x,[],e,snext,xnext,params{:});
hsnextnum = fjac(func,6,'h',s,x,[],e,snext,xnext,params{:});
hxnextnum = fjac(func,7,'h',s,x,[],e,snext,xnext,params{:});

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
