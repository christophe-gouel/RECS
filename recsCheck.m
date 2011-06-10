function [err,discrepancy] = recsCheck(model,s,x,z,e,snext,xnext)
% RECSCHECK Checks analytical derivatives against numerical ones
%
% RECSCHECK(MODEL,S,X,Z)
%
% RECSCHECK(MODEL,S,X,Z,E)
%
% RECSCHECK(MODEL,S,X,Z,E,SNEXT)
%
% RECSCHECK(MODEL,S,X,Z,E,SNEXT,XNEXT)
%
% ERR = RECSCHECK(MODEL,S,X,Z,...) returns ERR, a 1x6 vector containing the maximal
% differences between analytical and numerical derivatives.
%
% [ERR,DISCREPANCY] = RECSCHECK(MODEL,S,X,Z,...) returns DISCREPANCY, a structure
% containing the detailed differences between analytical and numerical derivatives.

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

if nargin < 7 || isempty(xnext), xnext = x;                end
if nargin < 6 || isempty(snext), snext = s;                end
if nargin < 5 || isempty(e)    , e     = model.w'*model.e; end

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

outputF  = struct('F',1,'Js',0,'Jx',0,'Jz',0,'Jsn',0,'Jxn',0,'hmult',0);
outputJ  = struct('F',0,'Js',0,'Jx',1,'Jz',1,'Jsn',1,'Jxn',1,'hmult',0);

% Analytical derivatives
[~,~,fx,fz] = func('f',s,x,z,[],[],[],params,outputJ);

[~,~,gx] = func('g',s,x,[],e,[],[],params,outputJ);

[~,~,hx,hsnext,hxnext] = func('h',s,x,[],e,snext,xnext,params,outputJ);

% Numerical derivatives
fxnum = numjac(@(X) func('f',s,X,z,[],[],[],params,outputF),x);
fznum = numjac(@(Z) func('f',s,x,Z,[],[],[],params,outputF),z);

gxnum = numjac(@(X) func('g',s,X,[],e,[],[],params,outputF),x);

hxnum = numjac(@(X) func('h',s,x,[],e,snext,xnext,params,outputF),x);
hsnextnum = numjac(@(SNEXT) func('h',s,x,[],e,SNEXT,xnext,params,outputF),snext);
hxnextnum = numjac(@(XNEXT) func('h',s,x,[],e,snext,XNEXT,params,outputF),xnext);

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
