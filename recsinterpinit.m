function [interp,s] = recsinterpinit(n,smin,smax,method,options)
% RECSINTERPINIT Prepares an interpolation structure for RECS
%
% [INTERP,S] = RECSINTERPINIT(N,SMIN,SMAX)
%
% [INTERP,S] = RECSINTERPINIT(N,SMIN,SMAX,METHOD)
%
% [INTERP,S] = RECSINTERPINIT(N,SMIN,SMAX,METHOD,OPTIONS)

% Copyright (C) 2011-2012 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

if nargin<4 || isempty(method), method = 'spli'; end

defaultopt = struct(  ...
    'extrapolate' , ~strcmp(method,'cheb'),...
    'order'       , 3);
if nargin<5
  options = defaultopt;
else
  warning('off','catstruct:DuplicatesFound')
  options = catstruct(defaultopt,options);
end

% Define approximation function
interp.fspace      = fundefn(method,n(:),smin(:),smax(:),options.order);

interp.extrapolate = options.extrapolate;

% State collocation nodes
s                  = gridmake(funnode(interp.fspace));
