function [interp,s] = recsinterpinit(n,smin,smax,method,options)
% RECSINTERPINIT Prepares an interpolation structure for RECS
%
% [INTERP,S] = RECSINTERPINIT(N,SMIN,SMAX) creates a structure of interpolation
% in which N designates the order of approximation (if scalar, the same order is
% applied for all dimensions), and SMIN and SMAX are size-d vectors of left and
% right endpoints of the state space. RECSINTERPINIT returns the structure
% INTERP with the following fields:
%    fspace       : a definition structure for the interpolation family
%    Phi          : the basis matrix at collocation nodes
% RECSINTERPINIT returns also the prod(N)-by-d matrix S, which represents the
% state variables on the grid.
%
% [INTERP,S] = RECSINTERPINIT(N,SMIN,SMAX,METHOD) use the string METHOD to
% define the interpolation method, either spline ('spli', default), or Chebyshev
% polynomials ('cheb').
%
% [INTERP,S] = RECSINTERPINIT(N,SMIN,SMAX,METHOD,OPTIONS)

% Copyright (C) 2011-2012 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
if nargin<4 || isempty(method), method = 'spli'; end

if strcmp(method,'lin')
  error('RECS does not support linear interpolation')
end

defaultopt = struct(  ...
    'extrapolate' , ~strcmp(method,'cheb'),...
    'order'       , 3);
if nargin<5
  options = defaultopt;
else
  warning('off','catstruct:DuplicatesFound')
  options = catstruct(defaultopt,options);
end

%% Define approximation function
d                  = length(smin);
if isscalar(n) && d>1, n = n*ones(d,1); end
interp.fspace      = fundefn(method,n(:),smin(:),smax(:),options.order);

% interp.extrapolate = options.extrapolate;

%% Basis matrix
interp.Phi         = funbasx(interp.fspace);

%% State collocation nodes
s                  = gridmake(funnode(interp.fspace));
