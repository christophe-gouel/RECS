function [interp,s] = recsinterpinit(n,smin,smax,method,options)
% RECSINTERPINIT Prepares an interpolation structure for RECS
%
% INTERP = RECSINTERPINIT(N,SMIN,SMAX) creates a structure of interpolation
% in which N designates the order of approximation (if scalar, the same order is
% applied for all dimensions), and SMIN and SMAX are size-d vectors of left and
% right endpoints of the state space. RECSINTERPINIT returns the structure
% INTERP with the following fields:
%    fspace       : a definition structure for the interpolation family
%    Phi          : the basis matrix at collocation nodes
%    s            : a prod(N)-by-d matrix that represents the state variables
%                   on the grid.
%
% INTERP = RECSINTERPINIT(N,SMIN,SMAX,METHOD) use the string METHOD to define
% the interpolation method, either spline ('spli', default), linear
% interpolation ('lin'), or Chebyshev polynomials ('cheb').
%
% INTERP = RECSINTERPINIT(N,SMIN,SMAX,METHOD,OPTIONS) uses the structure OPTIONS
% to create the interpolation structure. The fields of the structure are order :
% for a spline interpolation, it is the order of the spline (default: 3)
%
% [INTERP,S] = RECSINTERPINIT(N,SMIN,SMAX,...) returns the prod(N)-by-d matrix
% S, which represents the state variables on the grid.

% Copyright (C) 2011-2014 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
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

%% Define approximation function
d                  = length(smin);
if isscalar(n) && d>1, n = n*ones(d,1); end
interp.fspace      = fundefn(method,n(:),smin(:),smax(:),options.order);

% interp.extrapolate = options.extrapolate;

%% Basis matrix
interp.Phi         = funbasx(interp.fspace);

%% State collocation nodes
interp.s           = gridmake(funnode(interp.fspace));
if nargout==2, s   = interp.s; end              
