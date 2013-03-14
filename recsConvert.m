function [interp1,s,x] =  recsConvert(interp0,model,n,smin,smax,method,options)
% RECSCONVERT Converts the interpolation structure of a model to another form
%
% INTERP1 = RECSCONVERT(INTERP0,MODEL) convert the interpolation structure INTERP0
% corresponding to the model defined in the structure MODEL. It returns a new
% interpolation structure INTERP1. Without any more arguments this command does not
% change the interpolation structure. See below for entering new parameters.
%
% INTERP1 = RECSCONVERT(INTERP0,MODEL,N) changes the number of interpolation
% points along each dimension using N as the new number of grid points in each
% dimension. If N is a scalar, the same dimension is applied to every dimension.
%
% INTERP1 = RECSCONVERT(INTERP0,MODEL,N,SMIN) changes the minimum values of
% the interpolation grid using SMIN as the new minimum. The input N can be
% left empty (N=[]) is you want to change only the minimum values.
%
% INTERP1 = RECSCONVERT(INTERP0,MODEL,N,SMIN,SMAX) changes the maximum values
% of the interpolation grid using SMAX as the new maximum.
%
% INTERP1 = RECSCONVERT(INTERP0,MODEL,N,SMIN,SMAX,METHOD) changes the
% interpolation method, using the string METHOD.
%
% INTERP1 = RECSCONVERT(INTERP0,MODEL,N,SMIN,SMAX,METHOD,OPTIONS) converts the
% interpolation structure with the parameters defined by the structure OPTIONS.
% Most of the fields of the structure are those used in recsSimul, but it is
% also possible to define the following field:
%   order       : for a spline interpolation, it is the order of the spline
%                 (default: 3)
%   simulmethod : simulation method used to conert decision rules to the new
%                 interpolation structure, 'interpolation' or 'solve' (default)
%
% [INTERP1,S] = RECSCONVERT(INTERP0,MODEL,...) returns the matrix S containing
% the new grid of state variables.
%
% [INTERP1,S,X] = RECSCONVERT(INTERP0,MODEL,...) returns the matrix X containing
% the value of response variables on the new grid.
%
% See also RECSINTERPINIT, RECSSIMUL.

% Copyright (C) 2011-2013 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
if nargin<3 || isempty(n), n = interp0.fspace.n'; end
if nargin<4 || isempty(smin), smin = interp0.fspace.a'; end
if nargin<5 || isempty(smax), smax = interp0.fspace.b'; end
if nargin<6 || isempty(method), method = unique(interp0.fspace.bastype); end

defaultopt = struct(        ...
    'order'       , 3      ,...
    'simulmethod' , 'solve');
overridingopt = struct(...
    'accuracy'    , 0 ,...
    'stat'        , 0);
warning('off','catstruct:DuplicatesFound')
if nargin < 7
  options = defaultopt;
else
  options = catstruct(defaultopt,options);
end
options = catstruct(options,overridingopt);

if isscalar(n) && interp0.fspace.d>1
  n = n*ones(interp0.fspace.d,1); 
end

%% Define new interpolation structure
interp1.fspace = fundefn(method{1},n(:),smin(:),smax(:),options.order);
s              = gridmake(funnode(interp1.fspace));
interp1.Phi    = funbasx(interp1.fspace);

%% Simulate the model on the new grid
warning('off','RECS:Extrapolation');
[~,x] = recsSimul(model,interp0,s,1,[],options);
warning('on','RECS:Extrapolation');

%% Calculate the new interpolation coefficients
interp1.cx = funfitxy(interp1.fspace,interp1.Phi,x);

interp1.cz = funconv(interp0.cz,interp0.fspace,interp1.fspace);
if isfield(interp0,'ch')
  interp1.ch = funconv(interp0.ch,interp0.fspace,interp1.fspace);
end
