function [interp1,s,x] =  recsConvert(interp0,model,order,smin,smax,options)
% RECSCONVERT Converts the interpolation structure of a model to another form
%
% INTERP1 = RECSCONVERT(INTERP0,MODEL) convert the interpolation structure
% INTERP0 corresponding to the model defined in the structure MODEL. It returns
% a new interpolation structure INTERP1.
%
% INTERP1 = RECSCONVERT(INTERP0,MODEL,ORDER) changes the order of the
% interpolation structure using ORDER as the new order.
%
% INTERP1 = RECSCONVERT(INTERP0,MODEL,ORDER,SMIN) changes the minimum values of
% the interpolation grid using SMIN as the new minimum.
%
% INTERP1 = RECSCONVERT(INTERP0,MODEL,ORDER,SMAX) changes the maximum values of
% the interpolation grid using SMAX as the new maximum.
%
% INTERP1 = RECSCONVERT(INTERP0,MODEL,ORDER,SMAX,OPTIONS) converts the
% interpolation structure with the parameters defined by the structure OPTIONS.
% The fields of the structure are those used in recsSimul.
%
% [INTERP1,S] = RECSCONVERT(INTERP0,MODEL,...) returns the matrix S containing
% the new grid of state variables.
%
% [INTERP1,S,X] = RECSCONVERT(INTERP0,MODEL,...) returns the matrix X containing
% the value of response variables on the new grid.
%
% See also RECSSIMUL.

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
if nargin < 6
  options = struct([]);
  if nargin < 5
    smax = [];
    if nargin < 4
      smin = [];
    end
  end
end

if isempty(order), order = interp0.fspace.n'; end
if isempty(smax),  smax  = interp0.fspace.b'; end
if isempty(smin),  smin  = interp0.fspace.a'; end

defaultopt = struct(...
    'simulmethod' , 'solve',...
    'stat'        , 0);

options = catstruct(options,defaultopt);

%% Define new interpolation structure
interp1.fspace = fundefn('spli',order,smin,smax);
s              = gridmake(funnode(interp1.fspace));
interp1.Phi    = funbasx(interp1.fspace);

%% Simulate the model on the new grid
warning('off','RECS:Extrapolation');
[~,x] = recsSimul(model,interp0,s,0,[],options);
warning('on','RECS:Extrapolation');

%% Calculate the new interpolation coefficients
interp1.cx = funfitxy(interp1.fspace,interp1.Phi,x);

interp1.cz = funconv(interp0.cz,interp0.fspace,interp1.fspace);
if isfield(interp0,'ch')
  interp1.ch = funconv(interp0.ch,interp0.fspace,interp1.fspace);
end
