function [interp1,s,x] =  recsConvert(interp0,model,order,smin,smax,options)
% RECSCONVERT Converts the interpolation structure of a model to another form

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

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

interp1.fspace = fundefn('spli',order,smin,smax);                 % function space
s              = gridmake(funnode(interp1.fspace));
interp1.Phi    = funbasx(interp1.fspace);

warning('off','RECS:Extrapolation');
[~,x] = recsSimul(model,interp0,s,0,[],options);
warning('on','RECS:Extrapolation');

interp1.cx = funfitxy(interp1.fspace,interp1.Phi,x);

interp1.cz = funconv(interp0.cz,interp0.fspace,interp1.fspace);
if isfield(interp0,'ch')
  interp1.ch = funconv(interp0.ch,interp0.fspace,interp1.fspace);
end
