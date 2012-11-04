function [F,J,domerr] = pathtransform(x,jacflag)
% PATHTRANSFORM Allows to pass a model to the solver PATH
%
% The MCP solver PATH does not accept anonymous function or supplementary
% arguments, so the function PATHTRANSFORM is used as model argument for PATH
% and the original model, along with its supplementary arguments, is passed as
% an anonymous function in a global variable (eqtosolve).

% Copyright (C) 2011-2012 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

global eqtosolve

if (jacflag)
  [F,J]     = eqtosolve(x);
  isdomerrJ = any((~isfinite(J) | (imag(J)~=0) | isnan(J)),2);
else
  F         = eqtosolve(x);
  J         = [];
  isdomerrJ = zeros(size(F));
end

%%  Calculation of number of ex-post domain errors
isdomerrF = ~isfinite(F) | (imag(F)~=0) | isnan(F);
domerr    = sum(any([isdomerrF isdomerrJ],2));
