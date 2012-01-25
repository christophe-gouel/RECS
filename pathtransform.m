function [F,J,domerr] = pathtransform(x,jacflag)
% PATHTRANSFORM Allows to pass a model to the solver path
%
% Path does not accept anonymous function or supplementary arguments, so model
% parameters have to be passed in a global variable

% Copyright (C) 2011-2012 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

global par
model = par{1};

if (jacflag)
  [F,J] = model(x,par{2:end});
else
  F     = model(x,par{2:end});
  J     = [];
end

%%  Calculation of number of domain errors
domerr = 0;
% Domain errors should be calculated from the situation when x values exceed
% domain definition of the functions used. Instead, for simplicity, here the
% number of domain errors is calculated from the output: the number of infinite,
% complex or NaN numbers in F:
isdomerr = ~isfinite(F) | (imag(F)~=0) | isnan(F);
if any(isdomerr)
%  fprintf(1,'domerr = %i \n',sum(isdomerr));
  domerr = sum(isdomerr);
end