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
  [F,J]     = model(x,par{2:end});
  isdomerrJ = any((~isfinite(J) | (imag(J)~=0) | isnan(J)),2);
else
  F         = model(x,par{2:end});
  J         = [];
  isdomerrJ = zeros(size(F));
end

%%  Calculation of number of ex-post domain errors
isdomerrF = ~isfinite(F) | (imag(F)~=0) | isnan(F);
domerr    = sum(any([isdomerrF isdomerrJ],2));
% if domerr>0
%   fprintf(1,'domerr = %i \n',domerr);
% end