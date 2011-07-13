function [F,J,domerr] = pathtransform(x,jacflag)
% PATHTRANSFORM Allows to pass a model to the solver path
%
% Path does not accept anonymous function or supplementary arguments, so model
% parameters have to be passed in a global variable

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

global par
model = par{1};

domerr = 0;
if (jacflag)
  [F,J] = model(x,par{2:end});
else
  F     = model(x,par{2:end});
  J     = [];
end

% In test:
isdomerr = ~(isfinite(F)+(imag(F)==0));
if sum(isdomerr)>=1, fprintf(1,'domerr = %i \n',sum(isdomerr)); end