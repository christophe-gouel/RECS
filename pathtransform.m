function [F,J,domerr] = pathtransform(x,jacflag)
% PATHTRANSFORM Allows to pass a model to the solver path
%
% Path does not accept anonymous function or supplementary arguments, so model
% parameters have to be passed in a global variable

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

global par

domerr = 0;
if (jacflag)
  [F,J] = feval(par{1},x,par{2:end});
else
  F     = feval(par{1},x,par{2:end});
  J     = [];
end

