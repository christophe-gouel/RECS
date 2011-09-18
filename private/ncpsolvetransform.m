function [F,J] = ncpsolvetransform(x,model,varargin)
% NCPSOLVETRANSFORM Allows to pass a model to the solver ncpsolve
%
% With ncpsolve, equations are of opposite sign with respect to complementary
% variables

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%%
if nargout==2
  [F,J] = model(x,varargin{:});
  J     = -J;
else
  F     = model(x,varargin{:});
end
F       = -F;
