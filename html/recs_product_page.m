%% RECS toolbox
% Copyright (C) 2011-2012 by <http://www.christophegouel.com Christophe
% Gouel>. RECS is released under a free software license, please see the
% <LICENSE.txt software license> for more information. Source code is available
% at <https://github.com/christophe-gouel/recs>.
%
% RECS is a MATLAB solver for dynamic, stochastic, rational expectations
% equilibrium models. RECS stands for Rational Expectations Complementarity
% Solver, which emphasizes the fact that RECS is specifically developed to solve
% models that include complementarity equations, also known as models with
% occasionally binding constraints.
%
% RECS is designed to solve small-scale nonlinear and complementarity models,
% but not large-scale models. For solving large-scale problems, but without
% complementarity equations, see <http://www.dynare.org/ Dynare>.
%
% <html>
% <h2>Available documentation</h2>
% </html>
%
% RECS documentation is still not complete in particular with respect to
% advanced features, but it is sufficient to build standard models (many
% examples are provided in Demos).
%
% * <getting_started.html Getting started>
% * <user_guide.html User's Guide>
% * <recs_functions.html Functions>
% * <demos.html RECS Toolbox Demos>
%
% <html>
% <h2>Acknowledgments</h2>
% </html>
%
% This solver started as a reimplementation of the solvers |remsolve| and
% |resolve| from Miranda and Fackler (2002) and Fackler (2005). RECS would not
% exist without these earlier contributions. RECS benefited also of many inputs
% from Pablo Winant, especially with respect to models' parsing. This work was
% generously supported by the AgFoodTrade project, funded under the Seventh
% Framework Programme for Research and Development, DG-Research, European
% Commission.
%
% <html>
% <h2>References</h2>
% </html>
%
% <http://dx.doi.org/10.1007/s10614-005-1784-z Fackler, P. L. (2005). A MATLAB
% Solver for Nonlinear Rational Expectations Models. _Computational Economics_,
% 26(2), 173-181.>
%
% Miranda, M. J. and Fackler, P. L. (2002). _Applied Computational Economics and
% Finance_. Cambridge: MIT Press.
