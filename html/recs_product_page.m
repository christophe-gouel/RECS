%% RECS toolbox Version 0.5.1
% RECS is a MATLAB solver for dynamic, stochastic, rational expectations
% equilibrium models. RECS stands for "Rational Expectations Complementarity
% Solver". This name emphasizes that RECS has been developed specifically to
% solve models that include complementarity equations, also known as models with
% occasionally binding constraints.
%
% RECS is designed to solve small-scale nonlinear and complementarity models,
% but not large-scale models. For solving large-scale problems, but without
% complementarity equations, see <http://www.dynare.org/ Dynare> or similar
% toolboxes.
%
% <html>
% <h2>Available documentation</h2>
% </html>
%
% RECS documentation is not complete, but is sufficient to build standard
% models (many examples are provided in <demos.html Demos>):
%
% * <getting_started.html Getting started>
% * <user_guide.html User's Guide>
% * <recs_functions.html Functions>
% * <demos.html RECS Toolbox Demos>
%
% Release notes are available at <https://github.com/christophe-gouel/RECS/wiki/Release-Notes>.
%
% <html>
% <h2>License</h2>
% </html>
%
% Unless stated otherwise, all files in the RECS toolbox are licensed using the
% Expat license, a permissive free software license. Please see the <LICENSE.txt
% software license> for more information.
%
% <html>
% <h2>Source</h2>
% </html>
%
% RECS source code and development version is available at
% <https://github.com/christophe-gouel/recs>.
%
% Bugs should be reported at <https://github.com/christophe-gouel/RECS/issues>.
%
% <html>
% <h2>Acknowledgments</h2>
% </html>
%
% This solver started as a reimplementation of the solvers |remsolve| and
% |resolve| from Miranda and Fackler (2002) and Fackler (2005). RECS would not
% exist without these earlier contributions. RECS benefited also of many inputs
% from <http://www.mosphere.fr/ Pablo Winant>, especially with respect to
% models' parsing.
%
% This work was generously supported by
%
% * the European Union's Seventh Framework Programme FP7/2007-2011 under Grant
%   Agreements #212036 AgFoodTRAde and #290693 FOODSECURE (for solver development);
% * the Knowledge for Change (KCP) Trust Fund;
% * the AGRODEP Consortium (for documentation writing).
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
