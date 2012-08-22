%% Writing RECS model files

%% Structure of a rational expectation models
% Before starting to write the model file, you need to organize your equations
% in the way described in <def_sre.html Definition of a stochastic rational
% expectations problem>. This means that you need to identify and separate in
% your model the three following group of equations:
%
% *Equilibrium equations*
%
% $$\underline{x}(s) \le x \le \overline{x}(s) \perp f(s,x,z).$$
%
% *Expectations definition*
%
% $$z = \mathrm{E}_{\dot e} \left[h(s,x,\dot e,\dot s,\dot x)\right].$$
%
% *Transition equations*
%
% $$\dot s = g(s,x,\dot e).$$
%
% When defining your model, it is important to try to minimize the number of
% state variables. Currently, RECS is designed to solve small-scale models,
% which means that with more than 3 or 4 state variables you are in unknown
% territory. One way of doing this is, when it is possible, to sum together the
% predetermined variables than can be summed. You can see in demonstrations
% files (e.g., <cs1.html cs1> or <sto1.html sto1>) that states variables are
% often defined as sums of other predetermined variables.

%% Writing a model the easy way
% It is possible to write the model in a way that is quite similar to the
% original mathematical notations (as is proposed in most algebraic modeling
% languages). The model file that has to be created is called a yaml file
% (because it is written in <http://yaml.org YAML>). A yaml file is easily
% readable by human, but not by Matlab. So the file needs to be processed for
% Matlab to be able read it. This is done by the function
% <matlab:doc('recsmodelinit') |recsmodelinit|>, which calls a Python
% preprocessor, dolo-recs, to do the job.
%
% *Yaml file structure*
%
%  declarations:
%
%    states:
%
%    controls:
%
%    expectations:
%
%    shocks:
%
%    parameters:
%
%  equations:
%    
%    arbitrage:
%    
%    transition:
%
%    expectation:
%
%  calibration:
%
%    parameters:
%
%    steady_state:
%
% See <cs1.yaml> or <gro1.yaml> for examples of complete yaml files.
%
% *Yaml syntax conventions*
% 
% For the file to be readable, it is necessary to respect some syntax
% conventions:
%
% * *Indentation:* Inside each block the elements of same scope have to be
% indented using the same number of spaces (do not use tabs). For the above
% example, inside the |declarations| block, |states|, |controls|,
% |expectations|, |shocks|, and |parameters| have to be indented at the same
% level.
%
% * *Comments:* Comments are introduced by the character |#|.
%
% * *Lead/Lag:* Lead variables are indicated |X(1)| and lag variables |X(-1)|.
%
% * *Timing convention:* The following timing conventions are to be
% respected. Transition equations are written by defining the new state variable
% at current time as a function of lagged response and state variables: $s_t =
% g(s_{t-1},x_{t-1},e_t).$
%
% * Do not use |lambda| as a variable/parameter name, this is a restricted word.
%
% * Write equations in the same order than the order of variables declaration.
%
% * Yaml files are case sensitive.
%

%% Writing a model the hard way
%
% <html>
% <hr/>
% </html>
%
%%
% Copyright (C) 2011-2012 Christophe Gouel
%
% Licensed under the Expat license, see <LICENSE.txt>
