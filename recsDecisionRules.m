function result = recsDecisionRules(model,interp,states,s0,space,options)
% RECSDECISIONRULES Plots a model decision rules
%
% RESULT = RECSDECISIONRULES(MODEL,INTERP) plots a model decision rules
% successively with respect to each state variables, while holding others fixed
% at their steady-state values. The state space of interpolation is used to
% determine the range of variation and 100 points are used in each
% dimension. RECSDECISIONRULES produces a figure by state variable, with as many
% suplots as there are control variables. RESULT is a 1xd structure array with
% fields s and x, respectively the state variables values used to draw the
% decision rules and the control variables values.
%
% RESULT = RECSDECISIONRULES(MODEL,INTERP,STATES) plots the decision rules with
% respect to the state variables indexed by the vector STATES.
%
% RESULT = RECSDECISIONRULES(MODEL,INTERP,STATES,S0) plots the decision rules
% with respect to the state variables defined in STATES while holding the other
% state variables fixed at the values defined in S0. S0 is an n-by-d matrix, n
% indexing the different combinations of state variables on which the decision
% rules are plotted. If n>1, the different decision rules are plotted on the
% same figures.
%
% RESULT = RECSDECISIONRULES(MODEL,INTERP,STATES,S0,SPACE) plots the decision
% rules on the intervals and for a number of points defined in SPACE. SPACE is a
% matrix in which the first row contains the lower bounds of the intervals, the
% second row contains the upper bound, and the third row the number of points to
% use in each dimension. The columns correspond to the state variables indexed
% in STATES.
%
% RESULT = RECSDECISIONRULES(MODEL,INTERP,STATES,S0,SPACE,OPTIONS) plots the
% decision rules with the parameters defined by the structure OPTIONS. The
% fields of the structure are those used in recsSimul.
%
% See also RECSSIMUL.

% Copyright (C) 2011-2016 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
[d,m] = model.dim{1:2};
if nargin<3 || isempty(states), states = 1:d; end
if nargin<4 || isempty(s0)
  s0 = model.sss;
else
  validateattributes(s0,{'numeric'},{'size',[NaN,d]},4)
end
if nargin<5 || isempty(space)
  range = [interp.fspace.a(states)' interp.fspace.b(states)'];
  n     = 100*ones(1,length(states));
else
  validateattributes(space,{'numeric'},{'size',[3,length(states)]},5)
  range = space(1:2,:)';
  n     = space(3,:);
end
overridingopt = struct(...
    'accuracy', 0,...
    'stat'    , 0);
if nargin<6
  options = overridingopt;
else
  options = catstruct(options,overridingopt);
end

%% Calculate, plot and export decision rules
for i=1:length(states)
  %% Calculate DR
  s              = repmat(s0,1,n(i));
  s              = reshape(s',d,[])';
  s(:,states(i)) = repmat(linspace(range(i,1),range(i,2),n(i)),...
                          1,...
                          size(s0,1));

  [~,x] = recsSimul(model,interp,s,1,[],options);

  x = reshape(x,n(i),[],m);
  x = permute(x,[1 3 2]);

  s = reshape(s,n(i),[],d);
  s = permute(s,[1 3 2]);

  %% Plot DR
  figure
  for j=1:m
    subplot(ceil((m)/ceil(sqrt(m))),ceil(sqrt(m)),j)
    plot(s(:,states(i),1),squeeze(x(:,j,:)))
    xlabel(model.symbols.states{i});
    ylabel(model.symbols.controls{j});
  end

  %% Export DR
  result(i).s = s; %#ok<AGROW>
  result(i).x = x; %#ok<AGROW>
end
