function result = recsDecisionRules(model,interp,states,space,s0,options)
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
% RESULT = RECSDECISIONRULES(MODEL,INTERP,STATES)
%
% RESULT = RECSDECISIONRULES(MODEL,INTERP,STATES,SPACE)
%
% RESULT = RECSDECISIONRULES(MODEL,INTERP,STATES,SPACE,S0)
%
% RESULT = RECSDECISIONRULES(MODEL,INTERP,STATES,SPACE,S0,OPTIONS)
%
% See also RECSSIMUL.

% Copyright (C) 2011-2013 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
if nargin<3 || isempty(states), states = 1:interp.fspace.d; end
if nargin<4 || isempty(space)
  range = [interp.fspace.a(states)' interp.fspace.b(states)'];
  n     = 100*ones(1,length(states));
else
  range = space(1:2,:);
  n     = space(3,:);
end
if nargin<5 || isempty(s0), s0 = model.sss; end
overridingopt = struct(...
    'accuracy', 0,...
    'stat'    , 0);
if nargin<6
  options = overridingopt;
else
  warning('off','catstruct:DuplicatesFound')
  options = catstruct(options,overridingopt);
end

d = size(s0,2);

for i=1:length(states)
  s              = repmat(s0,1,n(i));
  s              = reshape(s',d,[])';
  s(:,states(i)) = repmat(linspace(range(i,1),range(i,2),n(i)),...
                          1,...
                          size(s0,1));

  [~,x] = recsSimul(model,interp,s,0,[],options);

  m = size(x,2);
  x = reshape(x,n(i),[],m);
  x = permute(x,[1 3 2]);
  
  s = reshape(s,n(i),[],d);
  s = permute(s,[1 3 2]);
  figure
  for j=1:m
    subplot(ceil((m)/ceil(sqrt(m))),ceil(sqrt(m)),j)
    plot(s(:,states(i),1),squeeze(x(:,j,:)))
  end
  
  result(i).s = s;
  result(i).x = x;
end
