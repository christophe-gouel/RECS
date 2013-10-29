function F = recsDeterministicPbSP(X,functions,s0,xss,e,params,M,P,D,ix,iz,is)
% RECSDETERMINISTICPBSP Evaluates the equations and Jacobian of the deterministic problem.
%
% RECSDETERMINISTICPBSP is called by recsSolveDeterministicPbSP. It is not meant to be
% called directly by the user.
%
% See also RECSSOLVEDETERMINISTICPBSP.

% Copyright (C) 2011-2013 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
n        = size(X,1);
T        = n/(M+P+D);
nperiods = length(e);

inext    = @(iperiod) (iperiod+1)*(iperiod<nperiods)+1*(iperiod==nperiods);

X     = reshape(X,M+P+D,T)';
x     = cellfun(@(iX) X(:,iX),ix,'UniformOutput', false);
z     = cellfun(@(iX) X(:,iX),iz,'UniformOutput', false);
snext = cellfun(@(iX) X(:,iX),is,'UniformOutput', false);

s        = snext;
s{1}     = [s0; s{1}(1:(end-1),:)];
xnext    = x;
xnext{1} = [x{1}(2:end,:); xss];

f = cell(nperiods,1);
h = cell(nperiods,1);
g = cell(nperiods,1);

%% Computation of equations and Jacobian
for i=1:nperiods
  f{i} = functions(i).f(s{i},x{i},z{i},params,[1 0 0 0]);
  h{i} = z{i}-functions(i).h(s{i},x{i},e{i},snext{inext(i)},xnext{inext(i)},params,[1 0 0 0 0 0]);
  g{i} = snext{inext(i)}-functions(i).g(s{i},x{i},e{i},params,[1 0 0 0]);
end

F = [f h g]';
F = cat(2,F{:});
F = reshape(F',n,1);