function [x,s,z,F,exitflag,N] = recsSolveDeterministicPbSP(model,s0,T,xss,zss,sss,options)
% RECSSOLVEDETERMINISTICPBSP Solves a perfect foresight problem
%
% X = RECSSOLVEDETERMINISTICPBSP(MODEL,S0,T,XSS,ZSS,SSS) tries to find the perfect
% foresight solution of the model defined in the object MODEL. The initial
% values of state variable are provided by the 1-by-d vector S0. Time horizon
% (number of periods) is given by the integer T. XSS, ZSS and SSS are,
% respectively, a 1-by-m vector containing the values of response variables at the
% deterministic steady state, a 1-by-p vector containing the values of expectations
% variables at steady state, and a 1-by-d vector containing the values of the state
% variables at steady state. RECSSOLVEDETERMINISTICPBSP returns X, a T-by-m matrix,
% the value of response variables over the time horizon.
%
% X = RECSSOLVEDETERMINISTICPBSP(MODEL,S0,T,XSS,ZSS,SSS,OPTIONS) solves the problem with the
% parameters defined by the structure OPTIONS. The fields of the structure are
%    eqsolver         : 'fsolve', 'lmmcp' (default), 'ncpsolve' or 'path'
%    eqsolveroptions  : options structure to be passed to eqsolver (default:
%                       empty structure)
%
% [X,S] = RECSSOLVEDETERMINISTICPBSP(MODEL,S0,T,XSS,ZSS,SSS,...) returns S, a T-by-d
% matrix, containing the value of state variables over the time horizon.
%
% [X,S,Z] = RECSSOLVEDETERMINISTICPBSP(MODEL,S0,T,XSS,ZSS,SSS,...) returns Z, a
% T-by-p matrix, containing the value of expectations variables over the time
% horizon.
%
% [X,S,Z,F] = RECSSOLVEDETERMINISTICPBSP(MODEL,S0,T,XSS,ZSS,SSS,...) returns F, a
% T-by-m matrix, containing the values of equilibrium equations over the time
% horizon.
%
% [X,S,Z,F,EXITFLAG] = RECSSOLVEDETERMINISTICPBSP(MODEL,S0,T,XSS,ZSS,SSS,...)
%
% [X,S,Z,F,EXITFLAG,N] = RECSSOLVEDETERMINISTICPBSP(MODEL,S0,T,XSS,ZSS,SSS,...)
%
% See also RECSFIRSTGUESS, RECSSOLVEREE, RECSSS, SCP.

% Copyright (C) 2011-2013 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
defaultopt = struct(                                      ...
    'checkfinalstate' , 0                                ,...
    'eqsolver'        , 'lmmcp'                          ,...
    'eqsolveroptions' , struct('Diagnostics'    , 'off' ,...
                               'DerivativeCheck', 'off' ,...
                               'Jacobian'       , 'off'));
if nargin <=6
  options = defaultopt;
else
  warning('off','catstruct:DuplicatesFound')
  if isfield(options,'eqsolveroptions')
    options.eqsolveroptions = catstruct(defaultopt.eqsolveroptions,options.eqsolveroptions);
  end
  options = catstruct(defaultopt,options);
end
eqsolver         = lower(options.eqsolver);
eqsolveroptions  = options.eqsolveroptions;

dim       = model.dim;
functions = model.functions;
nperiods  = model.nperiods;
params    = model.params;
e         = cell(nperiods,1);
for i=1:nperiods, e{i} = repmat(model.shocks{i}.w'*model.shocks{i}.e,T,1); end

inext    = @(iperiod) (iperiod+1)*(iperiod<nperiods)+1*(iperiod==nperiods);
vec      = @(X) X(:);

D       = sum(cell2mat(dim(:,1)));
M       = sum(cell2mat(dim(:,2)));
P       = sum(cell2mat(dim(:,3)));
n       = T*(D+M+P);

LBx         = cell(nperiods,1);
UBx         = cell(nperiods,1);
LB          = cell(nperiods,1);
UB          = cell(nperiods,1);
for i=1:nperiods, 
  [LBx{i},UBx{i}] = functions(i).b(sss{i},params);
  LB{i} = [LBx{i} -inf(1,dim{inext(i),1}+dim{i,3})];
  UB{i} = [UBx{i} +inf(1,dim{inext(i),1}+dim{i,3})];
end
LB = cat(2,LB{:});
UB = cat(2,UB{:});
LB = reshape(LB(ones(T,1),:)',n,1);
UB = reshape(UB(ones(T,1),:)',n,1);

X  = [xss'; zss'; [sss(2:end); sss{1}]'];
X  = cat(2, X{:});
X  = reshape( X(ones(T,1),:)',n,1);

%% Create indexes of variables' position
% Indexes of variables' position for one period
ix      = cell(nperiods,1);
iz      = cell(nperiods,1);
is      = cell(nperiods,1);
index   = 1;
for i=1:nperiods
  ix{i} = index:(index+dim{i,2}-1);
  index = index+dim{i,2};
  iz{i} = index:(index+dim{i,3}-1);
  index = index+dim{i,3};
  is{inext(i)} = index:(index+dim{inext(i),1}-1);
  index = index+dim{inext(i),1};
end

% Indexes of variables' position for all the horizon
iX2iXT = @(iX,dimX) vec((repmat(iX,T,1)+(D+M+P)*repmat((0:T-1)',1,dimX))');
ixT = cellfun(iX2iXT,ix,dim(:,2),'UniformOutput', false);
izT = cellfun(iX2iXT,iz,dim(:,3),'UniformOutput', false);
isT = cellfun(iX2iXT,is,dim(:,1),'UniformOutput', false);

%% Solve deterministic problem

% Precalculation of indexes for the sparse Jacobian
% tmp      = zeros(M+p+d,M+p+d,T);
% [~,grid] = blktridiag(tmp,tmp(:,:,1:end-1),tmp(:,:,1:end-1));

[X,F,exitflag] = runeqsolver(@recsDeterministicPbSP,X,LB,UB,eqsolver, ...
                             eqsolveroptions,functions,s0,xss{1},e, ...
                             params,M,P,D,ix,iz,is);

% SCPSubProblem = @(X0,S0) runeqsolver(@recsDeterministicPbSP,X0,LB,UB,eqsolver, ...
%                                      eqsolveroptions,functions,S0,xss{1},e, ...
%                                      params,M,P,D,ix,iz,is);

% % Simple continuation problem applied on a Newton solve
%[X,F,exitflag,N] = SCP(X,s0,sss{1},SCPSubProblem,1);
if exitflag~=1
  warning('RECS:FailureDeterministic',...
          'Failure to find the perfect foresight solution');
end

%% Prepare output
X = reshape(X,M+P+D,T)';
x = cellfun(@(iX) X(:,iX),ix,'UniformOutput', false);
z = cellfun(@(iX) X(:,iX),iz,'UniformOutput', false);
s = cellfun(@(iX) X(:,iX),is,'UniformOutput', false);

s{1} = [s0; s{1}];
if ~isempty(F), F = reshape(F,M+P+D,T)'; end

N = [];
%% Check is the final state is self-replicating (meaning that sT=sss and xT=xss)
% if options.checkfinalstate
%   deltafinal = max(abs([s(end,:)-sss x(end,:)-xss]));
%   if deltafinal>sqrt(eps)
%     warning('RECS:FinalState',...
%             'Final state is not self-replicating (max(|delta|)=%g)',...
%             deltafinal);
%   end
% end
