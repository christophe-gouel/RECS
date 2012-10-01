function [x,s,z,F] = recsSolveDeterministicPb(model,s0,T,xss,zss,sss,options)
% RECSSOLVEDETERMINISTICPB Solves a perfect foresight problem
%
% X = RECSSOLVEDETERMINISTICPB(MODEL,S0,T,XSS,ZSS,SSS) tries to find the perfect
% foresight solution of the model defined in the structure MODEL. The initial
% values of state variable are provided by the 1-by-d vector S0. Time horizon
% (number of periods) is given by the integer T. XSS, ZSS and SSS are,
% respectively, a 1-by-m vector containing the values of response variables at the
% deterministic steady state, a 1-by-p vector containing the values of expectations
% variables at steady state, and a 1-by-d vector containing the values of the state
% variables at steady state. RECSSOLVEDETERMINISTICPB returns X, a T-by-m matrix,
% the value of response variables over the time horizon.
%
% X = RECSSOLVEDETERMINISTICPB(MODEL,S0,T,XSS,ZSS,SSS,OPTIONS) solves the problem with the 
% parameters defined by the structure OPTIONS. The fields of the structure are
%    eqsolver         : 'fsolve', 'lmmcp' (default), 'ncpsolve' or 'path'
%    eqsolveroptions  : options structure to be passed to eqsolver (default:
%                       empty structure)
%
% [X,S] = RECSSOLVEDETERMINISTICPB(MODEL,S0,T,XSS,ZSS,SSS,...) returns S, a T-by-d
% matrix, containing the value of state variables over the time horizon.
%
% [X,S,Z] = RECSSOLVEDETERMINISTICPB(MODEL,S0,T,XSS,ZSS,SSS,...) returns Z, a
% T-by-p matrix, containing the value of expectations variables over the time
% horizon.
%
% [X,S,Z,F] = RECSSOLVEDETERMINISTICPB(MODEL,S0,T,XSS,ZSS,SSS,...) returns F, a
% T-by-m matrix, containing the values of equilibrium equations over the time
% horizon.
%
% See also RECSFIRSTGUESS, RECSSOLVEREE, RECSSS.

% Copyright (C) 2011-2012 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization

defaultopt = struct(                                     ...
    'eqsolver'         , 'lmmcp'                        ,...
    'eqsolveroptions'  , struct('DerivativeCheck','off'));
if nargin <=6
  options = defaultopt; 
else
  warning('off','catstruct:DuplicatesFound')
  options = catstruct(defaultopt,options);
end
eqsolver         = lower(options.eqsolver);
eqsolveroptions  = options.eqsolveroptions;

m = size(xss,2);
p = size(zss,2);
d = size(sss,2);
n = T*(m+p+d);

e = model.w'*model.e;
if isa(model.func,'char')
  func = str2func(model.func);
elseif isa(model.func,'function_handle')
  func = model.func;
else
  error('model.func must be either a string or a function handle')
end
params = model.params;

ix = [sum(numjac(@(S) Bounds(func,S,params,1,true(m,2)),sss)~=0,2,'native') ...
      sum(numjac(@(S) Bounds(func,S,params,2,true(m,2)),sss)~=0,2,'native')];
nx = int16(sum(ix,1));
M  = m+sum(nx);

n = n+T*sum(nx);

w = zeros(1,nx(1));
v = zeros(1,nx(2));

[LBx,UBx] = mcptransform(func,'b',sss,[],[],[],[],[],params,[],ix,nx);
LB = [LBx -inf(1,p+d)];
UB = [UBx +inf(1,p+d)];

LB = reshape(LB(ones(T,1),:)',n,1);
UB = reshape(UB(ones(T,1),:)',n,1);

X = [xss w v zss sss];
X = X(ones(T,1),:)';
X = reshape(X,n,1);

%% Solve deterministic problem
% Simple continuation problem applied on a Newton solve

SCPSubProblem = @(X0,S0) runeqsolver(@recsDeterministicPb,X0,LB,UB,eqsolver,...
                                     eqsolveroptions,func,S0,xss,p,e,params,ix,nx);

[X,F,exitflag] = SCP(X,s0,sss,SCPSubProblem,1);
if exitflag~=1
  warning('RECS:FailureDeterministic',...
          'Failure to find the perfect foresight solution');
end

%% Prepare output
X = reshape(X,M+p+d,T)';
x = X(:,1:m);
z = X(:,(M+1):(M+p));
s = [s0; X(1:end-1,(M+p+1):(M+p+d))];
if ~isempty(F), F = reshape(F,M+p+d,T)'; end
