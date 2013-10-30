function [x,s,z,F,exitflag,N] = recsSolveDeterministicPbSP(model,s0,istart,T,xss,zss,sss,options)
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
    'eqsolver'        , 'lmmcp'                          ,...
    'eqsolveroptions' , struct('Diagnostics'    , 'off' ,...
                               'DerivativeCheck', 'on' ,...
                               'Jacobian'       , 'on'));
if nargin <=7
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

functions = model.functions;
nperiods  = model.nperiods;
params    = model.params;
e         = cell(nperiods,1);
for i=1:nperiods
  if istart~=1 && i<=(istart-1)
    e{i} = repmat(model.shocks{i}.w'*model.shocks{i}.e,T-1,1); 
  else
    e{i} = repmat(model.shocks{i}.w'*model.shocks{i}.e,T,1); 
  end
end

inext    = @(iperiod) (iperiod+1).*(iperiod<nperiods)+ones(size(iperiod)).*(iperiod==nperiods);
vec      = @(X) X(:);

%% Dimensions
dim     = model.dim;
Dstart  = sum(cell2mat(dim(inext(istart:nperiods),1)));
Mstart  = sum(cell2mat(dim(istart:nperiods,2)));
Pstart  = sum(cell2mat(dim(istart:nperiods,3)));
D       = sum(cell2mat(dim(:,1)));
M       = sum(cell2mat(dim(:,2)));
P       = sum(cell2mat(dim(:,3)));

%% Maximum number of non-zero elements in the Jacobian
% This is an upper bound: it assumes that istart=1 and all derivatives are non-zero.
nnzJac = 0;
for i=1:nperiods
  nnzJac = nnzJac+T*(...
      dim{i,2}*(dim{i,2}+dim{i,3})+...               % Main diagonal blocks
      dim{i,3}*(1+dim{i,2}+dim{inext(i),1})+...
      dim{inext(i)}*(1+dim{i,2})+...
      dim{i,3}*dim{inext(i),2}+...                   % Superdiagonal blocks
      dim{i,1}*(dim{i,2}+dim{i,3}+dim{inext(i),1})); % Subdiagonal blocks
end
nnzJac = nnzJac...
         -dim{nperiods,3}*dim{1,2}...            % Superdiagonal blocks
         -dim{1,1}*(dim{1,2}+dim{1,3}+dim{2,1}); %   Subdiagonal blocks

%% Bounds
LBper = cell(nperiods,1);
UBper = cell(nperiods,1);
for i=1:nperiods, 
  [LBper{i},UBper{i}] = functions(i).b(sss{i},params);
  LBper{i} = [LBper{i} -inf(1,dim{inext(i),1}+dim{i,3})];
  UBper{i} = [UBper{i} +inf(1,dim{inext(i),1}+dim{i,3})];
end
LB      = cat(2, LBper{:});
UB      = cat(2, UBper{:});
LB      = [cat(2,LBper{istart:end})'; reshape(LB(ones(T-1,1),:)',(T-1)*(M+D+P),1)];
UB      = [cat(2,UBper{istart:end})'; reshape(UB(ones(T-1,1),:)',(T-1)*(M+D+P),1)];

%% First guess equal to steady state
X      = [xss'; zss'; [sss(2:end); sss{1}]'];
Xstart = X(:,istart:end);
X      = cat(2, X{:});
X      = [cat(2,Xstart{:})'; reshape(X(ones(T-1,1),:)',(T-1)*(M+D+P),1)];

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

% Indexes of variables at the first period
nskip = D+M+P-(Dstart+Mstart+Pstart);
ixstart = cellfun(@(dimX) zeros(0,dimX),dim(:,2),'UniformOutput', false);
izstart = cellfun(@(dimX) zeros(0,dimX),dim(:,3),'UniformOutput', false);
isstart = cellfun(@(dimX) zeros(0,dimX),dim(:,1),'UniformOutput', false);
for i=istart:nperiods
  ixstart{i}        = ix{i}-nskip;
  izstart{i}        = iz{i}-nskip;
  isstart{inext(i)} = is{inext(i)}-nskip;
end

% Indexes of variables' position for all the horizon
iX2iXT = @(istart,iX,dimX) [vec(istart);
                    vec((repmat(iX,T-1,1)+(Dstart+Mstart+Pstart)+(D+M+P)*repmat((0:T-2)',1,dimX))')];
ixT = cellfun(iX2iXT,ixstart,ix,dim(:,2),'UniformOutput', false);
izT = cellfun(iX2iXT,izstart,iz,dim(:,3),'UniformOutput', false);
isT = cellfun(iX2iXT,isstart,is,dim(:,1),'UniformOutput', false);

% Indexes of some variables to account for the difference between forward/predetermined and static variables
ixnext = cellfun(@(iX,dimX) vec((repmat(iX,T-1,1)+...
                                 (Dstart+Mstart+Pstart)+...
                                 (D+M+P)*repmat((0:T-2)',1,dimX))'),...
                 ix,dim(:,2),'UniformOutput', false);
izprev = [vec(izstart{4});
          vec((repmat(iz{4},T-2,1)+...
               (Dstart+Mstart+Pstart)+(D+M+P)*repmat((0:T-3)',1,dim{4,3}))')];
iznext = vec((repmat(iz{istart},T-1,1)+...
              (Dstart+Mstart+Pstart)+(D+M+P)*repmat((0:T-2)',1,dim{istart,3}))');

%% Solve deterministic problem
SCPSubProblem = @(X0,S0) runeqsolver(@recsDeterministicPbSP,X0,LB,UB,eqsolver, ...
                                     eqsolveroptions,functions,S0,xss{1},e, ...
                                     params,M,P,D,ixT,izT,isT, ixnext,izprev,iznext, ...
                                     nnzJac,dim,istart);

% Simple continuation problem applied on a Newton solve
[X,F,exitflag,N] = SCP(X,s0,sss{istart},SCPSubProblem,1);

if exitflag~=1
  warning('RECS:FailureDeterministic',...
          'Failure to find the perfect foresight solution');
end

%% Prepare output
X2xzs = @(iX,dimX) reshape(X(iX),dimX,[])';
x = cellfun(X2xzs,ixT,dim(:,2),'UniformOutput', false);
z = cellfun(X2xzs,izT,dim(:,3),'UniformOutput', false);
s = cellfun(X2xzs,isT,dim(:,1),'UniformOutput', false);
s{istart} = [s0; s{istart}];
