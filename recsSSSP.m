function [s,x,z,exitflag] = recsSSSP(model,s,x,options)
% RECSSSSP Solves for the deterministic steady state in rational expectations models
%
% S = RECSSSSP(MODEL) tries to find the non-stochastic steady state of the model
% defined in the object MODEL. This function call uses as first guess for
% steady-state state and response variable the values available in the
% properties sss and xss of the object MODEL. RECSSS returns the value of the
% state variables at steady state. MODEL is an object created by recsmodel.
%
% S = RECSSSSP(MODEL,S) uses the vector S as first guess for steady-state state
% variables.
%
% S = RECSSSSP(MODEL,S,X) uses the vector X as first guess for steady-state response
% variables.
%
% S = RECSSSSP(MODEL,S,X,OPTIONS) solves the problem with the parameters
% defined by the structure OPTIONS. The fields of the structure are
%    display          : 1 to display the steady state if found (default: 1)
%    eqsolver         : 'fsolve', 'lmmcp' (default), 'ncpsolve' or 'path'
%    eqsolveroptions  : options structure to be passed to eqsolver
%
% [S,X] = RECSSSSP(MODEL,...) returns the value of the response
% variables at steady state.
%
% [S,X,Z] = RECSSSSP(MODEL,...) returns the value of the
% expectations variable at steady state.
%
% [S,X,Z,EXITFLAG] = RECSSSSP(MODEL,...) returns EXITFLAG,
% which describes the exit conditions. Possible values are
%    1 : RECSSS converges to the deterministic steady state
%    0 : Failure to converge

% Copyright (C) 2011-2013 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
defaultopt = struct(...
    'display'         , 1                               ,...
    'eqsolver'        , 'lmmcp'                         ,...
    'eqsolveroptions' , struct('Diagnostics'    , 'off',...
                               'DerivativeCheck', 'off',...
                               'Jacobian'       , 'off'));
if nargin<4
  options = defaultopt;
else
  warning('off','catstruct:DuplicatesFound')
  if isfield(options,'eqsolveroptions')
    options.eqsolveroptions = catstruct(defaultopt.eqsolveroptions,options.eqsolveroptions);
  end
  options = catstruct(defaultopt,options);
end
eqsolver        = lower(options.eqsolver);
eqsolveroptions = options.eqsolveroptions;

nperiods  = model.nperiods;
params    = model.params;
dim       = model.dim;
functions = model.functions;
e         = cell(nperiods,1);
for i=1:nperiods, e{i} = model.shocks{i}.w'*model.shocks{i}.e; end

is      = cell(nperiods,1);
ix      = cell(nperiods,1);
is{1}   = 1:dim{1,1};
for i=2:nperiods
  is{i} = (is{i-1}(end)+1):sum(cell2mat(dim(1:i,1)));
end
ix{1} = (is{nperiods}(end)+1):(sum(cell2mat(dim(:,1)))+dim{1,2});
for i=2:nperiods
  ix{i} = (ix{i-1}(end)+1):(sum(cell2mat(dim(:,1)))+sum(cell2mat(dim(1:i,2))));
end

inext    = @(iperiod) (iperiod+1)*(iperiod<nperiods)+1*(iperiod==nperiods);

%% Solve for the deterministic steady state
X  = cat(2,s{:},x{:})';

LB = cell(nperiods,1);
UB = cell(nperiods,1);
for i=1:nperiods
  [LB{i},UB{i}] = functions(i).b(s{i},params);
end
LB = [-inf(size(cat(2,s{:}),2),1); cat(2,LB{:})'];
UB = [+inf(size(cat(2,s{:}),2),1); cat(2,UB{:})'];

[X,~,exitflag] = runeqsolver(@SSResidual,X,LB,UB,eqsolver,eqsolveroptions,...
                             functions,params,e,is,ix,inext,nperiods);

if exitflag~=1
  warning('RECS:SSNotFound','Failure to find a deterministic steady state');
end


%% Prepare outputs
s0     = s;
x0     = x;
s      = mat2cell(X(cat(2,is{:}))',1,cell2mat(dim(:,1)));
x      = mat2cell(X(cat(2,ix{:}))',1,cell2mat(dim(:,2)));
z = cell(nperiods,1);
for i=1:nperiods
 z{i} = functions(i).h(X(is{i})',X(ix{i})',e{i},X(is{inext(i)})',X(ix{inext(i)})',params);
end

%% Display steady state
if exitflag==1 && options.display==1
  deltass = max(abs(cat(2,s{:},x{:})-cat(2,s0{:},x0{:})));
  if deltass<sqrt(eps)
    fprintf(1,'Deterministic steady state (equal to first guess)\n')
  else
    fprintf(1,['Deterministic steady state (different from first guess, ' ...
               'max(|delta|)=%g)\n'],deltass)
  end
  fprintf(1,' State variables:\n\t\t')
  fprintf(1,'%0.4g\t',cat(2,s{:}))
  fprintf(1,'\n\n Response variables:\n\t\t')
  fprintf(1,'%0.4g\t',cat(2,x{:}))
  fprintf(1,'\n\n Expectations variables:\n\t\t')
  fprintf(1,'%0.4g\t',cat(2,z{:}))
  fprintf(1,'\n\n')
end


function F = SSResidual(X,functions,params,e,is,ix,inext,nperiods)
%% SSRESIDUAL evaluates the equations and Jacobians of the steady-state finding problem

f = cell(nperiods,1);
g = cell(nperiods,1);
s = cell(nperiods,1);
x = cell(nperiods,1);
z = cell(nperiods,1);
Fs = cell(nperiods,1);
for i=1:nperiods
  s{i} = X(is{i})';
  x{i} = X(ix{i})';
end

for i=1:nperiods
 z{i} = functions(i).h(s{i},x{i},e{i},s{inext(i)},x{inext(i)},params);
 g{i} = functions(i).g(s{i},x{i},e{i},params);
 f{i} = functions(i).f(s{i},x{i},z{i},params);
end

for i=1:nperiods
  Fs{i} = s{inext(i)}-g{i};
end

F = [cat(2,Fs{:}) cat(2,f{:})]';

return
