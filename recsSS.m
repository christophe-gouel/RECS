function [s,x,z] = recsSS(model,s,x,options)
% RECSSS Solves for the deterministic steady state in rational expectations models
%
% RECSSS cannot find the deterministic steady state of a functional
% equation problem since, in this case, the steady state depends on
% the model solution. It is not possible either to solve models in which bounds
% depend on state variables.
%
% S = RECSSS(MODEL,S,X) tries to find the non-stochastic steady
% state of the model defined in the structure MODEL, by using as
% first guess the vector of state and response variables S and
% X. RECSSS returns the value of the state variables at steady state.
% MODEL is a structure, which has to include the following fields:
%    [e,w] : discrete distribution with finite support with e the values and w the
%            probabilities (it could be also the discretisation of a continuous
%            distribution through quadrature or Monte Carlo drawings)
%    func   : function name or anonymous function that defines the model's equations
%    params : model's parameters, it is preferable to pass them as a cell array
%             (compulsory with the functional option) but other formats are
%             acceptable
%
% S = RECSSS(MODEL,S,X,OPTIONS) solves the problem with the parameters
% defined by the structure OPTIONS. The fields of the structure are
%    eqsolver         : 'fsolve', 'lmmcp' (default), 'ncpsolve' or 'path'
%    eqsolveroptions  : options structure to be passed to eqsolver
%
% [S,X] = RECSSS(MODEL,S,X,...) returns the value of the response
% variables at steady state.
%
% [S,X,Z] = RECSSS(MODEL,S,X,...) returns the value of the
% expectations variable at steady state.

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
if nargin<4, options = struct([]); end

defaultopt = struct(...
    'eqsolver'        , 'lmmcp',...
    'eqsolveroptions' , struct([]),...
    'functional'      , 0);
warning('off','catstruct:DuplicatesFound')

options = catstruct(defaultopt,options);

if options.functional
  error(['This program cannot solve for the deterministic steady state of a ' ...
         'functional equation problem']);
end
eqsolver        = lower(options.eqsolver);
eqsolveroptions = options.eqsolveroptions;

params = model.params;
e      = model.w'*model.e;
d      = length(s);
m      = length(x);
if isa(model.func,'char')
  func = str2func(model.func);
elseif isa(model.func,'function_handle')
  func   = model.func;
else
  error('model.func must be either a string or a function handle')
end

if norm(numjac(@(S) bounds(func,S,params),s),Inf)>eps
  error('Bounds should be independent of state variables')
end

%% Prepare input variables
X       = [s(:); x(:)];
[LB,UB] = func('b',s,[],[],[],[],[],params);
LB      = [-inf(size(s(:))); LB(:)];
UB      = [+inf(size(s(:))); UB(:)];

%% Solve for the deterministic steady state
[X,~,exitflag] = runeqsolver(@SSResidual,X,LB,UB,eqsolver,eqsolveroptions,...
                             func,params,e,d,m);

if exitflag~=1
  warning('RECS:SSNotFound','Failure to find a deterministic steady state');
end

%% Prepare outputs
s      = X(1:d)';
x      = X(d+1:d+m)';
output = struct('F',1,'Js',0,'Jx',0,'Jsn',0,'Jxn',0,'hmult',0);
z      = func('h',s,x,[],e,s,x,params,output);


function [F,J] = SSResidual(X,func,params,e,d,m)
%% SSRESIDUAL evaluates the equations and Jacobians of the steady-state finding problem

ss     = X(1:d)';
xx     = X(d+1:d+m)';

if nargout==2 % With Jacobian calculation

  output = struct('F',1,'Js',1,'Jx',1,'Jz',1,'Jsn',1,'Jxn',1,'hmult',1);
  if nargout(func)<6
    [zz,hs,hx,hsnext,hxnext]       = func('h',ss,xx,[],e ,ss,xx,params,output);
  else
    [h,hs,hx,hsnext,hxnext,hmult]  = func('h',ss,xx,[],e ,ss,xx,params,output);
    zz     = h.*hmult;
    hs     = hs.*hmult(:,:,ones(d,1));
    hx     = hx.*hmult(:,:,ones(m,1));
    hsnext = hsnext.*hmult(:,:,ones(d,1));
    hxnext = hxnext.*hmult(:,:,ones(m,1));
  end
  [f,fs,fx,fz] = func('f',ss,xx,zz,[],[],[],params,output);
  fz           = permute(fz,[2 3 1]);
  [g,gs,gx]    = func('g',ss,xx,[],e ,[],[],params,output);

  J                  = zeros(d+m,d+m);
  J(1:d,1:d)         = eye(d)-permute(gs,[2 3 1]);
  J(1:d,d+1:d+m)     = -permute(gx,[2 3 1]);
  J(d+1:d+m,1:d)     = permute(fs,[2 3 1])+fz*permute(hs+hsnext,[2 3 1]);
  J(d+1:d+m,d+1:d+m) = permute(fx,[2 3 1])+fz*permute(hx+hxnext,[2 3 1]);
  J                  = sparse(J);

  F      = [ss-g f]';

else % Without Jacobian calculation
  output = struct('F',1,'Js',0,'Jx',0,'Jz',0,'Jsn',0,'Jxn',0,'hmult',0);
  zz     = func('h',ss,xx,[],e ,ss,xx,params,output);
  g      = func('g',ss,xx,[],e ,[],[],params,output);
  f      = func('f',ss,xx,zz,[],[],[],params,output);
  F      = [ss-g f]';
end


function B = bounds(func,s0,params)
%% BOUNDS Concatenates lower and upper bounds to permit differentiation

[LB,UB]     = func('b',s0,[],[],[],[],[],params);
B           = [LB(:); UB(:)];
B(isinf(B)) = 0;

return
