function [s,x,z,exitflag] = recsSS(model,s,x,options)
% RECSSS Solves for the deterministic steady state in rational expectations models
%
% RECSSS cannot find the deterministic steady state of a functional
% equation problem since, in this case, the steady state depends on
% the model solution.
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
%
% [S,X,Z,EXITFLAG] = RECSSS(MODEL,S,X,...) returns EXITFLAG,
% which describes the exit conditions. Possible values are
%    1 : RECSSS converges to the deterministic steady state
%    0 : Failure to converge

% Copyright (C) 2011-2012 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
if nargin<2 || isempty(s), s = model.func('ss'); end
if nargin<3 || isempty(x), [~,x] = model.func('ss'); end

defaultopt = struct(...
    'eqsolver'        , 'lmmcp',...
    'eqsolveroptions' , struct('DerivativeCheck','off'),...
    'functional'      , 0);
if nargin<4
  options = defaultopt; 
else
  warning('off','catstruct:DuplicatesFound')
  options = catstruct(defaultopt,options);
end

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

%% Solve for the deterministic steady state

ix = [sum(numjac(@(S) Bounds(func,S,params,1,true(m,2)),s)~=0,2,'native') ...
      sum(numjac(@(S) Bounds(func,S,params,2,true(m,2)),s)~=0,2,'native')];
nx = int16(sum(ix,1));

w = zeros(nx(1),1);
v = zeros(nx(2),1);
X = [s(:); x(:); w; v];
[LBx,UBx] = mcptransform(func,'b',s,[],[],[],[],[],params,[],ix,nx);
LB = [-inf(size(s(:))); LBx(:)];
UB = [+inf(size(s(:))); UBx(:)];

[X,~,exitflag] = runeqsolver(@SSResidual,X,LB,UB,eqsolver,...
                             eqsolveroptions,func,params,e,d,m,ix,nx);

if exitflag~=1
  warning('RECS:SSNotFound','Failure to find a deterministic steady state');
end

%% Prepare outputs
s      = X(1:d)';
x      = X(d+1:d+m)';
output = struct('F',1,'Js',0,'Jx',0,'Jsn',0,'Jxn',0,'hmult',1);
if nargout(func)<6
  z = func('h',s,x,[],e,s,x,params,output);
else
  [h,~,~,~,~,hmult] = func('h',s,x,[],e,s,x,params,output);
  z = h.*hmult;
end


function [F,J] = SSResidual(X,func,params,e,d,m,ix,nx)
%% SSRESIDUAL evaluates the equations and Jacobians of the steady-state finding problem

ss     = X(1:d)';
xx     = X(d+1:d+m)';
ww     = X(d+m+1:d+m+nx(1))';
vv     = X(d+m+nx(1)+1:d+m+nx(1)+nx(2))';

M = m+nx(1)+nx(2);

if nargout==2 
  %% With Jacobian calculation
  output = struct('F',1,'Js',1,'Jx',1,'Jz',1,'Jsn',1,'Jxn',1);
  [zz,hs,hx,hsnext,hxnext] = mcptransform(func,'h',ss,xx,[],e ,ss,xx,params,output,ix,nx,[],[]);
  [f,fs,fx,fz]             = mcptransform(func,'f',ss,xx,zz,[],[],[],params,output,ix,nx,ww,vv);
  [g,gs,gx]                = mcptransform(func,'g',ss,xx,[],e ,[],[],params,output,ix,nx,[],[]);
  fz                       = permute(fz,[2 3 1]);

  J               = zeros(d+M,d+M);
  % With respect to s
  J(1:d     ,1:d) = eye(d)-permute(gs,[2 3 1]);
  J(d+1:d+M ,1:d) = permute(fs,[2 3 1])+fz*permute(hs+hsnext,[2 3 1]);
  % With respect to X
  J(1:d     ,d+1:d+M) = -permute(gx,[2 3 1]);
  J(d+1:d+M ,d+1:d+M) =  permute(fx,[2 3 1])+fz*permute(hx+hxnext,[2 3 1]);
  % Aggregation into a sparse matrix
  J = sparse(J);
  
else
  %% Without Jacobian calculation
  output = struct('F',1,'Js',0,'Jx',0,'Jz',0,'Jsn',0,'Jxn',0);
  zz = mcptransform(func,'h',ss,xx,[],e ,ss,xx,params,output,ix,nx,[],[]);
  f  = mcptransform(func,'f',ss,xx,zz,[],[],[],params,output,ix,nx,ww,vv);
  g  = mcptransform(func,'g',ss,xx,[],e ,[],[],params,output,ix,nx,[],[]);
end

F = [ss-g f]';

return
