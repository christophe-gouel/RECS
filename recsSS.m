function [s,x,z,exitflag] = recsSS(model,s,x,options)
% RECSSS Solves for the deterministic steady state in rational expectations models
%
% RECSSS cannot find the deterministic steady state of a functional
% equation problem since, in this case, the steady state depends on
% the model solution.
%
% S = RECSSS(MODEL) tries to find the non-stochastic steady state of the model
% defined in the structure MODEL. This function call uses as first guess for
% steady-state state and response variable the output of the call
% model.func('ss'). RECSSS returns the value of the state variables at steady state.
% MODEL is a structure, which has to include the following fields:
%    [e,w] : discrete distribution with finite support with e the values and w the
%            probabilities (it could be also the discretisation of a continuous
%            distribution through quadrature or Monte Carlo drawings)
%    func   : function name or anonymous function that defines the model's equations
%    params : model's parameters, it is preferable to pass them as a cell array
%             (compulsory with the functional option) but other formats are
%             acceptable
%
% S = RECSSS(MODEL,S) uses the vector S as first guess for steady-state state
% variables.
%
% S = RECSSS(MODEL,S,X) uses the vector X as first guess for steady-state response
% variables.
%
% S = RECSSS(MODEL,S,X,OPTIONS) solves the problem with the parameters
% defined by the structure OPTIONS. The fields of the structure are
%    display          : 1 to display the steady state if found (default: 1)
%    eqsolver         : 'fsolve', 'lmmcp' (default), 'ncpsolve' or 'path'
%    eqsolveroptions  : options structure to be passed to eqsolver
%
% [S,X] = RECSSS(MODEL,...) returns the value of the response
% variables at steady state.
%
% [S,X,Z] = RECSSS(MODEL,...) returns the value of the
% expectations variable at steady state.
%
% [S,X,Z,EXITFLAG] = RECSSS(MODEL,...) returns EXITFLAG,
% which describes the exit conditions. Possible values are
%    1 : RECSSS converges to the deterministic steady state
%    0 : Failure to converge

% Copyright (C) 2011-2013 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
if nargin<2 || isempty(s), s = model.func('ss'); end
if nargin<3 || isempty(x), [~,x] = model.func('ss'); end

defaultopt = struct(...
    'display'         , 1                               ,...
    'eqsolver'        , 'lmmcp'                         ,...
    'eqsolveroptions' , struct('DerivativeCheck', 'off',...
                               'Jacobian'       , 'on') ,...
    'functional'      , 0);
if nargin<4
  options = defaultopt;
else
  warning('off','catstruct:DuplicatesFound')
  if isfield(options,'eqsolveroptions')
    options.eqsolveroptions = catstruct(defaultopt.eqsolveroptions,options.eqsolveroptions);
  end
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

%% Solve for the deterministic steady state

[~,~,dLBxds,dUBxds] = model.b(s,params);
dLBxds = permute(dLBxds,[3 2 1]);
dUBxds = permute(dUBxds,[3 2 1]);
ix = [sum(dLBxds~=0,1,'native')' sum(dUBxds~=0,1,'native')'];
nx = int16(sum(ix,1));
w = zeros(nx(1),1);
v = zeros(nx(2),1);
X = [s(:); x(:); w; v];

model = mcptransform(model,ix,nx);
fp    = model.fp;
gp    = model.gp;
hp    = model.hp;

[LBx,UBx] = model.bp(s,params);
LB = [-inf(size(s(:))); LBx(:)];
UB = [+inf(size(s(:))); UBx(:)];

[X,~,exitflag] = runeqsolver(@SSResidual,X,LB,UB,eqsolver,eqsolveroptions,...
                             fp,gp,hp,params,e,d,m,nx);

if exitflag~=1
  warning('RECS:SSNotFound','Failure to find a deterministic steady state');
end

%% Prepare outputs
s0     = s;
x0     = x;
s      = X(1:d)';
x      = X(d+1:d+m)';
output = struct('F',1,'Js',0,'Jx',0,'Jsn',0,'Jxn',0,'hmult',1);
if nargout(model.h)<6
  z = model.h(s,x,e,s,x,params,output);
else
  [h,~,~,~,~,hmult] = model.h(s,x,e,s,x,params,output);
  z = h.*hmult;
end

%% Display steady state
if exitflag==1 && options.display==1
  deltass = max(abs([s x]-[s0 x0]));
  if deltass<sqrt(eps)
    fprintf(1,'Deterministic steady state (equal to first guess)\n')
  else
    fprintf(1,['Deterministic steady state (different from first guess, ' ...
               'max(|delta|)=%g)\n'],deltass)
  end
  fprintf(1,' State variables:\n\t\t')
  fprintf(1,'%0.4g\t',s)
  fprintf(1,'\n\n Response variables:\n\t\t')
  fprintf(1,'%0.4g\t',x)
  fprintf(1,'\n\n Expectations variables:\n\t\t')
  fprintf(1,'%0.4g\t',z)
  fprintf(1,'\n\n')
end


function [F,J] = SSResidual(X,fp,gp,hp,params,e,d,m,nx)
%% SSRESIDUAL evaluates the equations and Jacobians of the steady-state finding problem

ss     = X(1:d)';
xx     = X(d+1:d+m)';
ww     = X(d+m+1:d+m+nx(1))';
vv     = X(d+m+nx(1)+1:d+m+nx(1)+nx(2))';

M = m+nx(1)+nx(2);

if nargout==2
  %% With Jacobian calculation
  output = struct('F',1,'Js',1,'Jx',1,'Jz',1,'Jsn',1,'Jxn',1);
  [zz,hs,hx,hsnext,hxnext] = hp(ss,xx,e,ss,xx,params,output);
  [f,fs,fx,fz]             = fp(ss,xx,ww,vv,zz,params,output);
  [g,gs,gx]                = gp(ss,xx,e,params,output);
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
  zz = hp(ss,xx,e,ss,xx,params,output);
  f  = fp(ss,xx,ww,vv,zz,params,output);
  g  = gp(ss,xx,e,params,output);
end

F = [ss-g f]';

return
