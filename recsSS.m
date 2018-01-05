function [s,x,z,exitflag] = recsSS(model,s,x,options)
% RECSSS Solves for the deterministic steady state in rational expectations models
%
% RECSSS cannot find the deterministic steady state of a functional
% equation problem since, in this case, the steady state depends on
% the model solution.
%
% S = RECSSS(MODEL) tries to find the non-stochastic steady state of the model
% defined in the object MODEL. This function call uses as first guess for
% steady-state state and response variable the values available in the
% properties sss and xss of the object MODEL. RECSSS returns the value of the
% state variables at steady state. MODEL is an object created by recsmodel.
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

% Copyright (C) 2011-2018 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
if nargin<2 || isempty(s), s = model.sss; end
if nargin<3 || isempty(x), [~,x] = model.xss; end

defaultopt = struct(...
    'display'         , 1                               ,...
    'eqsolver'        , 'lmmcp'                         ,...
    'eqsolveroptions' , struct('Diagnostics'    , 'off',...
                               'DerivativeCheck', 'off',...
                               'Jacobian'       , 'on') ,...
    'functional'      , 0);
if nargin<4
  options = defaultopt;
else
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
e      = model.shocks.w'*model.shocks.e;
[d,m]  = model.dim{1:2};
fp     = model.functions.fp;
gp     = model.functions.gp;
hp     = model.functions.hp;

%% Solve for the deterministic steady state
nx = model.infos.nxvarbounds;
w  = zeros(nx(1),1);
v  = zeros(nx(2),1);
X  = [s(:); x(:); w; v];

[LBx,UBx] = model.functions.bp(s,params);
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
z = model.functions.h(s,x,e,s,x,params);

%% Display steady state
if exitflag==1 && options.display==1
  deltass = max(abs([s x]-[s0 x0]));
  if deltass<sqrt(eps)
    fprintf(1,'Deterministic steady state (equal to first guess)\n')
  else
    fprintf(1,['Deterministic steady state (different from first guess, ' ...
               'max(|delta|)=%g)\n'],deltass)
  end
  if exist('table','file')
    fprintf(1,' State variables:\n')
    disp(array2table(s,'VariableNames',model.symbols.states))
    fprintf(1,' Response variables:\n')
    disp(array2table(x,'VariableNames',model.symbols.controls))
    fprintf(1,' Expectations variables:\n')
    disp(array2table(z,'VariableNames',model.symbols.expectations))
  else
    fprintf(1,' State variables:\n\t\t')
    fprintf(1,'%0.4g\t',s)
    fprintf(1,'\n\n Response variables:\n\t\t')
    fprintf(1,'%0.4g\t',x)
    fprintf(1,'\n\n Expectations variables:\n\t\t')
    fprintf(1,'%0.4g\t',z)
    fprintf(1,'\n\n')
  end
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
  [zz,hs,hx,~,hsnext,hxnext] = hp(ss,xx,e,ss,xx,params,[1 1 1 0 1 1]);
  [f,fs,fx,fz]               = fp(ss,xx,ww,vv,zz,params,ones(4,1));
  [g,gs,gx]                  = gp(ss,xx,e,params,[1 1 1 0]);
  fz                         = permute(fz,[2 3 1]);

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
  zz = hp(ss,xx,e,ss,xx,params,[1 0 0 0 0 0]);
  f  = fp(ss,xx,ww,vv,zz,params,[1 0 0 0]);
  g  = gp(ss,xx,e,params,[1 0 0 0]);
end

F = [ss-g f]';

return
