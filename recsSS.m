function [s,x,z] = recsSS(model,s,x,options)
% RECSSS Solves for the deterministic steady state in rational expectations models

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

if nargin<4, options = struct([]); end

defaultopt = struct(...
    'eqsolver'        , 'lmmcp',...
    'eqsolveroptions' , struct([]));

options = catstruct(defaultopt,options);

eqsolver        = options.eqsolver;
eqsolveroptions = options.eqsolveroptions;

func   = model.func;
params = model.params;
e      = model.w'*model.e;
d      = length(s);
m      = length(x);

exitflag = 1;

X       = [s(:); x(:)];
[LB,UB] = feval(func,'b',s,[],[],[],[],[],params);
LB      = [-inf(size(s(:))); LB(:)];
UB      = [+inf(size(s(:))); UB(:)];

switch lower(eqsolver)
 case 'fsolve'
  opt = optimset('Display','off',...
                 'Jacobian','on');
  opt = optimset(opt,eqsolveroptions);
  [X,~,exitflag] = fsolve(@(Y) recsSSResidual(Y,func,params,e,d,m),...
                          X,...
                          opt);
 case 'lmmcp'
  [X,~,exitflag] = lmmcp(@(Y) recsSSResidual(Y,func,params,e,d,m),...
                         X,...
                         LB,...
                         UB,...
                         eqsolveroptions);
 case 'ncpsolve'
  X = ncpsolve('ncpsolvetransform',...
               LB,...
               UB,...
               X,...
               'recsSSResidual',...
               func,params,e,d,m);
 case 'path'
  global par %#ok<TLEV>
  par = {'recsSSResidual',func,params,e,d,m};
  X   = pathmcp(X,LB,UB,'pathtransform');
  clear global par
end

if exitflag~=1, disp('Failure to find a deterministic steady state'); end

s = X(1:d)';
x = X(d+1:d+m)';
z = feval(func,'h',s,x,[],e,s,x,params);