function [c,x,z,f,exitflag] = pabSolveREE(interp,model,s,xinit,options)
% RECSSOLVEREE finds the rational expectations equilibrium (REE) of a model for a
% given set of parameters
%
%  OPTIONS structure field:
%    display          : 1 to show iterations (default: 1)
%    eqsolver         : 'fsolve', 'lmmcp', 'ncpsolve' (default) or 'path'
%    loop_over_s      : 0 (default) to solve all grid points at once or 1 to loop
%    over each grid points
%    method           : 'expapprox' (default), 'expfunapprox', 'resapprox-simple'
%                       or 'resapprox-complete'
%    reesolver        : 'krylov' (default), 'mixed', 'SA' or 'fsolve' (in test)
%    reesolveroptions : options structure to be passed to the reesolver

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

x = xinit;

if nargin <=4, options = struct([]); end

defaultopt = struct(...
    'display'          , 1          ,...
    'eqsolver'         , 'ncpsolve' ,...
    'eqsolveroptions'  , struct([]) ,...
    'loop_over_s'      , 0          ,...
    'method'           , 'expapprox',...
    'reesolver'        , 'krylov'   ,...
    'reesolveroptions' , struct([]) ,...
    'residual_function', 1, ...
    'functional'       , 0);
warning('off','catstruct:DuplicatesFound')

options = catstruct(defaultopt,options);

eqsolver         = options.eqsolver;
eqsolveroptions  = options.eqsolveroptions;
loop_over_s      = options.loop_over_s;
method           = lower(options.method);
reesolver        = lower(options.reesolver);

reesolveroptions = catstruct(struct('showiters' , options.display,...
                                    'atol'      , sqrt(eps),...
                                    'lmeth'     , 3,...
                                    'rtol'      , sqrt(eps)),...
                             options.reesolveroptions);

e      = model.e;
params = model.params;
w      = model.w;
if isa(model.func,'char')
  func = str2func(model.func);
elseif isa(model.func,'function_handle')
  func   = model.func;
else
  error('model.func must be either a string or a function handle')
end

c = fit_coefficients(interp, xinit);

[n,m] = size(x);
p     = size(func('h',s(1,:),xinit(1,:),[],e(1,:),s(1,:),x(1,:),params),2);
K     = length(w);               % number of shock values

flatten = @(aa) aa(:);

Residual_Function = @(cc) flatten( pabResidual( reshape(cc, n, m )  ,n,m,p,s,model,interp,options) );

switch reesolver
 case 'mixed'
  reesolveroptions.maxit = 10;
  reesolveroptions.atol  = 1E-3;
  reesolveroptions.rtol  = 1E-4;

  c = SA(Residual_Function, c(:), reesolveroptions);

  reesolveroptions.maxit = 40;
  reesolveroptions.atol  = 1E-7;
  reesolveroptions.rtol  = 1E-25;
  [c,~,exitflag] = nsoli(Residual_Function, c(:), reesolveroptions);

  exitflag = 0;
  if exitflag==0, exitflag = 1; else exitflag = 0; end

 case 'krylov'
  [c,~,exitflag] = nsoli(Residual_Function, c(:), reesolveroptions);
  if exitflag==0, exitflag = 1; else exitflag = 0; end

 case 'sa'
    [c,~,exitflag] = SA(Residual_Function, c(:), reesolveroptions);

 case 'fsolve' % In test - Slow, because it uses numerical derivatives
  reesolveroptions.maxit = 10;
  reesolveroptions.atol  = 1E-2;
  reesolveroptions.rtol  = 1E-3;

  c = SA(Residual_Function, c(:), reesolveroptions);
  if options.display==1
    reesolveroptions = catstruct( optimset('display','iter-detailed','Diagnostics','on'),reesolveroptions)
  end
    disp(reesolveroptions)
  [c,~,exitflag] = fsolve(Residual_Function, c(:), reesolveroptions);

end

if exitflag~=1
  warning('recs:FailureREE','Failure to find a rational expectations equilibrium');
end


c = reshape(c,n,[]);
x = eval_function( interp, c, s );
z = [];
f = [];


return