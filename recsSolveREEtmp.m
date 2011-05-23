function [c,x,z,f,exitflag] = recsSolveREEtmp(interp,model,s,x,options)

if nargin <=4, options = struct([]); end

defaultopt = struct(                  ...
    'display'           , 1          ,...
    'eqsolver'          , 'ncpsolve' ,...
    'eqsolveroptions'   , struct([]) ,...
    'functional'        , 0          ,...
    'loop_over_s'       , 0          ,...
    'method'            , 'expapprox',...
    'reesolver'         , 'krylov'   ,...
    'reesolveroptions'  , struct([]) ,...
    'useapprox'         , 1);
warning('off','catstruct:DuplicatesFound')

options = catstruct(defaultopt,options);

functional         = options.functional;
method             = lower(options.method);
reesolver          = lower(options.reesolver);
reesolveroptions   = catstruct(struct('showiters' , options.display,...
                                      'atol'      , sqrt(eps),...
                                      'lmeth'     , 3,...
                                      'rtol'      , eps),...
                               options.reesolveroptions);
useapprox          = options.useapprox;

e      = model.e;
params = model.params;
w      = model.w;
if isa(model.func,'char')
  func = str2func(model.func);
elseif isa(model.func,'function_handle')
  func = model.func;
else
  error('model.func must be either a string or a function handle')
end

switch method
 case 'expapprox'
  c      = interp.cz;
 case 'expfunapprox'
  c      = interp.ch;
 otherwise
  c      = interp.cx;
end
fspace = interp.fspace;
Phi    = interp.Phi;
if functional, params = [params fspace c]; end

[ns,m] = size(x);
nc     = size(c,1);
output = struct('F',1,'Js',0,'Jx',0,'Jsn',0,'Jxn',0,'hmult',0);
p      = size(func('h',s(1,:),x(1,:),[],e(1,:),s(1,:),x(1,:),params,output),2);
K      = length(w);               % number of shock values
z      = zeros(ns,0);

[f,J] = recsSolveFullPb(s,x,func,params,c,e,w,fspace,Phi);

%spy(J)
