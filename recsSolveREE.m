function [c,x,z,f,exitflag] = recsSolveREE(interp,model,s,x,options)
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

if nargin <=4, options = struct([]); end

defaultopt = struct(...
    'display'          , 1          ,...
    'eqsolver'         , 'ncpsolve' ,...
    'eqsolveroptions'  , struct([]) ,...
    'loop_over_s'      , 0          ,...
    'method'           , 'expapprox',...
    'reesolver'        , 'krylov'   ,...
    'reesolveroptions' , struct([]) ,...
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
                                    'rtol'      , eps),...
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
if options.functional, params = [params fspace c]; end

[n,m] = size(x);
p     = size(func('h',s(1,:),x(1,:),[],e(1,:),s(1,:),x(1,:),params),2);
K     = length(w);               % number of shock values

switch reesolver
 case 'mixed'
  reesolveroptions.maxit = 10;
  reesolveroptions.atol  = 1E-2;
  reesolveroptions.rtol  = 1E-3;
  c = SA(@Residual_Function, c(:), reesolveroptions);

  reesolveroptions.maxit = 40;
  reesolveroptions.atol  = 1E-7;
  reesolveroptions.rtol  = 1E-25;
  [c,~,exitflag] = nsoli(@Residual_Function, c(:), reesolveroptions);
  if exitflag==0, exitflag = 1; else exitflag = 0; end

 case 'krylov'
  [c,~,exitflag] = nsoli(@Residual_Function, c(:), reesolveroptions);
  if exitflag==0, exitflag = 1; else exitflag = 0; end

 case 'sa'
  [c,~,exitflag] = SA(@Residual_Function, c(:), reesolveroptions);

 case 'fsolve' % In test - Slow, because it uses numerical derivatives
  if options.display==1
    reesolveroptions = optimset('display','iter-detailed','Diagnostics','on');
  end
  [c,~,exitflag] = fsolve(@Residual_Function, c(:), reesolveroptions);

 case 'kinsol'
  neq = numel(c);
  KINoptions  = KINSetOptions('Verbose',       false,...
                              'LinearSolver',  'GMRES',...
                              'ErrorMessages', false,...
                              'FuncNormTol',   reesolveroptions.atol);
  KINInit(@Residual_Function,neq,KINoptions);
  [status, c] = KINSol(c(:),'LineSearch',ones(neq,1),ones(neq,1));
  KINFree;
  if status==0 || status==1, exitflag = 1; else exitflag = 0; end

end

if exitflag~=1
  warning('recs:FailureREE','Failure to find a rational expectations equilibrium');
end

c = reshape(c,n,[]);

if isempty(z)
  ind   = (1:n);
  ind   = ind(ones(1,K),:);
  ss    = s(ind,:);
  xx    = x(ind,:);
  snext = func('g',ss,xx,[],e(repmat(1:K,1,n),:),[],[],params);
  switch method
   case 'expfunapprox'
    h   = funeval(c,fspace,snext);
    if nargout(func)==5
      [~,~,~,~,hmult] = func('h',[],[],[],e(repmat(1:K,1,n),:),snext,zeros(size(snext,1),m),params);
      h               = h.*hmult;
    end

   case 'resapprox-complete'
    [LB,UB] = func('b',snext,[],[],[],[],[],params);
    xnext   = min(max(funeval(c,fspace,snext),LB),UB);
    if nargout(func)<5
       h     = func('h',ss,xx,[],e(repmat(1:K,1,n),:),snext,xnext,params);
    else
      [h,~,~,~,hmult] = func('h',ss,xx,[],e(repmat(1:K,1,n),:),snext,xnext,params);
      h               = h.*hmult;
    end

  end
  z     = reshape(w'*reshape(h,K,n*p),n,p);
end

function [R,FLAG] = Residual_Function(cc)
% RESIDUAL_FUNCTION Calculates the residual of the model with regards to rational
% expectations

  switch method
    case 'expapprox'
     cc    = reshape(cc,n,p);

     z     = funeval(cc,fspace,Phi);
     [x,f] = recsSolveEquilibrium(s,x,z,func,params,cc,e,w,fspace,options);

     ind   = (1:n);
     ind   = ind(ones(1,K),:);
     ss    = s(ind,:);
     xx    = x(ind,:);
     snext = func('g',ss,xx,[],e(repmat(1:K,1,n),:),[],[],params);

     % xnext calculated by interpolation
     xnext = funeval(funfitxy(fspace,Phi,x),fspace,snext);
%{
     % xnext calculated by equation solve
     xnext = recsSolveEquilibrium(snext,...
                                  funeval(funfitxy(fspace,Phi,x),fspace,snext),...
                                  funeval(c,fspace,snext),...
                                  func,params,cc,e,w,fspace,options);
%}

     if nargout(func)<5
       h     = func('h',ss,xx,[],e(repmat(1:K,1,n),:),snext,xnext,params);
     else
       [h,~,~,~,hmult] = func('h',ss,xx,[],e(repmat(1:K,1,n),:),snext,xnext,params);
       h               = h.*hmult;
     end
     z     = reshape(w'*reshape(h,K,n*p),n,p);

     R     = funfitxy(fspace,Phi,z)-cc;

   case 'expfunapprox'
    cc    = reshape(cc,n,p);

    [x,f] = recsSolveEquilibrium(s,x,zeros(n,0),func,params,cc,e,w, fspace,options);
    R     = funfitxy(fspace,Phi,func('h',[],[],[],[],s,x,params))-cc;
    z     = [];

   case 'resapprox-simple'
    cc    = reshape(cc,n,m);

    [LB,UB] = func('b',s,[],[],[],[],[],params);
    x     = min(max(funeval(cc,fspace,Phi),LB),UB);

    ind   = (1:n);
    ind   = ind(ones(1,K),:);
    ss    = s(ind,:);
    xx    = x(ind,:);
    snext = func('g',ss,xx,[],e(repmat(1:K,1,n),:),[],[],params);

    [LB,UB] = func('b',snext,[],[],[],[],[],params);
    xnext = min(max(funeval(cc,fspace,snext),LB),UB);

    if nargout(func)<5
       h     = func('h',ss,xx,[],e(repmat(1:K,1,n),:),snext,xnext,params);
    else
      [h,~,~,~,hmult] = func('h',ss,xx,[],e(repmat(1:K,1,n),:),snext,xnext,params);
      h               = h.*hmult;
    end
    z     = reshape(w'*reshape(h,K,n*p),n,p);

    [x,f] = recsSolveEquilibrium(s,x,z,func,params,cc,e,w,fspace,options);
    R     = funfitxy(fspace,Phi,x)-cc;

   case 'resapprox-complete'
    cc    = reshape(cc,n,m);
    [x,f] = recsSolveEquilibrium(s,x,zeros(n,0),func,params,cc,e,w,fspace,options);
    R     = funfitxy(fspace,Phi,x)-cc;
    z     = [];
  end

  R       = R(:);
  FLAG    = 0;
end

end
