function [c,x,z,f,exitflag] = recsSolveREE(interp,model,s,x,options)
% RECSSOLVEREE finds the rational expectations equilibrium (REE) of a model
%
% RECSSOLVEREE implementes various approximation schemes, residual definition
% (collocation/least-squares), and equation solvers to find the REE of a model.
%
% C = RECSSOLVEREE(INTERP,MODEL,S,X) tries to find the rational expectations
% equilibrium of the model defined in the structure MODEL, by using the
% interpolation structure defined in the structure INTERP. The problem is solved
% on the grid of state variables provided in matrix S. Matrix X is used as a
% first guess of response variables on the grid. RECSSOLVEREE returns the
% coefficient of interpolation corresponding to the functions approximated.
% INTERP is a structure, which has to include the following fields:
%    ch, cx or cz : a coefficient matrix providing a first guess of the
%                   approximation of the expectations function for ch, of the
%                   response variables for cx, or of the expectations for cz
%    fspace       : a definition structure for the interpolation family (created
%                   by the function fundef)
%    Phi          : a basis structure defined on the grid S (created by funbas
%                   or funbasx)
% MODEL is a structure, which has to include the following fields:
%    [e,w] : discrete distribution with finite support with e the values and w the
%            probabilities (it could be also the discretisation of a continuous
%            distribution through quadrature or Monte Carlo drawings)
%    func   : function name or anonymous function that defines the model's equations
%    params : model's parameters, it is preferable to pass them as a cell array
%             (compulsory with the functional option) but other formats are
%             acceptable
%
% C = RECSSOLVEREE(INTERP,MODEL,S,X,OPTIONS) solves the problem with the parameters
% defined by the structure OPTIONS. The fields of the structure are
%    display          : 1 to show iterations (default: 1)
%    eqsolver         : 'fsolve', 'lmmcp', 'ncpsolve' (default) or 'path'
%    eqsolveroptions  : options structure to be passed to eqsolver
%    functional       : 1 if the equilibrium equations are a functional equation
%                       problem (default: 0)
%    loop_over_s      : 0 (default) to solve all grid points at once or 1 to loop
%                       over each grid points
%    method           : 'expapprox' (default), 'expfunapprox', 'resapprox-simple'
%                       or 'resapprox-complete'
%    reesolver        : 'krylov' (default), 'mixed', 'SA' or 'fsolve' (in test)
%    reesolveroptions : options structure to be passed to reesolver
%    useapprox        : (default: 1) behaviour dependent of the chosen method. If 0 and
%                       method is 'expapprox' then next-period responses are
%                       calculated by equations solve and not just
%                       interpolated. If 1 and method is 'resapprox', the guess
%                       of response variables is found with the new
%                       approximation structure
%
% [C,X] = RECSSOLVEREE(INTERP,MODEL,S,X,...) returns the value of the response
% variables on the grid.
%
% [C,X,Z] = RECSSOLVEREE(INTERP,MODEL,S,X,...) returns the value of the
% expectations variables on the grid.
%
% [C,X,Z,F] = RECSSOLVEREE(INTERP,MODEL,S,X,...) returns the value of the
% equilibrium equations on the grid.
%
% [C,X,Z,F,EXITFLAG] = RECSSOLVEREE(INTERP,MODEL,S,X,...) returns EXITFLAG,
% which describes the exit conditions. Possible values are
%    1 : RECSSOLVEREE converges to the REE
%    0 : Failure to converge
%
% See also FUNBAS, FUNBASX, FUNDEF, RECSCHECK, RECSSIMUL, RECSSS.

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

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

switch reesolver
 % Attention: the variables x and z are changed by the nested function
 % 'Residual_Function'
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

c = reshape(c,nc,[]);

% Calculation of z on the grid for output
if isempty(z)
  if functional, params{end} = c; end
  ind    = (1:ns);
  ind    = ind(ones(1,K),:);
  ss     = s(ind,:);
  xx     = x(ind,:);
  output = struct('F',1,'Js',0,'Jx',0);
  snext  = func('g',ss,xx,[],e(repmat(1:K,1,ns),:),[],[],params,output);
  switch method
   case 'expfunapprox'
    h   = funeval(c,fspace,snext);
    if nargout(func)==6
      output            = struct('F',0,'Js',0,'Jx',0,'Jsn',0,'Jxn',0,'hmult',1);
      [~,~,~,~,~,hmult] = func('h',[],[],[],e(repmat(1:K,1,ns),:),snext,zeros(size(snext,1),m),params,output);
      h                 = h.*hmult;
    end

   case 'resapprox-complete'
    [LB,UB] = func('b',snext,[],[],[],[],[],params);
    xnext   = min(max(funeval(c,fspace,snext),LB),UB);
    output  = struct('F',1,'Js',0,'Jx',0,'Jsn',0,'Jxn',0,'hmult',1);
    if nargout(func)<6
       h                = func('h',ss,xx,[],e(repmat(1:K,1,ns),:),snext,xnext,params,output);
    else
      [h,~,~,~,~,hmult] = func('h',ss,xx,[],e(repmat(1:K,1,ns),:),snext,xnext,params,output);
      h               = h.*hmult;
    end

  end
  z     = reshape(w'*reshape(h,K,ns*p),ns,p);
end

function [R,FLAG] = Residual_Function(cc)
% RESIDUAL_FUNCTION Calculates the residual of the model with regards to rational expectations

  switch method
    case 'expapprox'
     cc    = reshape(cc,nc,p);
     if functional, params{end} = cc; end

     z     = funeval(cc,fspace,Phi);
     [x,f] = recsSolveEquilibrium(s,x,z,func,params,cc,e,w,fspace,options);

     ind     = (1:ns);
     ind     = ind(ones(1,K),:);
     ss      = s(ind,:);
     xx      = x(ind,:);
     output  = struct('F',1,'Js',0,'Jx',0);
     snext   = func('g',ss,xx,[],e(repmat(1:K,1,ns),:),[],[],params,output);
     [LB,UB] = func('b',snext,[],[],[],[],[],params);

     % xnext calculated by interpolation
     xnext   = min(max(funeval(funfitxy(fspace,Phi,x),fspace,snext),LB),UB);
     if ~useapprox  % xnext calculated by equation solve
       xnext = recsSolveEquilibrium(snext,...
                                    xnext,...
                                    funeval(cc,fspace,snext),...
                                    func,params,cc,e,w,fspace,options);
     end

     output              = struct('F',1,'Js',0,'Jx',0,'Jsn',0,'Jxn',0,'hmult',1);
     if nargout(func)<6
       h                 = func('h',ss,xx,[],e(repmat(1:K,1,ns),:),snext,xnext,params,output);
     else
       [h,~,~,~,~,hmult] = func('h',ss,xx,[],e(repmat(1:K,1,ns),:),snext,xnext,params,output);
       h                 = h.*hmult;
     end
     z     = reshape(w'*reshape(h,K,ns*p),ns,p);

     R     = funfitxy(fspace,Phi,z)-cc;

   case 'expfunapprox'
    cc    = reshape(cc,nc,p);
    if functional, params{end} = cc; end

    [x,f]  = recsSolveEquilibrium(s,x,z,func,params,cc,e,w, fspace,options);
    output = struct('F',1,'Js',0,'Jx',0,'Jsn',0,'Jxn',0,'hmult',0);
    R      = funfitxy(fspace,Phi,func('h',[],[],[],[],s,x,params,output))-cc;

   case {'resapprox-complete','resapprox-simple'}
    cc    = reshape(cc,nc,m);
    if functional, params{end} = cc; end

    if useapprox % x calculated by interpolation
      [LB,UB] = func('b',s,[],[],[],[],[],params);
      x       = min(max(funeval(cc,fspace,Phi),LB),UB);
    end % if not previous x is used

    if strcmp(method,'resapprox-simple')
      ind    = (1:ns);
      ind    = ind(ones(1,K),:);
      ss     = s(ind,:);
      xx     = x(ind,:);
      output = struct('F',1,'Js',0,'Jx',0);
      snext  = func('g',ss,xx,[],e(repmat(1:K,1,ns),:),[],[],params,output);

      [LB,UB] = func('b',snext,[],[],[],[],[],params);
      xnext   = min(max(funeval(cc,fspace,snext),LB),UB);

      output              = struct('F',1,'Js',0,'Jx',0,'Jsn',0,'Jxn',0,'hmult',1);
      if nargout(func)<6
        h                 = func('h',ss,xx,[],e(repmat(1:K,1,ns),:),snext,xnext,params,output);
      else
        [h,~,~,~,~,hmult] = func('h',ss,xx,[],e(repmat(1:K,1,ns),:),snext,xnext,params,output);
        h                 = h.*hmult;
      end
      z     = reshape(w'*reshape(h,K,ns*p),ns,p);
    end

    [x,f] = recsSolveEquilibrium(s,x,z,func,params,cc,e,w,fspace,options);
    R     = funfitxy(fspace,Phi,x)-cc;
  end

  R       = R(:);
  FLAG    = 0;
end

end
