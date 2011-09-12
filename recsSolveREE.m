function [interp,x,z,f,exitflag] = recsSolveREE(interp,model,s,x,options)
% RECSSOLVEREE finds the rational expectations equilibrium (REE) of a model
%
% RECSSOLVEREE implementes various approximation schemes, and equation solvers to
% find the REE of a model.
%
% INTERP = RECSSOLVEREE(INTERP,MODEL,S,X) tries to find the rational expectations
% equilibrium of the model defined in the structure MODEL, by using the
% interpolation structure defined in the structure INTERP. The problem is solved
% on the grid of state variables provided in matrix S. Matrix X is used as a
% first guess of response variables on the grid. RECSSOLVEREE returns the
% interpolation structure containing the coefficient matrices cx
% and cz, and ch if this field was initially included in INTERP.
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
% INTERP = RECSSOLVEREE(INTERP,MODEL,S,X,OPTIONS) solves the problem with the 
% parameters defined by the structure OPTIONS. The fields of the structure are
%    display          : 1 to show iterations (default: 1)
%    eqsolver         : 'fsolve', 'lmmcp' (default), 'ncpsolve' or 'path'
%    eqsolveroptions  : options structure to be passed to eqsolver
%    extrapolate      : 1 if extrapolation is allowed outside the
%                       interpolation space or 0 to forbid it (default: 1)
%    funapprox        : 'expapprox', 'expfunapprox', 'resapprox-simple'
%                       or 'resapprox-complete' (default)
%    functional       : 1 if the equilibrium equations are a functional equation
%                       problem (default: 0)
%    loop_over_s      : 0 (default) to solve all grid points at once or 1 to loop
%                       over each grid points
%    reemethod        : 'iter' (default) or '1-step'  
%    reesolver        : 'krylov', 'mixed', 'SA' (default) or 'fsolve' (in test)
%    reesolveroptions : options structure to be passed to reesolver
%    useapprox        : (default: 1) behaviour dependent of the chosen function to
%                       approximate. If 0 and funapprox is 'expapprox' then
%                       next-period responses are calculated by equations solve and
%                       not just interpolated. If 1 and funapprox is 'resapprox', the
%                       guess of response variables is found with the new
%                       approximation structure
%
% [INTERP,X] = RECSSOLVEREE(INTERP,MODEL,S,X,...) returns the value of the response
% variables on the grid.
%
% [INTERP,X,Z] = RECSSOLVEREE(INTERP,MODEL,S,X,...) returns the value of the
% expectations variables on the grid.
%
% [INTERP,X,Z,F] = RECSSOLVEREE(INTERP,MODEL,S,X,...) returns the value of the
% equilibrium equations on the grid.
%
% [INTERP,X,Z,F,EXITFLAG] = RECSSOLVEREE(INTERP,MODEL,S,X,...) returns EXITFLAG,
% which describes the exit conditions. Possible values are
%    1 : RECSSOLVEREE converges to the REE
%    0 : Failure to converge
%
% See also RECSCHECK, RECSSIMUL, RECSSS.

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
if nargin <=4, options = struct([]); end

defaultopt = struct(                           ...
    'display'           , 1                   ,...
    'eqsolver'          , 'lmmcp'             ,...
    'eqsolveroptions'   , struct([])          ,...
    'extrapolate'       , 1                   ,...
    'functional'        , 0                   ,...
    'loop_over_s'       , 0                   ,...
    'funapprox'         , 'resapprox-complete',...
    'reemethod'         , 'iter'              ,...
    'reesolver'         , 'sa'                ,...
    'reesolveroptions'  , struct([])          ,...
    'useapprox'         , 1);
warning('off','catstruct:DuplicatesFound')

options = catstruct(defaultopt,options);

extrapolate        = options.extrapolate;
funapprox          = lower(options.funapprox);
functional         = options.functional;
reemethod          = lower(options.reemethod);

switch funapprox
 case 'expapprox'
  c      = interp.cz;
 case 'expfunapprox'
  c      = interp.ch;
 otherwise
  c      = interp.cx;
end
fspace = interp.fspace;
Phi    = interp.Phi;
if functional, model.params = [model.params fspace c]; end

e      = model.e;
params = model.params;
w      = model.w;
if isa(model.func,'char')
  model.func = str2func(model.func);
elseif ~isa(model.func,'function_handle')
  error('model.func must be either a string or a function handle')
end
func = model.func;

if strcmp(reemethod,'1-step') && (strcmp(funapprox,'expapprox') || strcmp(funapprox,'resapprox'))
  warning(['Solving the rational expectations problem is not implemented when ' ...
           'approximating this funtion. Switching to the iterative scheme.'])
end

[n,m]  = size(x);
output = struct('F',1,'Js',0,'Jx',0,'Jsn',0,'Jxn',0,'hmult',0);
p      = size(func('h',s(1,:),x(1,:),[],e(1,:),s(1,:),x(1,:),params,output),2);
k      = length(w);               % number of shock values
z      = zeros(n,0);

%% Solve for the rational expectations equilibrium
switch reemethod
  case 'iter'
    [c,x,z,f,exitflag] = recsSolveREEIter(interp,model,s,x,c,options);
  case '1-step'
    [c,x,f,exitflag] = recsSolveREEFull(interp,model,s,x,c,options);
end

if exitflag~=1
  warning('recs:FailureREE','Failure to find a rational expectations equilibrium');
end

%% Outputs calculations
% Interpolation coefficients
c = reshape(c,n,[]);
switch funapprox
 case 'expapprox'
  interp.cz = c;
  interp.cx = funfitxy(fspace,Phi,x); 
  if isfield(interp,'ch')
    output = struct('F',1,'Js',0,'Jx',0,'Jsn',0,'Jxn',0,'hmult',0);
    interp.ch = funfitxy(fspace,Phi,func('h',[],[],[],[],s,x,params,output));
  end
 case 'expfunapprox'
  interp.ch = c;
  interp.cx = funfitxy(fspace,Phi,x); 
 otherwise
  interp.cx = c;
  if ~isempty(z), interp.cz = funfitxy(fspace,Phi,z); end
  if isfield(interp,'ch')
    output = struct('F',1,'Js',0,'Jx',0,'Jsn',0,'Jxn',0,'hmult',0);
    interp.ch = funfitxy(fspace,Phi,func('h',[],[],[],[],s,x,params,output));
  end
end

% Calculation of z on the grid for output
if isempty(z)
  if functional, params{end} = c; end
  ind    = (1:n);
  ind    = ind(ones(1,k),:);
  ss     = s(ind,:);
  xx     = x(ind,:);
  output = struct('F',1,'Js',0,'Jx',0);
  snext  = func('g',ss,xx,[],e(repmat(1:k,1,n),:),[],[],params,output);
  if extrapolate, snextinterp = snext;
  else      
    snextinterp = max(min(snext,fspace.b(ones(n*k,1),:)),fspace.a(ones(n*k,1),:)); 
  end

  switch funapprox
   case 'expfunapprox'
    h   = funeval(c,fspace,snextinterp);
    if nargout(func)==6
      output            = struct('F',0,'Js',0,'Jx',0,'Jsn',0,'Jxn',0,'hmult',1);
      [~,~,~,~,~,hmult] = func('h',[],[],[],e(repmat(1:k,1,n),:),snext,zeros(size(snext,1),m),params,output);
      h                 = h.*hmult;
    end

   case 'resapprox-complete'
    [LB,UB] = func('b',snextinterp,[],[],[],[],[],params);
    xnext   = min(max(funeval(c,fspace,snextinterp),LB),UB);
    output  = struct('F',1,'Js',0,'Jx',0,'Jsn',0,'Jxn',0,'hmult',1);
    if nargout(func)<6
       h                = func('h',ss,xx,[],e(repmat(1:k,1,n),:),snext,xnext,params,output);
    else
      [h,~,~,~,~,hmult] = func('h',ss,xx,[],e(repmat(1:k,1,n),:),snext,xnext,params,output);
      h               = h.*hmult;
    end

  end
  z         = reshape(w'*reshape(h,k,n*p),n,p);
  interp.cz = funfitxy(fspace,Phi,z); 
end

