function [interp,x,z,fval,exitflag,output] = recsSolveREE(interp,model,s,x,options)
% RECSSOLVEREE finds the rational expectations equilibrium (REE) of a model
%
% RECSSOLVEREE implementes various approximation schemes, and equation solvers to
% find the REE of a model.
%
% INTERP = RECSSOLVEREE(INTERP,MODEL,S,X) tries to find the rational expectations
% equilibrium of the model defined in the object MODEL, by using the
% interpolation structure defined in the structure INTERP. The problem is solved
% on the grid of state variables provided in matrix S. Matrix X is used as a
% first guess of response variables on the grid. RECSSOLVEREE returns the
% interpolation structure containing the coefficient matrices cx
% and cz, and ch if this field was initially included in INTERP.
% INTERP is a structure, which has to include the following field:
%    fspace       : a definition structure for the interpolation family (created
%                   by the function fundef)
% Optionally INTERP can also include first guess for the coefficients of
% approximation. If absent, an approximation is made from X.
%    ch, cx or cz : a coefficient matrix providing a first guess of the
%                   approximation of the expectations function for ch, of the
%                   response variables for cx, or of the expectations for cz
% MODEL is an object created by recsmodel.
%
% INTERP = RECSSOLVEREE(INTERP,MODEL,S,X,OPTIONS) solves the problem with the
% parameters defined by the structure OPTIONS. The fields of the structure are
%    display          : 1 to show iterations (default: 1)
%    eqsolver         : 'fsolve', 'lmmcp' (default), 'ncpsolve' or 'path'
%    eqsolveroptions  : options structure to be passed to eqsolver
%    extrapolate      : 1 or 2 if extrapolation is allowed outside the
%                       interpolation space, 0 or -1 to forbid it (default: 1).
%                       For -1 and 2, RECSSOLVEREE displays a warning if state
%                       variables exceed the interpolation space.
%    funapprox        : 'expapprox', 'expfunapprox', or 'resapprox' (default)
%    functional       : 1 if the equilibrium equations are a functional equation
%                       problem (default: 0)
%    loop_over_s      : 0 (default) to solve all grid points at once, 1 to loop
%                       over each grid points, or n to loop over n blocks of
%                       grid points
%    reemethod        : 'iter' (default) or '1-step'
%    reesolver        : 'krylov', 'mixed', 'SA' (default) or 'fsolve'
%    reesolveroptions : options structure to be passed to reesolver
%    useapprox        : (default: 1) behaviour dependent of the chosen function to
%                       approximate. If 0 and funapprox is 'expapprox' then
%                       next-period responses are calculated by equations solve and
%                       not just interpolated. If 1 and funapprox is 'resapprox', the
%                       guess of response variables is found with the new
%                       approximation structure
%  UseParallel        : 'always' (default) to use parallel calculation (require
%                       Parallel Computing Toolbox)' or never'
%
% [INTERP,X] = RECSSOLVEREE(INTERP,MODEL,S,X,...) returns the value of the response
% variables on the grid.
%
% [INTERP,X,Z] = RECSSOLVEREE(INTERP,MODEL,S,X,...) returns the value of the
% expectations variables on the grid.
%
% [INTERP,X,Z,FVAL] = RECSSOLVEREE(INTERP,MODEL,S,X,...) returns the value of the
% equilibrium equations on the grid.
%
% [INTERP,X,Z,FVAL,EXITFLAG] = RECSSOLVEREE(INTERP,MODEL,S,X,...) returns EXITFLAG,
% which describes the exit conditions. Possible values are
%    1 : RECSSOLVEREE converges to the REE
%    0 : Failure to converge
%
% [INTERP,X,Z,FVAL,EXITFLAG,OUTPUT] = RECSSOLVEREE(INTERP,MODEL,S,X,...) returns
% OUTPUT, a structure containing the fields snextmin and snextmax, minimum and
% maximum of next-period state variables.
%
% See also RECSCHECK, RECSSIMUL, RECSSS.

% Copyright (C) 2011-2016 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
defaultopt = struct(                                        ...
    'ArrayProblem'      , false                            ,...
    'display'           , 1                                ,...
    'eqsolver'          , 'lmmcp'                          ,...
    'eqsolveroptions'   , struct('Diagnostics'    , 'off' ,...
                                 'DerivativeCheck', 'off' ,...
                                 'Jacobian'       , 'on')  ,...
    'explicit'          , 0                                ,...
    'extrapolate'       , 1                                ,...
    'functional'        , 0                                ,...
    'loop_over_s'       , 0                                ,...
    'funapprox'         , 'resapprox'                      ,...
    'reemethod'         , 'iter'                           ,...
    'reesolver'         , 'sa'                             ,...
    'reesolveroptions'  , struct([])                       ,...
    'useapprox'         , 1                                ,...
    'UseParallel'       , 'always');
if nargin <=4
  options = defaultopt;
else
  if isfield(options,'eqsolveroptions')
    options.eqsolveroptions = catstruct(defaultopt.eqsolveroptions,options.eqsolveroptions);
  end
  options = catstruct(defaultopt,options);
end

extrapolate        = options.extrapolate;
funapprox          = lower(options.funapprox);
functional         = options.functional;
reemethod          = lower(options.reemethod);

% Extract fields of model
b      = model.functions.b;
e      = model.shocks.e;
g      = model.functions.g;
h      = model.functions.h;
params = model.params;
w      = model.shocks.w;

if strcmp(reemethod,'1-step') && strcmp(funapprox,'expapprox')
  warning('RECS:Switching2Iterative',...
          ['Solving the rational expectations problem is not implemented when ' ...
           'approximating this function. Switching to the default options.'])
  funapprox         = 'resapprox';
  options.funapprox = funapprox;
  reemethod         = 'iter';
end

% Get s from interp structure
if nargin<=2 || isempty(s), s = interp.s; end

% Identify variables dimensions
n       = size(s,1);
[d,m,p] = model.dim{1:3};
k       = length(w);               % number of shock values

% Get x or generate it by interpolation if missing
if nargin<=3 || isempty(x)
  if isfield(interp,'x'), x = interp.x;
  elseif isfield(interp,'cx')
    [LB,UB] = b(s,params);
    x       = min(max(funeval(interp.cx,interp.fspace,s),LB),UB);
  else
    error(['Solving the rational expectations requires a first guess for ' ...
           'response variables'])
  end
end
validateattributes(x,{'numeric'},{'size',[n,m],'nonempty'},4)

% Precalculations
z      = zeros(n,0);
ind    = (1:n);
ind    = ind(ones(1,k),:);
ss     = s(ind,:);
ee     = e(repmat(1:k,1,n),:);

% Extract fields of interp
fspace     = interp.fspace;
Phi        = interp.Phi;
% If the coefficients of approximation are not present in interp, they are
% calculated from the first guess on x
switch funapprox
  case 'expapprox'
    if isfield(interp,'cz') && ~isempty(interp.cz)
      c      = interp.cz;
    elseif functional
      error(['With functional problems, a first guess has to be provided for ' ...
             'the approximation coefficients.'])
    else
      xx      = x(ind,:);
      snext   = g(ss,xx,ee,params);
      if extrapolate>=1, snextinterp = snext;
      else
        snextinterp = max(min(snext,fspace.b(ones(n*k,1),:)),fspace.a(ones(n*k,1),:));
      end
      [LBnext,UBnext] = b(snext,params);
      xnext   = min(max(funeval(funfitxy(fspace,Phi,x),fspace,snextinterp),LBnext),UBnext);
      hv      = h(ss,xx,ee,snext,xnext,params);
      z       = reshape(w'*reshape(hv,k,n*p),n,p);
      c       = funfitxy(fspace,Phi,z);
    end
  case 'expfunapprox'
    if isfield(interp,'ch') && ~isempty(interp.ch)
      c      = interp.ch;
    elseif functional
      error(['With functional problems, a first guess has to be provided for ' ...
             'the approximation coefficients.'])
    else
      c       = funfitxy(fspace,Phi,h(zeros(n,0),[],[],s,x,params));
    end
  otherwise
    if isfield(interp,'cx') && ~isempty(interp.cx)
      c = interp.cx;
    else
      c = funfitxy(fspace,Phi,x);
    end
end
if functional
  model.params = [model.params fspace c];
  params       = model.params;
end

%% Solve for the rational expectations equilibrium
switch reemethod
  case 'iter'
    [c,x,z,fval,exitflag] = recsSolveREEIter(interp,model,s,x,c,options);
  case '1-step'
    [c,x,fval,exitflag] = recsSolveREEFull(interp,model,s,x,c,options);
  otherwise
    error(['%s is not a valid value for reemethod. Valid values are ''iter'' ' ...
           '(default) and ''1-step''.'],reemethod)
end

if exitflag~=1
  warning('RECS:FailureREE','Failure to find a rational expectations equilibrium');
end

%% Outputs calculations

% Interpolation coefficients
c = reshape(c,n,[]);
if functional, params{end} = c; end
switch funapprox
 case 'expapprox'
  interp.cz = c;
  interp.cx = funfitxy(fspace,Phi,x);
  if strcmp(model.infos.model_type,'fgh1')
    interp.ch = funfitxy(fspace,Phi,h(zeros(n,0),[],[],s,x,params));
  end
 case 'expfunapprox'
  interp.ch = c;
  interp.cx = funfitxy(fspace,Phi,x);
 otherwise
  interp.cx = c;
  if ~isempty(z), interp.cz = funfitxy(fspace,Phi,z); end
  if strcmp(model.infos.model_type,'fgh1')
    interp.ch = funfitxy(fspace,Phi,h(zeros(n,0),[],[],s,x,params));
  end
end

% Check if state satisfies bounds
xx      = x(ind,:);
snext   = g(ss,xx,ee,params);
output  = struct('snextmin',min(snext),'snextmax',max(snext));
vari    = 1:d;
varmin  = vari(min(snext)<fspace.a);
varmax  = vari(max(snext)>fspace.b);
if extrapolate==2 || extrapolate==-1
  if ~isempty(varmin)
    warning('RECS:Extrapolation','State variables (%s) beyond smin',...
            int2str(varmin))
  end
  if ~isempty(varmax)
    warning('RECS:Extrapolation','State variables (%s) beyond smax',...
            int2str(varmax))
  end
end

% Calculation of z on the grid for output
if isempty(z)
  if extrapolate>=1, snextinterp = snext;
  else
    snextinterp = max(min(snext,fspace.b(ones(n*k,1),:)),fspace.a(ones(n*k,1),:));
  end

  switch funapprox
   case 'expfunapprox'
    hv   = funeval(c,fspace,snextinterp);

   case 'resapprox'
    [LBnext,UBnext] = b(snext,params);
    xnext           = min(max(funeval(c,fspace,snextinterp),LBnext),UBnext);
    hv              = h(ss,xx,ee,snext,xnext,params);

  end
  z         = reshape(w'*reshape(hv,k,n*p),n,p);
  interp.cz = funfitxy(fspace,Phi,z);
end

interp.x = x;
interp.z = z;
