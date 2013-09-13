function [ssim,xsim,esim,stat,fsim] = recsSimul(model,interp,s0,nper,shocks,options)
% RECSSIMUL Simulates a model from starting values given in s0 and for nper period
%
% SSIM = RECSSIMUL(MODEL,INTERP,S0,NPER) simulates the model defined in the
% object MODEL, by using the interpolation structure defined in the structure
% INTERP. The simulation starts from the initial state S0 and lasts NPER (scalar)
% periods. S0 is a nrep-by-d matrix with nrep the number of scenarios to simulate,
% and d the number of state variables. RECSSIMUL returns the nrep-by-d-by-nper
% array SSIM that contains the simulated state variables.
% MODEL is an object created by recsmodel.
% INTERP is a structure, which includes the following fields:
%    ch      : coefficient matrix of the interpolation of the expectations
%              function (optional, to be provided with method 'expfunapprox')
%    cx      : coefficient matrix of the interpolation of the response variables
%    cz      : coefficient matrix of the interpolation of the expectations variables
%    fspace  : a definition structure for the interpolation family (created by
%              the function fundef)
%
% SSIM = RECSSIMUL(MODEL,INTERP,S0,NPER,SHOCKS) uses the nrep-by-q-by-(nper-1) array
% SHOCKS to simulate the model instead of drawing random numbers. In this case
% size(SHOCKS,3) supersedes NPER.
%
% SSIM = RECSSIMUL(MODEL,INTERP,S0,NPER,SHOCKS,OPTIONS) simulates the model with
% the parameters defined by the structure OPTIONS. The fields of the structure are
%    accuracy         : 1 to check accuracy on the asymptotic distribution (default: 0)
%    display          : 1 (default) to display outputs
%    eqsolver         : 'fsolve', 'lmmcp' (default), 'ncpsolve' or 'path'
%    eqsolveroptions  : options structure to be passed to eqsolver
%    extrapolate      : 1 if extrapolation is allowed outside the
%                       interpolation space or 0 to forbid it (default: 1)
%    functional       : 1 if the equilibrium equations are a functional equation
%                       problem (default: 0)
%    loop_over_s      : 0 (default) to solve all grid points at once, 1 to loop
%                       over each grid points, or n to loop over n blocks of
%                       grid points
%    funapprox        : 'expapprox', 'expfunapprox', or 'resapprox' (default)
%    simulmethod      : 'interpolation' (default) or 'solve'
%    stat             : 1 to ouput summary statistics from the simulation
%                       (default: 0)
%
% [SSIM,XSIM] = RECSSIMUL(MODEL,INTERP,S0,NPER,...) returns the nrep-by-m-by-nper
% array XSIM that contains the simulated response variables.
%
% [SSIM,XSIM,ESIM] = RECSSIMUL(MODEL,INTERP,S0,NPER,...) returns the
% nrep-by-q-by-nper array ESIM that contains the shocks.
%
% [SSIM,XSIM,ESIM,STAT] = RECSSIMUL(MODEL,INTERP,S0,NPER,...) returns summary
% statistics as a structure that contains the number of observations (STAT.N),
% the moments (STAT.MOMENTS), the correlation between variables (STAT.COR), and
% the autocorrelation (STAT.ACOR). Asking RECSSIMUL to return STAT forces the
% OPTIONS.STAT to 1.
%
% [SSIM,XSIM,ESIM,STAT,FSIM] = RECSSIMUL(MODEL,INTERP,S0,NPER,...)  returns the
% nrep-by-m-by-nper array FSIM that contains the value of the equilibrium
% equations on the simulation.
%
% See also RECSACCURACY, RECSDECISIONRULES, RECSIRF.

% Copyright (C) 2011-2013 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
if nargin<5;
  shocks = [];
  if nargin<4
    nper = [];
    if nargin<3
      error('Not enough input arguments');
    end
  end
end

[d,m] = model.dim{1:2};
validateattributes(s0,{'numeric'},{'size',[NaN,d],'nonempty'},3)
nrep  = size(s0,1);

if ~isempty(shocks) && (numel(shocks)~=d || isempty(nper))
  nper = size(shocks,3)+1;
end

if isempty(nper) || nper==0, nper = 1; end

defaultopt = struct(...
    'accuracy'        , 0                                  ,...
    'display'         , 1                                  ,...
    'eqsolver'        , 'lmmcp'                            ,...
    'eqsolveroptions' , struct('Diagnostics'    , 'off'    ,...
                               'DerivativeCheck', 'off'    ,...
                               'Jacobian'       , 'on')    ,...
    'extrapolate'     , 1                                  ,...
    'funapprox'       , 'resapprox'                        ,...
    'functional'      , 0                                  ,...
    'loop_over_s'     , 0                                  ,...
    'simulmethod'     , 'interpolation'                    ,...
    'stat'            , 0);
if nargin<6
  options = defaultopt;
else
  warning('off','catstruct:DuplicatesFound')
  if isfield(options,'eqsolveroptions')
    options.eqsolveroptions = catstruct(defaultopt.eqsolveroptions,options.eqsolveroptions);
  end
  options = catstruct(defaultopt,options);
end

display         = options.display;
extrapolate     = options.extrapolate;
functional      = options.functional;
funapprox       = lower(options.funapprox);
simulmethod     = lower(options.simulmethod);
statdisplay     = options.stat;
if ~any(strcmp(simulmethod,{'interpolation','solve'}))
  warning('RECS:OptionError',['The simulmethod field can take only the values ' ...
                      '''interpolation'' or ''solve''. Simulations will '        ...
                      'be carried out using the default option, ''interpolation''.'])
end

b         = model.b;
e         = model.e;
f         = model.f;
g         = model.g;
h         = model.h;
ixforward = model.ixforward;
params    = model.params;
w         = model.w;
if isfield(model,'funrand') || (isobject(model) && ~isempty(model.funrand)) % Check if a random shocks generator function is provided
  funrand = model.funrand;
else % Use the discretisation to generate the shocks
  funrand = @(N) e(discrand(N,w),:); % could be implemented also with datasample
end
q        = size(funrand(1),2);

fspace  = interp.fspace;
if isfield(interp,'cX')
  cX = interp.cX;
  T  = size(interp.cX,3);
else
  cz      = interp.cz;
  cx      = interp.cx;
  if functional, params = [params fspace cx]; end
end

%% Generate shocks
output      = struct('F',1,'Js',0,'Jx',0,'Jz',0);
ssim        = zeros(nrep,d,nper);
ssim(:,:,1) = s0;
xsim        = zeros(nrep,m,nper);
esim        =   NaN(nrep,q,nper);
if nargout==5, fsim = zeros(nrep,m,nper); end
if isempty(shocks)
  for t=2:nper, esim(:,:,t) = funrand(nrep); end
elseif numel(shocks)==d
  esim(:,:,2:end) = shocks(ones(nrep,1),:,ones(nper,1));
else
  esim(:,:,2:end) = shocks;
end

%% Simulate the model
for t=1:nper
  if t>1, ssim(:,:,t) = g(ssim(:,:,t-1),xsim(:,:,t-1),esim(:,:,t),params,output); end
  if extrapolate, sinterp = ssim(:,:,t);
  else
    sinterp = max(min(ssim(:,:,t),fspace.b(ones(nrep,1),:)), ...
                  fspace.a(ones(nrep,1),:));
  end
  [LB,UB]    = b(ssim(:,:,t),params);
  Phi        = funbasx(fspace,sinterp);
  if exist('cX','var')
    xsim(:,:,t)       = min(max(funeval(cX(:,:,min(t,T)),fspace,Phi),LB),UB);
    if nargout==5, fval = NaN(nrep,m); end
  else
    xsim(:,:,t)       = min(max(funeval(cx,fspace,Phi),LB),UB);
    switch simulmethod
      case 'solve'
        switch funapprox
          case 'expapprox'
            [xsim(:,:,t),fval] = recsSolveEquilibrium(ssim(:,:,t),...
                                                      xsim(:,:,t),...
                                                      funeval(cz,fspace,Phi),...
                                                      b,f,g,h,...
                                                      params,...
                                                      [],[],[],[],...
                                                      ixforward,options);
          case 'resapprox'
            [xsim(:,:,t),fval] = recsSolveEquilibrium(ssim(:,:,t),...
                                                      xsim(:,:,t),...
                                                      zeros(nrep,0),...
                                                      b,f,g,h,...
                                                      params,...
                                                      cx(:,ixforward),e,w,fspace,...
                                                      ixforward,options);
          case 'expfunapprox'
            if functional, params{end} = interp.ch; end
            [xsim(:,:,t),fval] = recsSolveEquilibrium(ssim(:,:,t),...
                                                      xsim(:,:,t),...
                                                      zeros(nrep,0),...
                                                      b,f,g,h,...
                                                      params,...
                                                      interp.ch,e,w,fspace,...
                                                      ixforward,options);
        end % switch funapprox
      otherwise
        if nargout==5
          fval = f(ssim(:,:,t),xsim(:,:,t),funeval(cz,fspace,Phi),params,output);
        end
    end % switch simulmethod
  end % Finite or infinite horizon problem
  if nargout==5, fsim(:,:,t) = fval; end
end % for t

if exist('cX','var') && t>T
  warning('RECS:ExceedHorizon','Simulate beyond last period')
end

%% Check if state satisfies bounds
ssimlong = reshape(permute(ssim,[1 3 2]),[],d);
if any(min(ssimlong,[],1)<fspace.a)
  warning('RECS:Extrapolation','Extrapolating beyond smin')
end
if any(max(ssimlong,[],1)>fspace.b)
  warning('RECS:Extrapolation','Extrapolating beyond smax')
end
clear('ssimlong')

%% Compute some descriptive statistics
if (nargout>=4 || statdisplay) && (nper >= 40)
  X = cat(2,ssim,xsim);
  X = permute(X(:,:,21:end),[2 1 3]);
  X = reshape(X,d+m,[])';

  % Sample size
  stat.n = size(X,1);

  % Percent of time spent at the bounds
  [LB,UB] = b(X(:,1:d),params);
  pLB     = [NaN(1,d) mean(abs(X(:,d+1:d+m)-LB)<eps,1)*100];
  pUB     = [NaN(1,d) mean(abs(UB-X(:,d+1:d+m))<eps,1)*100];

  mX   = mean(X,1);
  y    = bsxfun(@minus,X,mX);
  varX = mean(y.*y,1);
  stat.moments = [mX' sqrt(varX)' (mean(y.^3,1)./varX.^1.5)' ...
                  (mean(y.^4,1)./(varX.*varX))' min(X)' max(X)' pLB' pUB'];
  if display==1
    disp('Statistics from simulated variables (excluding the first 20 observations):');
    disp(' Moments');
    disp('    Mean      Std. Dev. Skewness  Kurtosis  Min       Max       %LB       %UB');
    disp(stat.moments(1:d,1:end-2))
    disp(stat.moments(d+1:end,:))
  end

  stat.cor = corrcoef(X);
  if display==1
    disp(' Correlation');
    disp(stat.cor);

    figure
    for i=1:d+m
      subplot(ceil((d+m)/ceil(sqrt(d+m))),ceil(sqrt(d+m)),i)
      hist(X(:,i),log2(size(X,1))+1)
    end
  end

  X = cat(2,ssim,xsim);
  X = permute(X(:,:,21:end),[3 2 1]);
  acor = zeros(d+m,5);
  parfor n=1:nrep
    acor = acor+autocor(X(:,:,n))/nrep;
  end
  stat.acor = acor;
  if display==1
    disp(' Autocorrelation');
    disp('    1         2         3         4         5');
    disp(stat.acor(1:d,:));
    disp(stat.acor(d+1:end,:));
  end
end % stat

%% Check accuracy
if options.accuracy
  recsAccuracy(model,interp,ssim(:,:,21:end),options);
end


function acor = autocor(x)
%% AUTOCOR Computes the sample autocorrelation of a matrix of time-series x

maxlag = 5;
[N,d]  = size(x);
acov   = zeros(maxlag+1,d);
for n=0:maxlag
  acov(n+1,:) = mean(x(1+n:N,:).*x(1:N-n,:),1)-mean(x(1+n:N,:),1).*mean(x(1:N-n,:),1);
end

acor                  = (acov(2:maxlag+1,:)./acov(ones(maxlag,1),:))';

% If variance is smaller than precision, correlation is equal to 1
testprecision         = eps(mean(x))>acov(1,:);
acor(testprecision,:) = 1;

return