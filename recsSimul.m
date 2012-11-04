function [ssim,xsim,esim,fsim,stat] = recsSimul(model,interp,s0,nper,shocks,options)
% RECSSIMUL Simulates a model from starting values given in s0 and for nper period
%
% SSIM = RECSSIMUL(MODEL,INTERP,S0,NPER) simulates the model defined in the
% structure MODEL, by using the interpolation structure defined in the structure
% INTERP. The simulation starts from the initial state S0 and lasts NPER (scalar)
% periods. S0 is a nrep-by-d matrix with nrep the number of scenarios to simulate,
% and d the number of state variables. RECSSIMUL returns the nrep-by-d-by-nper+1
% array SSIM that contains the simulated state variables.
% MODEL is a structure, which includes the following fields:
%    [e,w]   : discrete distribution with finite support with e the values and w
%              the probabilities (it could be also the discretisation of a
%              continuous distribution through quadrature or Monte Carlo drawings)
%    func    : function name or anonymous function that defines the model's equations
%    funrand : random shocks generator function (optional). Function handle that
%              accepts an integer n as input and returns a n-by-q matrix of
%              shocks. If not provided, shocks will be drawn from the discrete
%              distribution [e,w].
%    params  : model's parameters, it is preferable to pass them as a cell array
%              (compulsory with the functional option) but other formats are
%              acceptable
% INTERP is a structure, which includes the following fields:
%    ch      : coefficient matrix of the interpolation of the expectations
%              function (optional, to be provided with method 'expfunapprox')
%    cx      : coefficient matrix of the interpolation of the response variables
%    cz      : coefficient matrix of the interpolation of the expectations variables
%    fspace  : a definition structure for the interpolation family (created by
%              the function fundef)
%
% SSIM = RECSSIMUL(MODEL,INTERP,S0,NPER,SHOCKS) uses the nrep-by-q-by-nper array
% SHOCKS to simulate the model instead of drawing random numbers. In this case
% size(SHOCKS,3) supersedes NPER.
%
% SSIM = RECSSIMUL(MODEL,INTERP,S0,NPER,SHOCKS,OPTIONS) simulates the model with
% the parameters defined by the structure OPTIONS. The fields of the structure are
%    eqsolver         : 'fsolve', 'lmmcp', 'ncpsolve' (default) or 'path'
%    eqsolveroptions  : options structure to be passed to eqsolver
%    extrapolate      : 1 if extrapolation is allowed outside the
%                       interpolation space or 0 to forbid it (default: 1)
%    functional       : 1 if the equilibrium equations are a functional equation
%                       problem (default: 0)
%    loop_over_s      : 0 (default) to solve all grid points at once or 1 to loop
%                       over each grid points
%    funapprox        : 'expapprox', 'expfunapprox', 'resapprox-simple'
%                       or 'resapprox-complete' (default)
%    simulmethod      : 'interpolation' (default) or 'solve'
%    stat             : 1 to ouput summary statistics from the simulation
%                       (default: 0)
%
% [SSIM,XSIM] = RECSSIMUL(MODEL,INTERP,S0,NPER,...) returns the nrep-by-m-by-nper+1
% array XSIM that contains the simulated response variables.
%
% [SSIM,XSIM,ESIM] = RECSSIMUL(MODEL,INTERP,S0,NPER,...) returns the
% nrep-by-q-by-nper+1 array ESIM that contains the shocks.
%
% [SSIM,XSIM,ESIM,FSIM] = RECSSIMUL(MODEL,INTERP,S0,NPER,...) returns the
% nrep-by-m-by-nper+1 array FSIM that contains the value of the equilibrium
% equations on the simulation.
%
% [SSIM,XSIM,ESIM,FSIM,STAT] = RECSSIMUL(MODEL,INTERP,S0,NPER,...) returns
% summary statistics as a structure that contains the number of observations
% (STAT.N), the moments (STAT.MOMENTS), the correlation between variables
% (STAT.COR), and the autocorrelation (STAT.ACOR). Asking RECSSIMUL to return
% STAT forces the OPTIONS.STAT to 1.
%
% See also RECSACCURACY.

% Copyright (C) 2011-2012 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
if nargin<5;
  shocks = [];
  if nargin<4
    error('Nor enough input arguments');
  end
end

[nrep,d] = size(s0);

if ~isempty(shocks) && numel(shocks)~=d, nper = size(shocks,3); end

defaultopt = struct(...
    'eqsolver'        , 'lmmcp'                            ,...
    'eqsolveroptions'   , struct('DerivativeCheck', 'off' ,...
                                 'Jacobian'       , 'on')  ,...
    'extrapolate'     , 1                                  ,...
    'funapprox'       , 'resapprox-complete'               ,...
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

extrapolate     = options.extrapolate;
functional      = options.functional;
funapprox       = lower(options.funapprox);
simulmethod     = lower(options.simulmethod);
statdisplay     = options.stat;
if ~strcmp(simulmethod,'interpolation') && ~strcmp(simulmethod,'solve')
  warning('RECS:OptionError',['The simulmethod field can take only the values ' ...
                      '''interpolation'' or ''solve''. Simulations will '        ...
                      'be carried out using the default option, ''interpolation''.'])
end

e       = model.e;
params  = model.params;
w       = model.w;
if isa(model.func,'char')
  func = str2func(model.func);
elseif isa(model.func,'function_handle')
  func   = model.func;
else
  error('model.func must be either a string or a function handle')
end
if isfield(model,'funrand') % Check if a random shocks generator function is provided
  funrand = model.funrand;
else % Use the discretisation to generate the shocks
  funrand = @(N) e(discrand(N,w),:); % could be implemented also with datasample
end
q        = size(funrand(1),2);

fspace  = interp.fspace;
if isfield(interp,'cX')
  cX = interp.cX;
  [~,m,T]  = size(interp.cX);
else
  cz      = interp.cz;
  cx      = interp.cx;
  if functional, params = [params fspace cx]; end
  m       = size(interp.cx,2);
end

%% Generate shocks
output = struct('F',1,'Js',0,'Jx',0,'Jz',0);
ssim   = zeros(nrep,d,nper+1);
xsim   = zeros(nrep,m,nper+1);
esim   = zeros(nrep,q,nper+1);
if nargout>=4, fsim = zeros(nrep,m,nper+1); end
if isempty(shocks)
  for t=2:nper+1, esim(:,:,t) = funrand(nrep); end
elseif numel(shocks)==d
  esim(:,:,2:end) = shocks(ones(nrep,1),:,ones(nper,1));
else
  esim(:,:,2:end) = shocks;
end

%% Simulate the model
for t=1:nper+1
  if t>1, s0 = func('g',s0,xx,[],esim(:,:,t),[],[],params,output); end
  if extrapolate, sinterp = s0;
  else sinterp = max(min(s0,fspace.b(ones(nrep,1),:)),fspace.a(ones(nrep,1),:)); end
  [LB,UB]    = func('b',sinterp,[],[],[],[],[],params);
  Phi        = funbasx(fspace,sinterp);
  if exist('cX','var')
    xx       = min(max(funeval(cX(:,:,min(t,T)),fspace,Phi),LB),UB);
    if nargout>=4, f = NaN(nrep,m); end
  else
    xx       = min(max(funeval(cx,fspace,Phi),LB),UB);
    switch simulmethod
      case 'solve'
        switch funapprox
          case {'expapprox','resapprox-simple'}
            [xx,f] = recsSolveEquilibrium(s0,...
                                          xx,...
                                          funeval(cz,fspace,Phi),...
                                          func,...
                                          params,...
                                          [],[],[],[],options);
          case 'resapprox-complete'
            [xx,f] = recsSolveEquilibrium(s0,...
                                          xx,...
                                          zeros(nrep,0),...
                                          func,...
                                          params,...
                                          cx,e,w,fspace,options);
          case 'expfunapprox'
            if functional, params{end} = interp.ch; end
            [xx,f] = recsSolveEquilibrium(s0,...
                                          xx,...
                                          zeros(nrep,0),...
                                          func,...
                                          params,...
                                          interp.ch,e,w,fspace,options);
        end
      otherwise
        if nargout>=4
          f = func('f',s0,xx,funeval(cz,fspace,Phi),[],[],[],params,output);
        end
    end
  end
  ssim(:,:,t) = s0;
  xsim(:,:,t) = xx;
  if nargout>=4, fsim(:,:,t) = f; end
end

if exist('cX','var') && t>T
  warning('RECS:ExceedHorizon','Simulate beyond last period')
end

%% Check if state satisfies bounds
snmin = min(reshape(permute(ssim,[1 3 2]),[],d));
snmax = max(reshape(permute(ssim,[1 3 2]),[],d));
if any(snmin<fspace.a)
  warning('RECS:Extrapolation','Extrapolating beyond smin')
end
if any(snmax>fspace.b)
  warning('RECS:Extrapolation','Extrapolating beyond smax')
end

%% Compute some descriptive statistics
if (nargout==5 || statdisplay) && (nper > 40)
  X = cat(2,ssim,xsim);
  X = permute(X(:,:,20:end),[2 1 3]);
  X = reshape(X,d+m,[])';
  
  % Sample size
  stat.n = size(X,2);
  
  % Percent of time spent at the bounds
  [LB,UB] = func('b',X(:,1:d),[],[],[],[],[],params);
  pLB     = [zeros(1,d) mean(abs(X(:,d+1:d+m)-LB)<eps,1)*100];
  pUB     = [zeros(1,d) mean(abs(UB-X(:,d+1:d+m))<eps,1)*100];

  mX   = mean(X,1);
  y    = bsxfun(@minus,X,mX);
  varX = mean(y.*y,1);
  stat.moments = [mX' sqrt(varX)' (mean(y.^3,1)./varX.^1.5)' ...
                  (mean(y.^4,1)./(varX.*varX))' min(X)' max(X)' pLB' pUB'];
  disp('Statistics from simulated variables (excluding the first 20 observations):');
  disp(' Moments');
  disp('    Mean      Std. Dev. Skewness  Kurtosis  Min       Max       %LB       %UB');
  disp(stat.moments(1:d,:))
  disp(stat.moments(d+1:end,:))

  stat.cor = corrcoef(X);
  disp(' Correlation');
  disp(stat.cor);

  figure
  for i=1:d+m
    subplot(ceil((d+m)/ceil(sqrt(d+m))),ceil(sqrt(d+m)),i)
    hist(X(:,i),log2(size(X,1))+1)
  end

  X = cat(2,ssim,xsim);
  X = permute(X(:,:,20:end),[3 2 1]);
  stat.acor = zeros(d+m,5);
  for n=1:nrep
    stat.acor = stat.acor+autocor(X(:,:,n))/nrep;
  end
  disp(' Autocorrelation');
  disp('    1         2         3         4         5');
  disp(stat.acor(1:d,:));
  disp(stat.acor(d+1:end,:));
end


function acor = autocor(x)
%% AUTOCOR Computes the sample autocorrelation of a matrix of time-series x

maxlag = 5;
[N,d]  = size(x);
acov   = zeros(maxlag+1,d);
for n=0:maxlag
  acov(n+1,:) = mean(x(1+n:N,:).*x(1:N-n,:),1)-mean(x(1+n:N,:),1).*mean(x(1:N-n,:),1);
end
acor   = (acov(2:maxlag+1,:)./acov(ones(maxlag,1),:))';

return