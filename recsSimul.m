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
%    UseParallel      : 'always' (default) to use parallel calculation (require
%                       Parallel Computing Toolbox)' or never'
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

% Copyright (C) 2011-2022 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
if nargin<5
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
    'ArrayProblem'    , false                              ,...
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
    'stat'            , 0                                  ,...
    'Tburn'           , 20                                 ,...
    'UseParallel'     , 'always');
if nargin<6
  options = defaultopt;
else
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
Tburn           = options.Tburn;
if ~any(strcmp(simulmethod,{'interpolation','solve'}))
  warning('RECS:OptionError',['The simulmethod field can take only the values ' ...
                      '''interpolation'' or ''solve''. Simulations will '        ...
                      'be carried out using the default option, ''interpolation''.'])
end
switch lower(options.UseParallel)
  case 'never'
    UseParallel = 0;
  case 'always'
    UseParallel = nrep;
end

b         = model.functions.b;
e         = model.shocks.e;
f         = model.functions.f;
g         = model.functions.g;
h         = model.functions.h;
ixforward = model.infos.ixforward;
params    = model.params;
w         = model.shocks.w;
if isfield(model.shocks,'funrand') % Check if a random shocks generator function is provided
  funrand = model.shocks.funrand;
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
ssim        = zeros(nrep,d,nper);
ssim(:,:,1) = s0;
xsim        = zeros(nrep,m,nper);
esim        =   NaN(nrep,q,nper);
if nargout==5, fsim = zeros(nrep,m,nper); end
if isempty(shocks)
  for t=2:nper, esim(:,:,t) = funrand(nrep); end
elseif numel(shocks)==q
  esim(:,:,2:end) = shocks(ones(nrep,1),:,ones(nper-1,1));
else
  esim(:,:,2:end) = shocks;
end

%% Simulate the model
for t=1:nper
  if t>1, ssim(:,:,t) = g(ssim(:,:,t-1),xsim(:,:,t-1),esim(:,:,t),params); end
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
          fval = f(ssim(:,:,t),xsim(:,:,t),funeval(cz,fspace,Phi),params);
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
vari     = 1:d;
varmin   = vari(min(ssimlong,[],1)<fspace.a);
varmax   = vari(max(ssimlong,[],1)>fspace.b);
if ~isempty(varmin)
  warning('RECS:Extrapolation','Extrapolating state variables (%s) beyond smin',...
          int2str(varmin))
end
if ~isempty(varmax)
  warning('RECS:Extrapolation','Extrapolating state variables (%s) beyond smax',...
          int2str(varmax))
end
clear('ssimlong')

%% Compute some descriptive statistics
if nargout>=4 || statdisplay
  if nper >= Tburn+20
    if exist('table','file'), tabularform = true;
    else                      tabularform = false;
    end
    symbols = [model.symbols.states model.symbols.controls];

    X = cat(2,ssim,xsim);
    X = permute(X(:,:,(Tburn+1):end),[2 1 3]);
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

    if tabularform
      stat.moments = array2table(stat.moments,...
                                 'RowNames',symbols,...
                                 'VariableNames',...
                                 {'Mean' 'StdDev' 'Skewness' 'Kurtosis' 'Min' 'Max' 'pLB' 'pUB'});
    end
    if display==1
      fprintf(1,'Statistics from simulated variables (excluding the first %i observations):\n',Tburn);
      disp(' Moments');
      if tabularform
        disp(stat.moments)
      else
        disp('    Mean      Std. Dev. Skewness  Kurtosis  Min       Max       %LB       %UB');
        disp(stat.moments(1:d,1:end-2))
        disp(stat.moments(d+1:end,:))
      end
    end

    stat.cor = corrcoef(X);
    if tabularform
      stat.cor = array2table(stat.cor,'RowNames',symbols,'VariableNames',symbols);
    end
    if display==1
      disp(' Correlation');
      disp(stat.cor);

      figure
      for i=1:d+m
        subplot(ceil((d+m)/ceil(sqrt(d+m))),ceil(sqrt(d+m)),i)
        hist(X(:,i),log2(size(X,1))+1)
        xlabel(symbols{i});
      end
    end

    X = cat(2,ssim,xsim);
    X = permute(X(:,:,(Tburn+1):end),[3 2 1]);
    acor = zeros(d+m,5);
    parfor (n=1:nrep, UseParallel)
      acor = acor+autocor(X(:,:,n))/nrep;
    end
    if tabularform
      stat.acor = array2table(acor,'RowNames',symbols,...
                              'VariableNames',{'T1' 'T2' 'T3' 'T4' 'T5'});
    else
      stat.acor = acor;
    end
    if display==1
      disp(' Autocorrelation');
      if tabularform
        disp(stat.acor);
      else
        disp('    1         2         3         4         5');
        disp(stat.acor(1:d,:));
        disp(stat.acor(d+1:end,:));
      end
    end
  else
    warning('Insufficient number of observations after burn-in period')
    stat = [];
  end
end % stat

%% Check accuracy
if options.accuracy
  recsAccuracy(model,interp,ssim(:,:,(Tburn+1):end),options);
end
