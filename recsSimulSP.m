function [ssim,xsim,esim,stat,fsim] = recsSimulSP(model,interp,s0,nper,esim,options)
% RECSSIMULSP Simulates a model from starting values given in s0 and for nper period
%
% SSIM = RECSSIMULSP(MODEL,INTERP,S0,NPER) simulates the model defined in the
% object MODEL, by using the interpolation structure defined in the structure
% INTERP. The simulation starts from the initial state S0 and lasts NPER (scalar)
% periods. S0 is a nrep-by-d matrix with nrep the number of scenarios to simulate,
% and d the number of state variables. RECSSIMULSP returns the nrep-by-d-by-nper
% array SSIM that contains the simulated state variables.
% MODEL is an object created by recsmodelsp.
% INTERP is a structure, which includes the following fields:
%    ch      : coefficient matrix of the interpolation of the expectations
%              function (optional, to be provided with method 'expfunapprox')
%    cx      : coefficient matrix of the interpolation of the response variables
%    cz      : coefficient matrix of the interpolation of the expectations variables
%    fspace  : a definition structure for the interpolation family (created by
%              the function fundef)
%
% SSIM = RECSSIMULSP(MODEL,INTERP,S0,NPER,OPTIONS) simulates the model with
% the parameters defined by the structure OPTIONS. The fields of the structure are
%    accuracy         : 1 to check accuracy on the asymptotic distribution (default: 0)
%    display          : 1 (default) to display outputs
%    eqsolver         : 'fsolve', 'lmmcp' (default), 'ncpsolve' or 'path'
%    eqsolveroptions  : options structure to be passed to eqsolver
%    extrapolate      : 1 if extrapolation is allowed outside the
%                       interpolation space or 0 to forbid it (default: 1)
%    loop_over_s      : 0 (default) to solve all grid points at once, 1 to loop
%                       over each grid points, or n to loop over n blocks of
%                       grid points
%    simulmethod      : 'interpolation' (default) or 'solve'
%    stat             : 1 to ouput summary statistics from the simulation
%                       (default: 0)
%
% [SSIM,XSIM] = RECSSIMULSP(MODEL,INTERP,S0,NPER,...) returns the nrep-by-m-by-nper
% array XSIM that contains the simulated response variables.
%
% [SSIM,XSIM,ESIM] = RECSSIMULSP(MODEL,INTERP,S0,NPER,...) returns the
% nrep-by-q-by-nper array ESIM that contains the shocks.
%
% [SSIM,XSIM,ESIM,STAT] = RECSSIMULSP(MODEL,INTERP,S0,NPER,...) returns summary
% statistics as a structure that contains the number of observations (STAT.N),
% the moments (STAT.MOMENTS), the correlation between variables (STAT.COR), and
% the autocorrelation (STAT.ACOR). Asking RECSSIMULSP to return STAT forces the
% OPTIONS.STAT to 1.
%
% [SSIM,XSIM,ESIM,STAT,FSIM] = RECSSIMULSP(MODEL,INTERP,S0,NPER,...)  returns the
% nrep-by-m-by-nper array FSIM that contains the value of the equilibrium
% equations on the simulation.
%
% See also RECSACCURACY, RECSDECISIONRULES, RECSIRF.

% Copyright (C) 2011-2018 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
if nargin<5
  esim = [];
  if nargin<4
    nper = [];
    if nargin<3
      error('Not enough input arguments');
    end
  end
end

validateattributes(s0,{'numeric'},{'size',[NaN,model.dim{1,1}],'nonempty'},3)
nrep  = size(s0,1);

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
simulmethod     = lower(options.simulmethod);
statdisplay     = options.stat;
Tburn           = options.Tburn;
if ~any(strcmp(simulmethod,{'interpolation','solve'}))
  warning('RECS:OptionError',['The simulmethod field can take only the values ' ...
                      '''interpolation'' or ''solve''. Simulations will '        ...
                      'be carried out using the default option, ''interpolation''.'])
end

% Extract fields of model
nperiods  = model.nperiods;
params    = model.params;
shocks    = model.shocks;
dim       = model.dim;
functions = model.functions;
ixforward = cell(nperiods,1);
for i=1:nperiods, ixforward{i} = model.infos(i).ixforward; end

D = sum(cell2mat(dim(:,1)));
M = sum(cell2mat(dim(:,2)));

fspace  = interp.fspace;
cX      = interp.cX;

inext    = @(iperiod) (iperiod+1)*(iperiod<nperiods)+1*(iperiod==nperiods);

%% Generate shocks
ssim = cell(nperiods,1); xsim = ssim;
for i=1:nperiods
  ssim{i} = zeros(nrep,dim{i,1},nper); 
  xsim{i} = zeros(nrep,dim{i,2},nper);
end
if isempty(esim)
  esim = cell(nperiods,1);
  for i=1:nperiods
    esim{i} =   NaN(nrep,dim{i,4},nper);
    for t=1:nper, esim{i}(:,:,t) = shocks{i}.funrand(nrep); end
  end
end

if nargout==5
  fsim = cell(nperiods,1);
  for i=1:nperiods, fsim{i} = zeros(nrep,dim{i,2},nper); end
end

%% Simulate the model
for t=1:nper
  for i=1:nperiods
    if i~=1
      ssim{i}(:,:,t) = functions(i-1).g(ssim{i-1}(:,:,t),xsim{i-1}(:,:,t),...
                                        esim{i}(:,:,t),params); 
    elseif t>1
      ssim{1}(:,:,t) = functions(nperiods).g(ssim{nperiods}(:,:,t-1),...
                                             xsim{nperiods}(:,:,t-1),...
                                             esim{1}(:,:,t),params); 
    else
      ssim{1}(:,:,1) = s0;
    end
    if extrapolate, sinterp = ssim{i}(:,:,t);
    else
      sinterp = max(min(ssim{i}(:,:,t),fspace{i}.b(ones(nrep,1),:)),...
                    fspace{i}.a(ones(nrep,1),:));
    end
    [LB,UB]    = functions(i).b(ssim{i}(:,:,t),params);
    Phi        = funbasx(fspace{i},sinterp);
    xsim{i}(:,:,t) = min(max(funeval(cX{i},fspace{i},Phi),LB),UB);
    switch simulmethod
      case 'solve'
        [xsim{i}(:,:,t),fval] = recsSolveEquilibrium(ssim{i}(:,:,t),...
                                                     xsim{i}(:,:,t),...
                                                     zeros(nrep,0),...
                                                     functions(inext(i)).b,...
                                                     functions(i).f,...
                                                     functions(i).g,...
                                                     functions(i).h,...
                                                     params,...
                                                     cX{inext(i)}(:,ixforward{i}),...
                                                     shocks{inext(i)}.e,shocks{inext(i)}.w,...
                                                     fspace{inext(i)},...
                                                     ixforward{i},options,LB,UB);
      otherwise
        if nargout==5
          fval = NaN(nrep,dim{i,2});
        end
    end % switch simulmethod
    if nargout==5, fsim{i}(:,:,t) = fval; end
  end % for i=1:nperiods
end % for t

%% Check if state satisfies bounds
for i=1:nperiods
  ssimlong = reshape(permute(ssim{i},[1 3 2]),[],dim{i,1});
  vari     = 1:dim{i,1};
  varmin   = vari(min(ssimlong,[],1)<fspace{i}.a);
  varmax   = vari(max(ssimlong,[],1)>fspace{i}.b);
  if ~isempty(varmin)
    warning('RECS:Extrapolation','Extrapolating state variables (%s) in subperiod %i beyond smin',...
            int2str(varmin),i)
  end
  if ~isempty(varmax)
    warning('RECS:Extrapolation','Extrapolating state variables (%s) in subperiod %i beyond smax',...
            int2str(varmax),i)
  end
end
clear('ssimlong')

%% Compute some descriptive statistics
if nargout>=4 || statdisplay
  if exist('table','file'), tabularform = true;
  else                      tabularform = false;
  end
  if nper >= Tburn+20
    X = cat(2,ssim{:},xsim{:});
    X = permute(X(:,:,(Tburn+1):end),[2 1 3]);
    X = reshape(X,D+M,[])';
    symbols = struct2cell(model.symbols);
    symbols = cat(2,symbols{1,:,:},symbols{2,:,:});

    % Sample size
    stat.n = size(X,1);
    
    % Percent of time spent at the bounds
    LB = cell(nperiods,1);
    UB = cell(nperiods,1);
    is = 1;
    for i=1:nperiods
      [LB{i},UB{i}] = functions(i).b(X(:,is:(is-1+dim{i,1})),params);
      is = is+dim{i,1};
    end
    LB = cat(2,LB{:});
    UB = cat(2,UB{:});
    pLB = [NaN(1,D) mean(abs(X(:,D+1:D+M)-LB)<eps,1)*100];
    pUB = [NaN(1,D) mean(abs(UB-X(:,D+1:D+M))<eps,1)*100];
    
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
      if ~tabularform
        disp('    Mean      Std. Dev. Skewness  Kurtosis  Min       Max       %LB       %UB');
      end
      disp(stat.moments)
    end
    
    stat.cor = corrcoef(X);
    if tabularform
      stat.cor = array2table(stat.cor,'RowNames',symbols,'VariableNames',symbols);
    end
    if display==1
      disp(' Correlation');
      disp(stat.cor);
      
      figure
      for i=1:D+M
        subplot(ceil((D+M)/ceil(sqrt(D+M))),ceil(sqrt(D+M)),i)
        hist(X(:,i),log2(size(X,1))+1)
        xlabel(symbols{i});
      end
    end
    
    X = cat(2,ssim{:},xsim{:});
    X = permute(X(:,:,(Tburn+1):end),[3 2 1]);
    acor = zeros(D+M,5);
    parfor n=1:nrep
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
      if ~tabularform
        disp('    1         2         3         4         5');
      end
      disp(stat.acor);
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
