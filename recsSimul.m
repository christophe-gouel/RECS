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
%    functional       : 1 if the equilibrium equations are a functional equation
%                       problem (default: 0)
%    loop_over_s      : 0 (default) to solve all grid points at once or 1 to loop
%                       over each grid points
%    method           : 'expapprox' (default), 'expfunapprox', 'resapprox-simple'
%                       or 'resapprox-complete'
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
% [SSIM,XSIM,ESIM,FSIM,STAT] = RECSSIMUL(MODEL,INTERP,S0,NPER,...) returns summary
% statistics as a structure that contains the moments (STAT.MOMENTS), the
% correlation between variables (STAT.COR), and the autocorrelation
% (STAT.ACOR). Asking RECSSIMUL to return STAT forces the OPTIONS.STAT to 1.
%
% See also RECSACCURACY.

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

if nargin<6
  options = struct([]);
  if nargin<5;
    shocks = [];
    if nargin<4
      error('Nor enough input arguments');
    end
  end
end

if ~isempty(shocks), nper = size(shocks,3); end

defaultopt = struct(...
    'eqsolver'        , 'ncpsolve'     ,...
    'eqsolveroptions' , struct([])     ,...
    'extrapolate'     , 1              ,...
    'functional'      , 0              ,...
    'loop_over_s'     , 0              ,...
    'method'          , 'expapprox'    ,...
    'simulmethod'     , 'interpolation',...
    'stat'            , 0);
warning('off','catstruct:DuplicatesFound')

options         = catstruct(defaultopt,options);

functional      = options.functional;
method          = lower(options.method);
simulmethod     = lower(options.simulmethod);
statdisplay     = options.stat;

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
  funrand = @(N) e(discrand(N,w),:);
end

fspace  = interp.fspace;
cz      = interp.cz;
cx      = interp.cx;
if functional, params = [params fspace cx]; end

[nrep,d] = size(s0);
m        = size(interp.cx,2);
q        = size(funrand(1),2);

output = struct('F',1,'Js',0,'Jx',0,'Jz',0);
ssim   = zeros(nrep,d,nper+1);
xsim   = zeros(nrep,m,nper+1);
esim   = zeros(nrep,q,nper+1);
if nargout>=4, fsim = zeros(nrep,m,nper+1); end
if isempty(shocks)
  for t=2:nper+1, esim(:,:,t) = funrand(nrep); end
else
  esim(:,:,2:end) = shocks;
end

for t=1:nper+1
  if t>1, s0 = func('g',s0,xx,[],esim(:,:,t),[],[],params,output); end
  [LB,UB]    = func('b',s0,[],[],[],[],[],params);
  Phi        = funbasx(fspace,s0);
  xx         = min(max(funeval(cx,fspace,Phi),LB),UB);
  switch simulmethod
   case 'solve'
    switch method
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
   case 'interpolation'
    if nargout>=4, f = func('f',s0,xx,funeval(cz,fspace,Phi),[],[],[],params,output); end
  end
  ssim(:,:,t) = s0;
  xsim(:,:,t) = xx;
  if nargout>=4, fsim(:,:,t) = f; end
end

% Check if state satisfies bounds
snmin = min(reshape(permute(ssim,[1 3 2]),[],d));
snmax = max(reshape(permute(ssim,[1 3 2]),[],d));
if any(snmin<fspace.a), warning('recs:Extrapolation','Extrapolating beyond smin'), end;
if any(snmax>fspace.b), warning('recs:Extrapolation','Extrapolating beyond smax'), end;

% Compute some descriptive statistics
if (nargout==5 || statdisplay) && (nper > 40)
  X = cat(2,ssim,xsim);
  X = permute(X(:,:,20:end),[2 1 3]);
  X = reshape(X,d+m,[])';

  % Percent of time spent at the bounds
  [LB,UB] = func('b',X(:,1:d),[],[],[],[],[],params);
  pLB     = [zeros(1,d) mean(abs(X(:,d+1:d+m)-LB)<eps)*100];
  pUB     = [zeros(1,d) mean(abs(UB-X(:,d+1:d+m))<eps)*100];

  mX   = mean(X);
  y    = bsxfun(@minus,X,mX);
  varX = mean(y.*y);
  stat.moments = [mX' sqrt(varX)' (mean(y.^3)./varX.^1.5)' (mean(y.^4)./(varX.*varX))' ...
                  min(X)'    max(X)' pLB' pUB'];
% $$$   stat.moments = [mean(X,1)' std(X,0,1)' skewness(X)' kurtosis(X,1,1)' ...
% $$$                   min(X)'    max(X)' pLB' pUB'];
  disp('Statistics from simulated variables (excluding the first 20 observations):');
  disp(' Moments');
  disp('    Mean      Std. Dev. Skewness  Kurtosis  Min       Max       %LB       %UB');
  disp(stat.moments)

  stat.cor = corrcoef(X);
  disp(' Correlation');
  disp(stat.cor);

  figure
  for i=1:d+m
    subplot(ceil((d+m)/ceil(sqrt(d+m))),ceil(sqrt(d+m)),i)
    hist(X(:,i),log2(size(X,1))+1)
  end

  X = cat(2,ssim,xsim);
  X = permute(X(:,:,20:end),[3 1 2]);
  X = reshape(X,[],nrep*(d+m));

  stat.acor = autocor(X);
  stat.acor = reshape(stat.acor,nrep,d+m,[]);
  stat.acor = squeeze(mean(stat.acor,1));
  disp(' Autocorrelation');
  disp('    1         2         3         4         5');
  disp(stat.acor);
end


function acor = autocor(x)
% AUTOCOR Computes the sample autocorrelation of a matrix of time-series x

maxlag = 5;
[N,d]  = size(x);
acov   = zeros(maxlag+1,d);
for n=0:maxlag
  acov(n+1,:) = mean(x(1+n:N,:).*x(1:N-n,:))-mean(x(1+n:N,:)).*mean(x(1:N-n,:));
end
acor   = (acov(2:maxlag+1,:)./acov(ones(maxlag,1),:))';

return