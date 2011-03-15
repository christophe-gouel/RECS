function [ssim,xsim,esim,fsim,stat] = recsSimul(model,interp,s0,nper,shocks,options)
% RECSSIMUL Simulates a model from starting values given in s0 and for nper period

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
    'loop_over_s'     , 0              ,...
    'method'          , 'expapprox'    ,...
    'simulmethod'     , 'interpolation',...
    'stat'            , 0              ,...
    'functional'      , 0);
warning('off','catstruct:DuplicatesFound')

options         = catstruct(defaultopt,options);

eqsolver        = options.eqsolver;
eqsolveroptions = options.eqsolveroptions;
loop_over_s     = options.loop_over_s;
method          = lower(options.method);
simulmethod     = lower(options.simulmethod);
statdisplay     = options.stat;

e       = model.e;
params  = model.params;
funrand = model.funrand;
w       = model.w;
if isa(model.func,'char')
  func = str2func(model.func);
elseif isa(model.func,'function_handle')
  func   = model.func;
else
  error('model.func must be either a string or a function handle')
end

fspace  = interp.fspace;
cz      = interp.cz;
cx      = interp.cx;
if options.functional, params = [params fspace cx]; end

[nrep,d] = size(s0);
m        = size(interp.cx,2);
q        = size(funrand(1),2);

ssim = zeros(nrep,d,nper+1);
xsim = zeros(nrep,m,nper+1);
esim = zeros(nrep,q,nper+1);
if nargout>=4, fsim = zeros(nrep,m,nper+1); end
if isempty(shocks)
  for t=2:nper+1, esim(:,:,t) = funrand(nrep); end
else
  esim(:,:,2:end) = shocks;
end

for t=1:nper+1
  if t>1, s0 = func('g',s0,xx,[],esim(:,:,t),[],[],params); end
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
      [xx,f] = recsSolveEquilibrium(s0,...
                                    xx,...
                                    zeros(nrep,0),...
                                    func,...
                                    params,...
                                    interp.ch,e,w,fspace,options);
    end
   case 'interpolation'
    if nargout>=4, f = func('f',s0,xx,funeval(cz,fspace,Phi),[],[],[],params); end
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