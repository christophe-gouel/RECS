function [si,xi] = recsIRF(model,interp,shock,nrep,nper,options)
% RECSIRF Computes Impulse Response Functions
% 
% SI = RECSIRF(MODEL,INTERP,SHOCK,NREP,NPER)
%
% RECSIRF(MODEL,INTERP,SHOCK,NREP,NPER,OPTIONS)
%
% [SI,XI] = RECSIRF(MODEL,INTERP,SHOCK,NREP,NPER,...)
%
% References
% [1] Koop, G., Pesaran, M. H. and Potter, S. M. (1996), Impulse response
% analysis in nonlinear multivariate models, Journal of Econometrics, 74(1), 119-147.
% [2] Coeurdacier, N., Rey, H. and Winant, P. (2011), The Risky Steady State,
% American Economic Review - Papers and Proceedings, 101(3), 398-401.  
% DOI: <a href="http://dx.doi.org/10.1257/aer.101.3.398">10.1257/aer.101.3.398</a>
%
% See also RECSSIMUL.

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
if nargin < 6
  options = struct([]);
  if nargin < 5
    error('Nor enough input arguments');
  end
end

defaultopt = struct(...
    'stat'        , 0);

options    = catstruct(options,defaultopt);

e       = model.e;
funrand = model.funrand;
w       = model.w;

fspace = interp.fspace;

s0    = (fspace.a+fspace.b)/2;
snorm = 1;
it    = 0;

%% Find risky steady state
while(snorm > eps && it < 1000)
  it    = it+1;
  srss  = recsSimul(model,interp,s0,[],w'*e,options);
  snorm = norm(srss(:,:,end)-s0);
  s0    = srss(:,:,end);
end

srss = s0;

%% Calculate IRF
% Calculate reference and shocked simulations
d       = size(s0,2);
q       = size(funrand(1),2);
m       = size(interp.cx,2);
shocks0 = zeros(nrep,q,nper+1);
shocks1 = zeros(nrep,q,nper+1);
for t=1:nper+1
  shocks0(:,:,t) = funrand(nrep);
end

shocks1(:,:,2:end) = shocks0(:,:,2:end);
shocks1(:,:,1)     = repmat(shock,[nrep 1 1]);

[s0,x0] = recsSimul(model,interp,srss(ones(nrep,1),:),[],shocks0,options);
[s1,x1] = recsSimul(model,interp,srss(ones(nrep,1),:),[],shocks1,options);

% Difference between shocked and reference
si = permute(mean(s1(:,:,2:end)-s0(:,:,2:end),1),[3 2 1]);
xi = permute(mean(x1(:,:,2:end)-x0(:,:,2:end),1),[3 2 1]);

%% Plot IRF
figure
it = 0;
for i=1:d
  it = it+1;
  subplot(ceil((d+m)/ceil(sqrt(d+m))),ceil(sqrt(d+m)),it)
  plot(si(:,i))
end
for i=1:m
  it = it+1;
  subplot(ceil((d+m)/ceil(sqrt(d+m))),ceil(sqrt(d+m)),it)
  plot(xi(:,i))
end
