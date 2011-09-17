function [si,xi] = recsIRF(model,interp,shock,nrep,nper,options)
% RECSIRF Computes Impulse Response Functions (IRF)
%
% RECSIRF computes IRF based on the approach of Koop et al. [1].
%
% SI = RECSIRF(MODEL,INTERP,SHOCK,NREP,NPER) computes the IRF for the model defined
% in the structure MODEL, by using the interpolation structure defined in the
% structure INTERP. The IRF are calculated on NREP scenarios of NPER periods
% each. The impulse is given by the 1-by-q vector SHOCK which is applied at the
% initial period while the model is on its risky steady state (for a definition,
% see Coeurdacier et al. [2]). RECSIRF produces plots of the IRF and returns the
% NPER-by-d matrix SI, containing the difference between the reference simulation
% and the shock simulation (e.g., the IRF) for the state variables.
%
% RECSIRF(MODEL,INTERP,SHOCK,NREP,NPER,OPTIONS) computes the IRF with the
% parameters defined by the structure OPTIONS. The fields of the structure are
% those used in recsSimul.
%
% [SI,XI] = RECSIRF(MODEL,INTERP,SHOCK,NREP,NPER,...) returns the NPER-by-m matrix
% XI, containing the IRF for the response variables.
%
% References
% [1] Koop, G., Pesaran, M. H. and Potter, S. M. (1996), Impulse response
% analysis in nonlinear multivariate models, Journal of Econometrics, 74(1), 119-147.
% DOI: <a href="http://dx.doi.org/10.1016/0304-4076(95)01753-4">10.1016/0304-4076(95)01753-4</a>
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
