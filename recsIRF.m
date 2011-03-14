function [si,xi] = recsIRF(model,interp,shock,nrep,nper,options)
% RECSIRF Computes Impulse Response Functions

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

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

% Find stochastic steady state
while(snorm > eps && it < 1000)
  it    = it+1;
  ssss  = recsSimul(model,interp,s0,[],w'*e,options);
  snorm = norm(ssss(:,:,end)-s0);
  s0    = ssss(:,:,end);
end

ssss = s0;

d       = size(s0,2);
q       = size(feval(funrand,1),2);
m       = size(interp.cx,2);
shocks0 = zeros(nrep,q,nper+1);
shocks1 = zeros(nrep,q,nper+1);
for t=1:nper+1
  shocks0(:,:,t) = feval(funrand,nrep);
end

shocks1(:,:,2:end) = shocks0(:,:,2:end);
shocks1(:,:,1)     = repmat(shock,[nrep 1 1]);

[s0,x0] = recsSimul(model,interp,ssss(ones(nrep,1),:),[],shocks0,options);
[s1,x1] = recsSimul(model,interp,ssss(ones(nrep,1),:),[],shocks1,options);

si = permute(mean(s1(:,:,2:end)-s0(:,:,2:end),1),[3 2 1]);
xi = permute(mean(x1(:,:,2:end)-x0(:,:,2:end),1),[3 2 1]);

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
