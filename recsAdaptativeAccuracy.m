function res = recsAdaptativeAccuracy(interp,model,snodes,options)
% RECSADAPTATIVEACCURACY 

% Copyright (C) 2011 Christophe Gouel
% Licensed under the Expat license, see LICENSE.txt

%% Initialization
if nargin <=3, options = struct([]); end

defaultopt = struct('simulmethod','solve',...
                    'stat'       ,0);
warning('off','catstruct:DuplicatesFound')

options = catstruct(options,defaultopt);

res = zeros(size(snodes,2),1);             
for i=1:size(snodes,2)
  MidPoints = blktridiag(1,1,0,length(snodes{i}))/2;
  snodei = (snodes{i}'*MidPoints(:,1:end-1))';
  
  snodesnew    = snodes;
  snodesnew{i} = snodei;
  s            = gridmake(snodesnew);
  [~,x] = recsSimul(model,interp,s,0,[],options);
  Phi   = funbasx(interp.fspace,s);
  R = recsResidual(s,x,model.func,model.params,interp.cx,interp.fspace,...
                   'resapprox-complete',Phi);
  res(i) = norm(R,Inf);
end